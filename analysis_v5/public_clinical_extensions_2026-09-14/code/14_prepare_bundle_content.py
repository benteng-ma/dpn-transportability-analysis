import os
from pathlib import Path
import hashlib, json, shutil
import numpy as np
import pandas as pd
from openpyxl import load_workbook

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; INT=RUN/"06_integration"
def sha(f):
    h=hashlib.sha256()
    with open(f,"rb") as z:
        for b in iter(lambda:z.read(8<<20),b""):h.update(b)
    return h.hexdigest()

# Avoid duplicate formal F2 numbering; retain the assay-scale diagnostic clearly.
old=RUN/"03_GSE14806x/figures/M03_F2_RRBS_effects.png"; new=RUN/"03_GSE14806x/figures/M03_DIAGNOSTIC_BETA_VS_M.png"
if old.exists() and not new.exists(): old.rename(new)

figure_sources={
 "M01_F1_structure.png":["01_GSE302658/results/GSE302658_CLINICAL_ANALYSIS_WIDE.tsv"],
 "M01_F2_baseline.png":["01_GSE302658/results/BASELINE_SYMPTOM_ASSOCIATIONS.tsv"],
 "M01_F3_longitudinal.png":["01_GSE302658/results/LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv"],
 "M01_F4_interaction.png":["01_GSE302658/results/TREATMENT_INTERACTION_OMNIBUS.tsv"],
 "M02_F1_cohorts.png":["02_GSE286347/results/M02_PUBLIC_METADATA_AUDIT.tsv"],
 "M02_F2_MODEL_EFFECT_COMPARISON.png":["02_GSE286347/results/M02_PROPGER_EWAS_SEX_BATCH_ALL.tsv.gz","02_GSE286347/results/M02_PROPENG_EWAS_SEX_BATCH_ALL.tsv.gz"],
 "M02_F3_fixed_regions.png":["02_GSE286347/results/M02_FIXED_REGION_EFFECT_COMPARISON.tsv"],
 "M02_F4_burden.png":["02_GSE286347/results/M02_PROPENG_BURDEN_ASSOCIATIONS.tsv"],
 "M03_F1_groups.png":["03_GSE14806x/results/M03_RRBS_SAMPLE_AUDIT.tsv"],
 "M03_F2_RNA_PROMOTER_EFFECTS.png":["03_GSE14806x/results/M03_CROSSOMIC_GENE_EFFECTS.tsv.gz"],
 "M03_F3_PROGRAM_EVIDENCE_MATRIX.png":["03_GSE14806x/results/M03_PROGRAM_CROSSOMIC_SUMMARY.tsv"],
 "M03_F4_GLOBAL_SAMPLE_QC.png":["03_GSE14806x/results/M03_RRBS_GLOBAL_SAMPLE_QC.tsv"],
 "M03_DIAGNOSTIC_BETA_VS_M.png":["03_GSE14806x/results/M03_RRBS_DIFFERENTIAL_RESULTS.tsv.gz"],
 "M04_F1_PUBLIC_FIELD_AVAILABILITY.png":["04_JCI184075/results/M04_DIRECT_KEY_AUDIT.tsv"],
 "M05_F1_real_coordinates.png":["05_GSE295206/results/M05_NAGEOTTE_BURDEN_BY_SECTION.tsv"],
 "M05_F2_burden.png":["05_GSE295206/results/M05_NAGEOTTE_BURDEN_BY_DONOR.tsv"],
 "M05_F3_program_burden.png":["05_GSE295206/results/M05_NAGEOTTE_BURDEN_PROGRAM_ASSOCIATIONS.tsv"],
 "M05_F4_ARM_SPECIFIC_BURDEN.png":["05_GSE295206/results/M05_NAGEOTTE_BURDEN_ARM_ASSOCIATIONS.tsv"],
 "M05_F5_MARKER_PROXY_ROI_DIFFERENCES.png":["05_GSE295206/results/M05_SPATIAL_MARKER_PROXY_SUMMARY.tsv"],
 "M05_F6_TOPOLOGICAL_RING_TRENDS.png":["05_GSE295206/results/M05_TOPOLOGICAL_RING_ANALYSIS.tsv"]}
rows=[]
for fig,srcs in figure_sources.items():
    found=list(RUN.rglob(fig)); rows.append({"figure":fig,"figure_path":str(found[0].relative_to(RUN)) if found else "MISSING","source_tables":"|".join(srcs),"source_exists":all((RUN/s).exists() for s in srcs),"display_role":"formal" if "DIAGNOSTIC" not in fig else "diagnostic","data_status":"actual public-data result"})
pd.DataFrame(rows).to_csv(INT/"FIGURE_SOURCE_REGISTER.tsv",sep="\t",index=False)

# Only direct deposited keys are recorded. No cross-source linkage is inferred.
cross=[]
d=pd.read_csv(RUN/"01_GSE302658/results/GSE302658_CLINICAL_ANALYSIS_WIDE.tsv",sep="\t")
for x in d.itertuples(): cross.append({"module":"M01","source":"GSE302658","participant_key":x.participant_id,"sample_key":x.gsm,"timepoint":x.visit,"group":x.treatment,"hierarchy":"participant>visit sample","cross_source_link":"NONE"})
d=pd.read_csv(RUN/"02_GSE286347/results/M02_PUBLIC_METADATA_AUDIT.tsv",sep="\t")
for x in d.itertuples(): cross.append({"module":"M02","source":"GSE286347","participant_key":x.geo_accession,"sample_key":x.geo_accession,"timepoint":"cross-sectional","group":getattr(x,"phenotype_geo",None) or getattr(x,"phenotype",None),"hierarchy":"public biospecimen; participant-level study record","cross_source_link":"NONE"})
for module,file,source in [("M03_RNA","M03_RNA_SAMPLE_AUDIT.tsv","GSE148059"),("M03_RRBS","M03_RRBS_SAMPLE_AUDIT.tsv","GSE148060")]:
    d=pd.read_csv(RUN/"03_GSE14806x/results"/file,sep="\t")
    for x in d.itertuples(): cross.append({"module":module,"source":source,"participant_key":x.sample_id,"sample_key":x.sample_id,"timepoint":"cross-sectional","group":getattr(x,"group",None) or getattr(x,"dmfd_class",None),"hierarchy":"deposited sample","cross_source_link":"BLOCKED_NO_VERIFIED_RNA_RRBS_KEY"})
d=pd.read_csv(RUN/"05_GSE295206/results/M05_NAGEOTTE_BURDEN_BY_SECTION.tsv",sep="\t")
for x in d.itertuples(): cross.append({"module":"M05","source":"GSE295206","participant_key":x.donor_id,"sample_key":x.geo_accession,"timepoint":"cross-sectional","group":"diabetes history; DPN not asserted","hierarchy":f"donor>{x.DRG_ID}>section","cross_source_link":"NONE"})
pd.DataFrame(cross).to_csv(INT/"SAMPLE_CROSSWALK.tsv",sep="\t",index=False)

sources=pd.DataFrame([
 ["M01","GSE302658","public processed transcript counts and GEO metadata","REUSED_LOCAL","no FASTQ download","104 participants; 205 visits"],
 ["M02","GSE286347","public beta matrix, GEO metadata and IDAT filename metadata","REUSED_LOCAL","IDAT files not downloaded for cell estimation","315 samples; 231 DPN"],
 ["M03","GSE148059","public processed rlog RNA and metadata","REUSED_LOCAL","no RNA-RRBS pairing","77 RNA samples"],
 ["M03","GSE148060","53 public processed CpG-percent files","REUSED_LOCAL","hg19; no RNA pairing","53 RRBS samples"],
 ["M03","GSE148061","parent-series public lineage metadata","AUDITED","no patient crosswalk inferred","parent accession"],
 ["M04","JCI184075","all locally retained public supplements/source values","REUSED_LOCAL","continuous pathology absent","28 files; 13 XLSX"],
 ["M05","GSE295206","public processed spatial H5 and author ROI annotations","REUSED_LOCAL","controlled sequencing not accessed","16 sections; 7 DRGs; 6 donors"]],columns=["module","accession","public_input","availability","boundary","denominator"])
sources.to_csv(INT/"SOURCE_REGISTER.tsv",sep="\t",index=False)

# Consolidate every recorded module input plus required large/public sources.
audit=[]
for man in RUN.glob("0[1-5]_*/inputs_manifest/INPUTS.tsv"):
    x=pd.read_csv(man,sep="\t")
    for r in x.to_dict("records"): audit.append({"module":man.parts[-3].split('_')[0],"role":"recorded module input","path":r.get("path"),"bytes":r.get("bytes"),"sha256":r.get("sha256")})
extra=[
 ("M02","public beta matrix",P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz"),
 ("M02","official GEO series matrix",P/"analysis_v4/all_extensions_2026-09-11/01_inputs/public_downloads/GSE286347/GSE286347_series_matrix.txt.gz"),
 ("ADMIN","original user CodeKit ZIP",Path(os.environ.get("DPN_CODEKIT_ZIP", "DPN_Public_Clinical_Extensions_Codex_Kit.zip"))),
 ("ADMIN","locked analysis plan",RUN/"00_admin/ANALYSIS_LOCK.json")]
known={str(x["path"]) for x in audit}
for mod,role,f in extra:
    if str(f) not in known: audit.append({"module":mod,"role":role,"path":str(f),"bytes":f.stat().st_size,"sha256":sha(f)})
pd.DataFrame(audit).drop_duplicates("path").to_csv(INT/"INPUT_AUDIT.tsv",sep="\t",index=False)

readme="""# DPN public clinical extensions review bundle

This archive contains a public-data-only, post-correction exploratory extension of the fixed v3.1 DPN programs. Start with `06_integration/CLINICAL_EXTENSION_REVIEW_REPORT.md`, `06_integration/MODULE_STATUS.tsv`, and `06_integration/CLINICAL_EVIDENCE_MATRIX_LONG.tsv`.

The live manuscript was not modified. No author contact, controlled access, inferred patient crosswalk, program redefinition, coverage relaxation, GitHub/Zenodo update, DOI change, or submission occurred.

`MANIFEST_SHA256.tsv` covers every payload file except itself. The outer ZIP has a separate adjacent SHA256 file.
"""
(RUN/"README_REVIEW_BUNDLE.md").write_text(readme,encoding="utf-8")
print(json.dumps({"figure_rows":len(rows),"crosswalk_rows":len(cross),"input_audit_rows":len(pd.DataFrame(audit).drop_duplicates('path'))},indent=2))
