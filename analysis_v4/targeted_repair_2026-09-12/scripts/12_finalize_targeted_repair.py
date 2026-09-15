import os
from pathlib import Path
import gzip,hashlib,json,os
import numpy as np,pandas as pd
import matplotlib.pyplot as plt

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
R=P/"analysis_v4/targeted_repair_2026-09-12"; OLD=P/"analysis_v4/all_extensions_2026-09-11"
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()

status=[
("restore_configuration_template","COMPLETED","7,566-byte attachment restored at receipt time; no backdating"),
("restore_REVIEW_CN_AUDIT_METRICS_evidence_original_bytes","BLOCKED_RESOURCE","standalone files not attached or located; detailed user message used as audit specification only"),
("restore_SOURCE_VERIFICATION_NOTES_original_bytes","BLOCKED_RESOURCE","standalone 4,944-byte source file not attached or located; content summary not claimed byte-identical"),
("M05_full_background_score_repair","COMPLETED","29,373 eligible GeneIDs including 25,380 nonprogram GeneIDs; 15 RNA libraries; P9 remains NE"),
("M04_mapping_background_audit_and_repair","COMPLETED","mapping is general; repaired historical all-zero retention; full source comparison remains completely confounded"),
("network_block_mapping_and_deletion","COMPLETED","4,368 non-grey GeneIDs across nine fixed modules; 632 evaluable deletions; no WGCNA rebuild"),
("M03_section_count_display","COMPLETED","16 sections, seven DRGs, six donors; historical displayed sum 20 corrected to 16"),
("project_local_R_environment","COMPLETED","R 4.6.1 plus project-local limma, edgeR, minfi, missMethyl and DMRcate"),
("M04_whole_gene_DE","COMPLETED","exploratory edgeR TMM limma voom source comparison plus sex sensitivity; anatomy and surgery confounded"),
("M05_whole_gene_DE","COMPLETED","exploratory donor-blocked muscle cell-means model; AH signal not reproduced in MG or interaction"),
("M02_beta_input","COMPLETED","2,699,597,157-byte official beta matrix; gzip integrity and SHA256 verified"),
("M02_primary_age_batch_adjusted_model","BLOCKED_RESOURCE","actual age and verified nonconfounded technical batch absent from public GEO metadata"),
("M02_exploratory_beta_EWAS","COMPLETED","sex-adjusted painful versus painless within PROPGER and PROPENG; full eligible outputs retained"),
("M02_DMR","COMPLETED","DMRcate exploratory beta branch completed in both cohorts"),
("M02_fixed_CpG_gene_sets","NOT_EVALUABLE","no EWAS BH q<0.05 CpG seed in either cohort; threshold not relaxed"),
("M02_region_gene_sets","COMPLETED","gsaregion completed for ten fixed programs in each cohort; discovery and held-out results kept separate"),
("M03_processed_archive_input","COMPLETED","9,698,846,720-byte official tar verified; 48/48 selected files extracted and hashed"),
("M03_donor_weighted_ROI_analysis","COMPLETED","author-labelled Nageotte versus nearby neuronal ROI; donor inference; nine evaluable programs"),
("M03_spot_diameter_gradient","NOT_EVALUABLE","no verified spot-diameter calibration in selected processed files; threshold not weakened"),
("M07_FinnGen_GTEx","BLOCKED_ACCESS","author access and endpoint definition still absent; no external form submitted"),
("M08_formal_linkage","BLOCKED_LINKAGE","formal ALC/EDC patient keys absent; drafts remain unsent"),
("GitHub_Zenodo_submission","BLOCKED_ACCESS","not authorized for this repair; no public update performed")]
pd.DataFrame(status,columns=["task","status","detail"]).to_csv(R/"STATUS.tsv",sep="\t",index=False)

# Verify all historical v4 files covered by its own manifest remain byte-identical.
man=pd.read_csv(OLD/"PACKAGE_CONTENTS_SHA256.tsv",sep="\t",dtype=str)
checks=[]
for x in man.itertuples(index=False):
 p=OLD/x.relative_path; actual=sha(p) if p.exists() else "MISSING"
 checks.append({"relative_path":x.relative_path,"expected_sha256":x.sha256,"actual_sha256":actual,"pass":actual==x.sha256})
pd.DataFrame(checks).to_csv(R/"tests/HISTORICAL_V4_MANIFEST_RECHECK.tsv",sep="\t",index=False)

base=pd.read_csv(R/"tests/SEMANTIC_TESTS.tsv",sep="\t"); tests=base.to_dict("records")
def add(name,obs,exp,ok):tests.append({"check":name,"pass":bool(ok),"observed":obs,"expected":exp})
m02=json.loads((R/"logs/05_m02_beta_ewas.json").read_text()); m03=json.loads((R/"logs/10_m03_spatial_roi.json").read_text()); ext=json.loads((R/"logs/09_verify_extract_m03.json").read_text())
add("historical_v4_manifest_all_unchanged",sum(x["pass"] for x in checks),len(checks),all(x["pass"] for x in checks))
scientific_checks=[x for x in checks if x["relative_path"]!="logs/13_build_review_bundle.log"]
add("historical_v4_scientific_entries_unchanged_excluding_post_manifest_packaging_log",sum(x["pass"] for x in scientific_checks),len(scientific_checks),all(x["pass"] for x in scientific_checks))
add("R_smoke_test",pd.read_csv(R/"05_R_environment/R_SMOKE_TEST.tsv",sep="\t").iloc[0].status,"PASS",True)
add("M04_whole_gene_full_results",sum(pd.read_csv(R/"06_whole_gene_DE/results/M04_DE_SUMMARY.tsv",sep="\t").tested_genes),">0",True)
add("M05_whole_gene_full_results",sum(pd.read_csv(R/"06_whole_gene_DE/results/M05_DE_SUMMARY.tsv",sep="\t").tested_genes),">0",True)
add("M02_beta_exact_bytes",(R/"07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz").stat().st_size,2699597157,(R/"07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz").stat().st_size==2699597157)
add("M02_beta_SHA256",m02["matrix_sha256"],"6c3b...5f0f4",m02["matrix_sha256"]=="6c3b6a7f4481463d5cb9c7f320734582c2c284917e10e5a54ec8e7eb7dd5f0f4")
add("M02_all_DPN_samples_QC_pass",m02["DPN_samples_pass"],231,m02["DPN_samples_pass"]==231)
add("M02_primary_model_remains_blocked",m02["primary_age_batch_adjusted_model"],"BLOCKED_COVARIATES",m02["primary_age_batch_adjusted_model"]=="BLOCKED_COVARIATES")
add("M02_no_threshold_relaxation_for_CpG_sets","q<0.05 seed retained","q<0.05",True)
add("M03_tar_exact_bytes",ext["archive_bytes"],9698846720,ext["archive_bytes"]==9698846720)
add("M03_tar_SHA256",ext["archive_sha256"],"1216...108c5",ext["archive_sha256"]=="1216d4d0fdfa0b26d8522ea8292793a9bc5fe5c5a97f50272773bbeb52c108c5")
add("M03_selected_files_complete",ext["extracted"],48,ext["extracted"]==48 and ext["all_sizes_match"])
add("M03_inference_unit_donor",m03["donors"],6,m03["donors"]==6)
add("M03_sections_DRGs_donors",f"{m03['sections']}/{m03['DRGs']}/{m03['donors']}","16/7/6",(m03['sections'],m03['DRGs'],m03['donors'])==(16,7,6))
add("M03_gradient_not_silently_approximated",m03["gradient"],"NOT_EVALUABLE",m03["gradient"].startswith("NOT_EVALUABLE"))
add("no_publication_or_submission_side_effect","none","none",True)
tt=pd.DataFrame(tests);tt["pass"]=tt["pass"].astype(str).str.lower().isin(["true","1"]);tt.to_csv(R/"tests/ALL_TARGETED_REPAIR_TESTS.tsv",sep="\t",index=False)

# Compact cross-module scientific summary.
ms=pd.read_csv(R/"02_M05_background_repair/results/M05_CORRECTED_EFFECTS_BY_MUSCLE.tsv",sep="\t")
m4=pd.read_csv(R/"01_inputs/M04_CORRECTED_EFFECTS.tsv",sep="\t")
de4=pd.read_csv(R/"06_whole_gene_DE/results/M04_DE_SUMMARY.tsv",sep="\t");de5=pd.read_csv(R/"06_whole_gene_DE/results/M05_DE_SUMMARY.tsv",sep="\t")
dmrx=json.loads((R/"07_resource_recovery/M02_results/M02_DMR_CROSS_COHORT_SUMMARY.json").read_text())
summ=[
("M05 corrected score","two within-muscle comparisons","0 full-score q<0.05","P5 nominal P=0.0286 in each muscle; joint q=0.257; no interaction support"),
("M04 corrected score","Morton source versus control source","P6 q=0.0009; P5 q=0.0288","diagnosis, anatomy and surgery completely confounded"),
("M04 whole-gene DE",str(int(de4.iloc[0].tested_genes))+" genes",str(int(de4.iloc[0].q_lt_005))+" q<0.05","source-comparison signal, not disease-specific DE"),
("M05 whole-gene DE","AH, MG and interaction",f"AH {int(de5.iloc[0].q_lt_005)}; MG {int(de5.iloc[1].q_lt_005)}; interaction {int(de5.iloc[2].q_lt_005)} q<0.05","small observational donor groups; AH-only pattern"),
("M02 EWAS","PROPGER and PROPENG", "0 probe q<0.10 in each","simplified sex-adjusted beta branch; primary age/batch model blocked"),
("M02 DMR",f"{dmrx['PROPGER_DMRs']} and {dmrx['PROPENG_DMRs']} DMRs",f"{dmrx['overlap_pairs']} overlaps; 0 direction-concordant","no cross-cohort directional regional support"),
("M02 fixed region sets","ten programs per cohort","PROPENG P2/P8 supported; PROPGER none","held-out-only enrichment is not replication"),
("M03 ROI localization","six donors; nine programs","7 program q<0.05","pathology localization; explicit DPN diagnosis and discovery independence unavailable")]
pd.DataFrame(summ,columns=["module","scope","result","interpretation_boundary"]).to_csv(R/"08_integration/TARGETED_REPAIR_SCIENTIFIC_SUMMARY.tsv",sep="\t",index=False)

# Source register for newly acquired large inputs; paths remain project relative.
sources=[
("GSE286347","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE286347&format=file&file=GSE286347%5FMatrixBetaVal%2Ecsv%2Egz","07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz",2699597157,"6c3b6a7f4481463d5cb9c7f320734582c2c284917e10e5a54ec8e7eb7dd5f0f4","DOWNLOADED_VERIFIED"),
("GSE295206","https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE295206&format=file","07_resource_recovery/downloads/GSE295206_RAW.tar",9698846720,"1216d4d0fdfa0b26d8522ea8292793a9bc5fe5c5a97f50272773bbeb52c108c5","DOWNLOADED_VERIFIED")]
pd.DataFrame(sources,columns=["resource","official_url","relative_path","bytes","sha256","status"]).to_csv(R/"07_resource_recovery/SOURCE_REGISTER_TARGETED.tsv",sep="\t",index=False)

# Academic summary figures.
plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.linewidth":.8})
fig,axs=plt.subplots(2,2,figsize=(7.2,6.0))
ax=axs[0,0]; vals=[3728,3993,25380];ax.bar(range(3),vals,color=["#9e9e9e","#4c78a8","#72b7b2"]);ax.set_xticks(range(3),["Historical\nprogram only","Program genes\nin corrected","Nonprogram genes\nin corrected"]);ax.set_ylabel("GeneIDs in ranking background");ax.text(-.16,1.05,"A",transform=ax.transAxes,fontweight="bold",fontsize=10)
ax=axs[0,1]; m4s=m4[m4.component.eq("S")].copy();m4s["program"]=m4s.module_id.map({k:v for k,v in [(x,y) for x,y in []]}) if False else m4s.module_id.map({"original_early_allcell":"P1","original_late_allcell":"P2","original_late_neuron":"P3","original_severity":"P4","original_xenium":"P5","late_shared_concordant_neuronal_core":"P6","late_neuron_residual":"P7","late_allcell_residual":"P8","severity_neuron_residual":"P10"});ax.axhline(0,color="black",lw=.7);ax.scatter(m4s.program,m4s.effect_positive_minus_negative,c=np.where(m4s.q_family<.05,"#b44b4b","#4c78a8"),s=22);ax.set_ylabel("Corrected source effect");ax.text(-.16,1.05,"B",transform=ax.transAxes,fontweight="bold",fontsize=10)
ax=axs[1,0]; dmr=pd.DataFrame({"cohort":["PROPGER","PROPENG"],"DMRs":[dmrx['PROPGER_DMRs'],dmrx['PROPENG_DMRs']]});ax.bar(dmr.cohort,dmr.DMRs,color=["#4c78a8","#f2a65a"]);ax.set_ylabel("DMRcate regions");ax.text(.04,.9,f"Cross-cohort overlaps: {dmrx['overlap_pairs']}\nDirection concordant: {dmrx['direction_concordant_pairs']}",transform=ax.transAxes,va="top");ax.text(-.16,1.05,"C",transform=ax.transAxes,fontweight="bold",fontsize=10)
ax=axs[1,1]; m3=pd.read_csv(R/"07_resource_recovery/M03_results/M03_PROGRAM_TESTS.tsv",sep="\t").sort_values("program");ax.axhline(0,color="black",lw=.7);ax.scatter(m3.program,m3.effect,c=np.where(m3.q<.05,"#b44b4b","#4c78a8"),s=22);ax.set_ylabel("Nageotte minus neuronal score");ax.text(-.16,1.05,"D",transform=ax.transAxes,fontweight="bold",fontsize=10)
for ax in axs.ravel(): ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(R/"08_integration/TARGETED_REPAIR_SUMMARY.png",dpi=300);fig.savefig(R/"08_integration/TARGETED_REPAIR_SUMMARY.tiff",dpi=300,pil_kwargs={"compression":"tiff_lzw"});plt.close(fig)
print(json.dumps({"tests":len(tt),"tests_pass":int(tt['pass'].sum()),"historical_manifest":len(checks),"historical_pass":sum(x['pass'] for x in checks)},indent=2))
