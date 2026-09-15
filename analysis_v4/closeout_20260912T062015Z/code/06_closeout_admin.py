from __future__ import annotations

import os
import hashlib,json,platform,sys
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
C=P/"analysis_v4/closeout_20260912T062015Z"; T=P/"analysis_v4/targeted_repair_2026-09-12"
def sha(p):
 h=hashlib.sha256();
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def w(df,name): df.to_csv(C/"00_admin"/name,sep="\t",index=False,na_rep="NA")

status=pd.DataFrame([
 ["C00","CodeKit verification","NOT_EVALUABLE","Only ten loose files were supplied; seal_bundle.py, tests/, smoke_fixed_region.R, docx_three_line.py and REPORTED_BASELINE.md were absent. R scripts parsed; Python helper operations executed on real adapters."],
 ["C01","M04 feature flow and annotation cache","COMPLETED","58,174 raw rows; 35,749 unique mapped rows; 27,506 eligible rank GeneIDs. Fixed NCBI snapshot reproduces all 38,551 cache pairs; report denominator corrected, no matrix change."],
 ["C02","M02 DMR runtime and generation audit","COMPLETED","Installed R runtime audited. Exact rerun reproduced 167 PROPGER and 660 PROPENG coordinates; maximum numeric difference <6e-16."],
 ["C03","PROPGER fixed regions tested in PROPENG","COMPLETED","Known-results retrospective test of all 167 fixed regions in 92 participants; 167 evaluable, four nominal P<0.05, none BH q<0.10; minimum q=0.3027."],
 ["C04","M02 ancillary audit","PARTIAL","Probe QC, M-value scale, sex-adjusted design, DMR and region-set provenance verified. Age/batch primary model, cell composition, clinical linkage and pre-existing arm/promoter analyses remain unavailable or not executed; no proxy substituted."],
 ["C05","M03 spatial ROI closeout","COMPLETED","16 sections, 7 DRG, 6 donors; nine programs evaluable, seven q=0.04018. Real coordinates, two arms, gene contributions and LODO produced. Gradient NE because no image, scale or boundary was available."],
 ["C06","Fixed-module deletion summary","COMPLETED","All 632 planned, evaluable deletions retained with effects, sign, coverage and remaining arm counts. No WGCNA reconstruction and no favorable-block selection."],
 ["C07","M04 and M05 whole-gene closeout","COMPLETED","M04 primary and sex-adjusted families kept distinct; M05 AH signal not called an interaction. Pain/pathology keys unavailable and remain NE."],
 ["C08","FinnGen and ALC/EDC access boundary","BLOCKED_ACCESS","No access agreement, email, patient key or controlled data were obtained. No form was completed, no message sent and no patient mapping guessed."],
])
status.columns=["task_id","task","status","evidence_or_reason"]; w(status,"CLOSEOUT_STATUS.tsv")

claims=pd.DataFrame([
 ["M04 denominator","GSE250152 contains 58,174 raw features, 35,749 uniquely mapped rows and 27,506 eligible canonical GeneIDs.","01_M04_flow/helper_flow_validation_retry/COUNT_FLOW.tsv","all rows","Supplementary Fig. S11; Table S17","report correction only"],
 ["M04 cache","The mapping cache is a general annotation product from the fixed NCBI GeneInfo snapshot, not a blood-selected expression universe.","01_M04_flow/ANNOTATION_CACHE_VS_FIXED_SNAPSHOT.tsv","all 38,551 pairs","Table S17","no score or DE rerun required"],
 ["M02 DMR provenance","The installed workflow exactly reproduced 167 PROPGER and 660 PROPENG DMR coordinates.","02_M02_audit/DMR_exact_reproduction_v2/DMR_REPRODUCTION_AUDIT.tsv","PROPGER and PROPENG","Supplementary Table S18","computational reproduction, not validation"],
 ["M02 fixed regions","No one of 167 fixed PROPGER regions met BH q<0.10 in PROPENG.","03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv","all 167 regions","Supplementary Fig. S10; Table S18","known-results retrospective average-region test"],
 ["M03 ROI","Seven of nine evaluable programs differed between author Nageotte and neuronal ROIs at q=0.04018.","04_M03_spatial/helper_spatial_oracle/SPATIAL_AUDIT.tsv","P1-P10 including P9 NE","Figure 7; Table S19","diabetes history, not explicit DPN diagnosis"],
 ["M03 arms","The two score arms contributed differently across programs.","04_M03_spatial/M03_TWO_ARM_TESTS.tsv","all 18 eligible tests","Figure 7; Table S20","separate family; not independent validation"],
 ["M03 coordinates","Displayed coordinates are real array coordinates with author ROI labels.","04_M03_spatial/M03_REAL_COORDINATES_WITH_AUTHOR_ROI.tsv.gz","all 16 sections","Figure 7; Supplementary Fig. S13","not a histology image; distance gradient NE"],
 ["Block deletion","632 fixed-module deletions were evaluable; sign preservation varies by endpoint and program.","05_block_summary/helper_block_summary/BLOCK_SUMMARY.tsv","all rows","Supplementary Fig. S12; Table S21","evaluability is not robustness"],
 ["M05 DE","AH had 1,288 q<0.05 genes, MG and diagnosis-by-muscle interaction had none.","06_remaining_status/M05_DE_SUMMARY.tsv","all planned contrasts","Table S22","AH significance does not establish interaction"],
 ["Access boundary","FinnGen and ALC/EDC tasks have no authorized data result.","00_admin/CLOSEOUT_STATUS.tsv","C08","Table S23","blocked access, not biological negative"],
])
claims.columns=["claim_id","claim","result_file","record_scope","display","boundary"]; w(claims,"CLAIM_EVIDENCE_REGISTER.tsv")

for name in ("M04_DE_SUMMARY.tsv","M05_DE_SUMMARY.tsv"):
    d=pd.read_csv(T/"06_whole_gene_DE/results"/name,sep="\t")
    d.to_csv(C/"06_remaining_status"/name,sep="\t",index=False,na_rep="NA")

baseline=("# Reconstructed reported baseline\n\nThis file was reconstructed contemporaneously on 12 September 2026 because the CodeKit supplied to the workspace did not contain its referenced `REPORTED_BASELINE.md`. It is not the missing original and is not backdated.\n\n"
"Authoritative baseline: `analysis_v4/targeted_repair_2026-09-12`, the final v3.1 strict background, `strict_background_addendum/consolidated_results`, and `FINAL_VIEW_INDEX.tsv`. Historical v2, Batch01–04, the first repair implementation and the failed closeout attempts remain read only. P9 remains ineligible for a full RNA score because its repaired source-down arm contains seven genes.\n\n"
"The targeted-repair branch already corrected the M05 rank universe, repaired fixed module mapping, completed M04/M05 whole-gene analyses, recovered GSE286347 and GSE295206 public inputs, and produced the existing ROI and DMR results. The closeout reruns only provenance checks and analyses explicitly required by C01–C08.\n")
(C/"00_admin/REPORTED_BASELINE_RECONSTRUCTED.md").write_text(baseline,encoding="utf-8")

cfg={"project":"phase0_6_human_dpn_stage_projection","closeout_root":str(C),"locked_utc":datetime.now(timezone.utc).isoformat(),"baseline":str(T),"strict_members":str(P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"),"python":sys.version,"platform":platform.platform(),"scientific_boundaries":{"no_new_selection":True,"P9_full_RNA":"NOT_EVALUABLE","fixed_region_role":"RETROSPECTIVE_FIXED_REGION_RETEST","controlled_access":"NOT_AUTHORIZED","public_release":"NOT_AUTHORIZED"}}
(C/"00_admin/LOCKED_CLOSEOUT_CONFIG.json").write_text(json.dumps(cfg,indent=2),encoding="utf-8")

readme=("# DPN v4 targeted closeout\n\nThis directory is a bounded closeout of the existing repaired analysis. It does not restart v4 or overwrite historical results.\n\n"
"## Main outcomes\n\n- C01 corrected the M04 report denominator; the fixed annotation cache was reproduced exactly and no M04 matrix change was indicated.\n- C02 exactly reproduced the 167/660 DMR coordinate lists under the installed implementation.\n- C03 tested all 167 fixed German regions in British participants; none met BH q<0.10.\n- C05 completed donor-level spatial ROI display, arms, contributions, coordinates and leave-one-donor-out summaries; the distance gradient remains NE.\n- C06 summarizes every one of 632 fixed-module deletions.\n- C08 remains blocked by access and linkage.\n\nSee `CLOSEOUT_STATUS.tsv` and `CLAIM_EVIDENCE_REGISTER.tsv` for exact boundaries.\n")
(C/"README.md").write_text(readme,encoding="utf-8")
print("admin files written")
