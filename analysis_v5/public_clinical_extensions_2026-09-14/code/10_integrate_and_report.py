import os
from pathlib import Path
import hashlib, json, shutil
import numpy as np
import pandas as pd

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; OUT=RUN/"06_integration"; OUT.mkdir(exist_ok=True)
PROGRAMS=[f"P{i}" for i in range(1,11)]
LABELS={"original_early_allcell":"P1","original_late_allcell":"P2","original_late_neuron":"P3","original_severity":"P4","original_xenium":"P5","late_shared_concordant_neuronal_core":"P6","late_neuron_residual":"P7","late_allcell_residual":"P8","severity_neuron_shared_concordant_core":"P9","severity_neuron_residual":"P10"}
def direction(x): return "positive" if pd.notna(x) and x>0 else ("negative" if pd.notna(x) and x<0 else "NA")
def cell(status,effect=None,q=None,n=None,measurement="",independence="",boundary=""):
    return json.dumps({"status":status,"effect_direction":direction(effect),"effect_size":None if pd.isna(effect) else float(effect),"q":None if pd.isna(q) else float(q),"n_people":n,"measurement":measurement,"independence_note":independence,"claim_boundary":boundary},separators=(",",":"))

hist=pd.read_csv(P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/results/ASSOCIATION_IMPACT_ALL.tsv",sep="\t");hist["program"]=hist.module_id.map(LABELS)
roi=pd.read_csv(P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/M03_results/M03_PROGRAM_TESTS.tsv",sep="\t")
m1=RUN/"01_GSE302658/results"; base=pd.read_csv(m1/"BASELINE_SYMPTOM_ASSOCIATIONS.tsv",sep="\t");long=pd.read_csv(m1/"LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv",sep="\t");inter=pd.read_csv(m1/"TREATMENT_INTERACTION_OMNIBUS.tsv",sep="\t")
m2=pd.read_csv(RUN/"02_GSE286347/results/M02_PROPENG_BURDEN_ASSOCIATIONS.tsv",sep="\t")
m3=pd.read_csv(RUN/"03_GSE14806x/results/M03_PROGRAM_CROSSOMIC_SUMMARY.tsv",sep="\t")
m5=pd.read_csv(RUN/"05_GSE295206/results/M05_NAGEOTTE_BURDEN_PROGRAM_ASSOCIATIONS.tsv",sep="\t")

wide=[]; longrows=[]
for p in PROGRAMS:
    row={"program":p}
    z=hist[(hist.endpoint=="original_JCI_sural_DPN_control")&(hist.program==p)]
    if len(z) and z.status.iloc[0]=="COMPLETED":
        x=z.iloc[0];st="SUPPORTED" if x.q_BH_common_method<.05 else "EXPLORATORY" if x.q_BH_common_method<.10 else "UNSUPPORTED";row["sural_DPN_case_association"]=cell(st,x.mean_difference,x.q_BH_common_method,int(x.n_total),"fixed RNA percentile-rank score","same 12 JCI participants as historical corrected analysis","association, not diagnosis or causality")
    else: row["sural_DPN_case_association"]=cell("NE",measurement="fixed RNA score",boundary="fewer than 10 measured genes per arm")
    z=roi[roi.program==p]
    if len(z):
        x=z.iloc[0];st="EXPLORATORY" if x.q<.05 else "UNSUPPORTED";row["DRG_Nageotte_ROI_association"]=cell(st,x.effect,x.q,int(x.n_donors),"Nageotte minus adjacent neuronal ROI score","six donors; sections not independent people","diabetes history; explicit DPN diagnosis not established")
    else: row["DRG_Nageotte_ROI_association"]=cell("NE",measurement="spatial ROI score",boundary="P9 coverage failure")
    z=base[base.program==p]
    if len(z):
        x=z.sort_values("BH_q").iloc[0];row["GSE302658_baseline_symptom"]=cell("UNSUPPORTED",x.beta_program,x.BH_q,int(x.n),f"best-q among 8 locked domains: {x.symptom}","trial participants","no baseline domain survived BH")
    else: row["GSE302658_baseline_symptom"]=cell("NE",boundary="program coverage below locked two-arm threshold")
    z=long[long.program==p]
    if len(z):
        x=z.sort_values("BH_q").iloc[0];n_sig=int((z.BH_q<.05).sum());st="SUPPORTED" if n_sig else "UNSUPPORTED";row["GSE302658_within_patient_symptom_change"]=cell(st,x.beta_delta_program,x.BH_q,int(x.n_paired),f"change association; {n_sig}/8 domains q<0.05; best {x.symptom}","101 paired trial participants","contemporaneous covariation, not prediction or treatment mediation")
    else: row["GSE302658_within_patient_symptom_change"]=cell("NE",boundary="program coverage below locked threshold")
    z=inter[inter.program==p]
    if len(z):
        x=z.sort_values("BH_q_omnibus").iloc[0];row["GSE302658_randomized_treatment_interaction"]=cell("UNSUPPORTED",None,x.BH_q_omnibus,int(x.n),f"2-df treatment interaction; best {x.symptom}","randomized trial participants","specific AZD2423 trial; no general treatment conclusion")
    else: row["GSE302658_randomized_treatment_interaction"]=cell("NE",boundary="program coverage below locked threshold")
    z=m2[m2.score==p]
    if len(z) and pd.notna(z.P.iloc[0]):
        x=z.iloc[0];row["GSE286347_painful_painless_methylation_burden"]=cell("UNSUPPORTED",x.effect_painful_minus_painless,x.BH_q,int(x.n),f"held-out signed burden from {int(x.n_regions)} fixed PROPGER regions","PROPENG participants distinct from discovery cohort","methylation burden, not RNA-score validation")
    else: row["GSE286347_painful_painless_methylation_burden"]=cell("NE",measurement="program-linked fixed-region burden",boundary="fewer than 5 linked fixed regions")
    z=m3[m3.program==p]
    if len(z) and pd.notna(z.RNA_extreme_Hedges_g.iloc[0]):
        x=z.iloc[0];row["GSE14806x_RNA_pathology_grouping"]=cell("UNSUPPORTED",x.RNA_extreme_Hedges_g,x.RNA_extreme_q,int(x.RNA_extreme_n),"GSE148059 Degenerator versus Regenerator fixed RNA score","RNA and RRBS participant overlap unknown","processed-rlog exploratory association")
    else: row["GSE14806x_RNA_pathology_grouping"]=cell("NE",measurement="fixed RNA score",boundary="P9 two-arm coverage failure")
    row["GSE14806x_RRBS_pathology_grouping"]=cell("NE",measurement="program promoter-methylation enrichment",independence="RNA-RRBS pairing not used",boundary="bias-corrected gene-set method unavailable; assay-wide CpG trend is not program support")
    row["JCI_continuous_pathology"]=cell("NE",measurement="donor-linked continuous pathology",independence="JCI supplement audit",boundary="only categorical axonal-loss grades publicly linked; already analyzed")
    z=m5[m5.program==p]
    if len(z):
        x=z.sort_values("BH_q").iloc[0];row["GSE295206_Nageotte_burden"]=cell("UNSUPPORTED",x.rho,x.BH_q,int(x.n_donors),f"Spearman with annotated-spot burden; best {x.outcome}","six donors; equal donor weight","spot fraction is not histologic area; no prediction claim")
    else: row["GSE295206_Nageotte_burden"]=cell("NE",measurement="donor burden association",boundary="P9 coverage failure")
    wide.append(row)
    for question,value in row.items():
        if question!="program": longrows.append({"program":p,"evidence_question":question,**json.loads(value)})
pd.DataFrame(wide).to_csv(OUT/"CLINICAL_EVIDENCE_MATRIX.tsv",sep="\t",index=False)
pd.DataFrame(longrows).to_csv(OUT/"CLINICAL_EVIDENCE_MATRIX_LONG.tsv",sep="\t",index=False)

# Multiplicity and status registers.
families=[
 ["M01_baseline",64,"BH","8 programs x 8 baseline domains"],["M01_longitudinal",64,"BH","8 programs x 8 change domains"],["M01_interaction",64,"BH","8 programs x 8 omnibus interactions"],
 ["M02_PROPGER_sex_batch",840040,"BH","all eligible CpGs"],["M02_PROPENG_sex_batch",839953,"BH","all eligible CpGs"],["M02_PROPENG_fixed_regions_batch",167,"BH","all frozen regions"],["M02_program_burden",2,"BH","only P2/P8 meeting >=5 regions; global separate"],
 ["M03_RRBS_extreme",1261118,"BH","all finite eligible CpGs"],["M03_RRBS_ordered",1446560,"BH","all finite eligible CpGs"],["M03_program_methylation",0,"NE","bias-corrected method unavailable"],
 ["M05_burden_each_outcome",9,"BH","separate family for whole/adjacent/contrast"],["M05_arms",18,"BH","all eligible fixed arms"]]
pd.DataFrame(families,columns=["family","tests","adjustment","scope"]).to_csv(OUT/"MULTIPLICITY_FAMILIES.tsv",sep="\t",index=False)
status=pd.DataFrame([
 ["M01","COMPLETED","104 participants; baseline, longitudinal and randomized interaction families completed","P5/P9 NE from fixed coverage"],
 ["M02","COMPLETED_WITH_NE","315 public samples audited; 231 DPN analyzed; batch sensitivity and fixed burden completed","cell-composition sensitivity NE"],
 ["M03","COMPLETED_WITH_NE","53 RRBS samples and corrected RNA results compared at group level","patient pairing blocked; methylation program enrichment NE"],
 ["M04","COMPLETED_AUDIT_ONLY","28 public supplements audited","continuous/composition/paired pathology NE"],
 ["M05","COMPLETED","16 sections, 7 DRGs, 6 donors; burden, arm, marker proxy and topological summaries completed","P9 NE; no corrected burden association"],
 ["M06","COMPLETED","program-by-question evidence matrix and bounded draft completed","no P-value pooling or evidence score"]],columns=["module","status","completed","boundary"])
status.to_csv(OUT/"MODULE_STATUS.tsv",sep="\t",index=False)

# Independence and sample-key register (only documented keys; no inferred crosswalk).
ind=pd.DataFrame([
 ["GSE302658",104,"participant_id","within-study baseline/follow-up direct key","independent of other sources unknown"],
 ["GSE286347_PROPGER",139,"GEO sample/public cohort key","distinct discovery cohort","not merged with PROPENG"],["GSE286347_PROPENG",92,"GEO sample/public cohort key","distinct held-out cohort","not merged with PROPGER"],
 ["GSE148059_RNA",77,"GEO sample key","group-level","overlap with RRBS unknown"],["GSE148060_RRBS",53,"GEO sample key","group-level","no RNA pairing"],
 ["JCI184075_sural",12,"documented sample label","historical corrected cohort","same participants for case and some morphology analyses"],
 ["GSE295206",6,"documented donor hierarchy","16 sections/7 DRGs collapsed to donors","DPN diagnosis not asserted"]],columns=["source","n_people_or_donors","key_basis","within_source_relation","cross_source_boundary"])
ind.to_csv(OUT/"INDEPENDENCE_REGISTER.tsv",sep="\t",index=False)

# Trend magnitude audit preserves the locked inferential result while exposing its
# two-scale effect distribution; no threshold is used to redefine significance.
tr=pd.read_csv(RUN/"03_GSE14806x/results/M03_RRBS_ORDERED_TREND.tsv.gz",sep="\t");sig=tr[tr.BH_q<.05]
pd.DataFrame([{"absolute_M_slope_cut":c,"q_lt_0_05_total":len(sig),"n_at_or_above_cut":int((sig.M_slope_per_pathology_step.abs()>=c).sum()),"use":"descriptive_magnitude_audit_not_filter"} for c in [1e-6,1e-5,1e-4,1e-3,.01,.1,.5,1]]).to_csv(OUT/"M03_TREND_EFFECT_MAGNITUDE_AUDIT.tsv",sep="\t",index=False)

claims=pd.DataFrame([
 ["C01","Baseline fixed RNA programs did not associate with any of eight symptom domains after BH correction.","01_GSE302658/results/BASELINE_SYMPTOM_ASSOCIATIONS.tsv","UNSUPPORTED","No equivalence claim."],
 ["C02","Within-person changes in P1, P4 and P10 covaried with selected symptom changes after BH correction.","01_GSE302658/results/LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv","SUPPORTED","Contemporaneous association; not prediction, mediation or causality."],
 ["C03","No fixed program modified AZD2423 treatment effects after BH correction.","01_GSE302658/results/TREATMENT_INTERACTION_OMNIBUS.tsv","UNSUPPORTED","Specific trial and analyzable programs only."],
 ["C04","One PROPGER CpG emerged only after Sentrix adjustment and did not replicate in PROPENG.","02_GSE286347/results/M02_PROPGER_EWAS_SEX_BATCH_ALL.tsv.gz","EXPLORATORY","Batch-sensitive, single-source signal; no biomarker claim."],
 ["C05","No frozen PROPGER region or held-out program burden was supported in PROPENG.","02_GSE286347/results/M02_FIXED_REGION_BATCH_ADJUSTED_RETEST.tsv","UNSUPPORTED","P1/P3-P7/P9/P10 burdens mostly NE by fixed >=5 rule."],
 ["C06","RRBS showed a broad ordered group trend but no extreme-group CpG after BH; fixed-program methylation enrichment was NE.","03_GSE14806x/results/M03_RRBS_ORDERED_TREND.tsv.gz","EXPLORATORY","No patient pairing; strong directional/numerical structure; not program validation."],
 ["C07","Public JCI supplements contain categorical morphology but no new linked continuous pathology or composition endpoint.","04_JCI184075/results/M04_DIRECT_KEY_AUDIT.tsv","NE","Historical categorical analysis retained."],
 ["C08","Donor-level Nageotte annotated-spot burden did not track any fixed program after BH.","05_GSE295206/results/M05_NAGEOTTE_BURDEN_PROGRAM_ASSOCIATIONS.tsv","UNSUPPORTED","n=6; spot fraction not area."],
 ["C09","Some marker-expression proxy correlations were nominally large but remain descriptive and unadjusted as a family.","05_GSE295206/results/M05_SPATIAL_MARKER_PROXY_SUMMARY.tsv","EXPLORATORY","Expression proxies are not cell fractions."]],columns=["claim_id","claim","result_file","status","boundary"])
claims.to_csv(OUT/"CLINICAL_EXTENSION_CLAIM_REGISTER.tsv",sep="\t",index=False)

panel=[]
for m in ["01_GSE302658","02_GSE286347","03_GSE14806x","04_JCI184075","05_GSE295206"]:
    for f in sorted((RUN/m/"figures").glob("*.png")): panel.append({"figure":f.name,"module":m.split('_')[0],"source_directory":str((RUN/m/"results").relative_to(RUN)),"status":"actual output","boundary":"See module results and claim register"})
pd.DataFrame(panel).to_csv(OUT/"PANEL_SOURCE_REGISTER.tsv",sep="\t",index=False)

report="""# DPN public clinical extensions review report

## Scope and execution

This is a post-correction, public-data-only exploratory extension of the fixed v3.1 programs. It does not alter program membership, direction, the locked coverage rule, or historical results. No author was contacted, no controlled-access data were requested, no patient crosswalk was inferred, and no P values were pooled across sources.

## Direct answers to the eight clinical questions

1. **Baseline symptoms:** No. Across 64 program-domain tests in 102 baseline participants, no association survived BH correction (minimum q = 0.9987).
2. **Within-patient symptom change:** Yes, within the prespecified exploratory family. Twelve of 64 program-domain associations had q < 0.10 and eight had q < 0.05. These involved P1, P4 and P10 and the NPSI total, superficial burning, deep pressing and paresthesia/dysesthesia domains. The strongest family q was 0.0231. Leave-one-participant-out effects retained direction in every q < 0.10 result. These are contemporaneous molecular-score/symptom-change associations, not baseline prediction, treatment mediation, or causality.
3. **Treatment-effect modification:** No. None of the 64 two-degree-of-freedom AZD2423 interaction tests survived correction (minimum q = 0.9981). This does not generalize beyond this trial or the eight coverage-eligible programs.
4. **Adjusted painful-versus-painless methylation:** The expanded adjustment did not establish a reproducible clinical methylation result. Sentrix adjustment yielded one PROPGER CpG, cg00563832 (M-value difference -0.1772; P = 4.94e-8; q = 0.0415), whereas the same CpG was unsupported in PROPENG (effect -0.0444; P = 0.249). No PROPENG CpG, none of the 167 fixed regions, and neither evaluable program-linked burden (P2/P8) survived correction. Cell-composition sensitivity was NE because a validated EPIC reference implementation was unavailable locally.
5. **RRBS pathology support:** Only at a broad assay-level group-trend tier, not as fixed-program support. Among 1,261,118 finite extreme-group tests, no CpG survived BH (minimum q = 0.8058). The ordered three-class model produced 101,009 q < 0.05 loci, mostly positive, but 98,194 of those had an absolute M-value slope below 1e-5 and 2,815 were at least 0.01. Sample-level median methylation did not show an ordinal shift (rho = 0.113; descriptive P = 0.419), so the extensive one-direction pattern requires technical and locus-level caution. RNA-versus-promoter methylation effect magnitudes were essentially uncorrelated (rho = -0.0032); gene-wise correlation P was deliberately not assigned. Bias-corrected methylation program enrichment was not available and remains NE.
6. **JCI continuous pathology:** No. All 28 public supplements (13 XLSX) were audited. Direct sample labels, age/sex and categorical axonal-loss grades are public, but no new donor-linked continuous pathology, matched composition measure, or quantitative pathology in both paired nerves was found. The corrected categorical moderate-versus-severe analysis was retained without rerunning it for a different P value.
7. **Nageotte burden:** No after correction. In six donors, none of nine programs related to author-labelled Nageotte spot burden across whole-section, adjacent or ROI-contrast scores after BH. The smallest exact P was 0.0583 for P3 adjacent score (rho = 0.829; q = 0.525). The burden is an annotated-spot fraction, not histologic area. Marker-expression proxy correlations are contextual only and are not cell fractions.
8. **Clinical relevance versus context:** The within-person symptom-change result is the clearest new clinical-relevance increment, but it is exploratory and nonpredictive. The painful-versus-painless methylation analysis directly addresses a clinical phenotype yet remains single-source/batch-sensitive and unreplicated. RRBS, spatial ROI, marker proxies and topology mainly add pathology or tissue context. They do not create a diagnostic test, treatment-selection biomarker, prognostic model, or causal mechanism.

## Overall interpretation

The extension improves the manuscript's clinical anchoring by showing that selected fixed neural programs covary with symptom change within a randomized-trial cohort while failing to predict baseline symptom burden or modify the tested treatment response. The wider public methylation and spatial analyses narrow the claims: blood methylation signals were not reproducible across cohorts, fixed methylation burdens were unsupported, JCI lacks a new continuous public endpoint, and Nageotte burden associations were not corrected-significant. The evidence therefore remains a bounded translational association study rather than a clinically deployable model.
"""
(OUT/"CLINICAL_EXTENSION_REVIEW_REPORT.md").write_text(report,encoding="utf-8")

draft="""# Manuscript clinical integration draft

## Proposed Results insertion

We next evaluated whether fixed sensory-ganglion programs related to patient-reported neuropathic pain in the public GSE302658 trial cohort. None of 64 baseline program-by-domain associations was supported after false-discovery-rate correction (minimum q=0.999; 102 participants). In contrast, within-participant changes in three programs (P1, P4, and P10) covaried with changes in NPSI total, superficial burning, deep pressing, and paresthesia/dysesthesia scores. Eight of 64 associations met q<0.05 and 12 met q<0.10 (minimum q=0.023; 101 paired participants), with positive coefficients indicating that increasing program scores accompanied increasing symptom scores. All q<0.10 effects retained their direction in leave-one-participant-out analyses. These associations were contemporaneous and do not establish baseline prediction or causal mediation. No program-by-treatment interaction was supported for AZD2423 across the same eight symptom domains (minimum q=0.998).

Public blood methylation provided limited and nonreproducible painful-DPN evidence. Adjustment for sex and directly derived Sentrix identifiers identified one CpG in PROPGER (cg00563832, q=0.041), but the locus was not supported in PROPENG. None of 167 discovery-fixed regions survived correction in PROPENG, and the global fixed-region burden and the two coverage-eligible program-linked burdens were also unsupported. Cell-composition sensitivity could not be evaluated with a validated EPIC reference implementation in the locked local environment.

Group-level comparison of GSE148059 RNA and GSE148060 RRBS products did not assume patient pairing. No extreme Degenerator-versus-Regenerator CpG survived correction. An extensive, predominantly positive ordered methylation trend was observed across the three deposited pathology classes; however, most corrected-significant slopes were numerically minute, sample-level median methylation did not show a corresponding ordinal shift, and promoter-methylation and RNA effect magnitudes were essentially uncorrelated. Because a locus-representation-aware program test could not be implemented, RRBS did not provide fixed-program validation.

Finally, all public JCI184075 supplements were examined for linked continuous pathology and cell-composition variables. Only categorical axonal-loss grades suitable for the previously reported analysis were available. In GSE295206, author-labelled Nageotte spot burden was quantified across 16 sections, seven DRGs, and six donors, but no whole-section, adjacent-region, or ROI-contrast program association survived correction. Fixed marker-expression proxies and graph-step gradients are retained as descriptive tissue context rather than estimates of cell fractions or physical distance.

## Proposed Discussion insertion

The within-person symptom-change associations improve clinical anchoring but should be interpreted conservatively. They were identified after the fixed program definitions were known, they quantify concurrent change rather than prediction, and they were not accompanied by baseline symptom or randomized-treatment-interaction support. The methylation and spatial extensions further argue against a simple transferable biomarker interpretation: a batch-sensitive CpG was not reproduced across painful-DPN cohorts, fixed methylation burdens were negative or not evaluable, and Nageotte burden did not track the programs after correction. Together, the results support selective clinical covariation while defining substantial assay, tissue, and endpoint boundaries.

## Text that should not be used

- The programs predict neuropathic pain progression.
- The programs identify responders to AZD2423 or other treatments.
- RRBS independently validates the RNA programs.
- Nageotte burden proves a DPN mechanism.
- Marker-expression proxies are cell fractions.
"""
(OUT/"MANUSCRIPT_CLINICAL_INTEGRATION_DRAFT.md").write_text(draft,encoding="utf-8")
print(json.dumps({"programs":len(wide),"cells":len(longrows),"claims":len(claims),"panels":len(panel)},indent=2))
