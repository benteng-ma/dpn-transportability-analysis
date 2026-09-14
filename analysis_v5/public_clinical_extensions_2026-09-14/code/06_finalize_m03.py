import os
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"
OUT=RUN/"03_GSE14806x/results"; LOG=RUN/"03_GSE14806x/logs"
DE_RNA=P/"analysis_v3/batch02/results/de/DE_ALL_GENES_C02.tsv.gz"

def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);finite=np.isfinite(p)
    pf=p[finite];o=np.argsort(pf);z=np.minimum.accumulate((pf[o]*len(pf)/np.arange(1,len(pf)+1))[::-1])[::-1]
    qf=np.empty(len(pf));qf[o]=np.minimum(z,1);q[finite]=qf;return q

# Correct only the BH layer after detecting non-finite Welch tests; raw effects
# and P values are retained exactly. The family is every finite eligible CpG.
de=pd.read_csv(OUT/"M03_RRBS_DIFFERENTIAL_RESULTS.tsv.gz",sep="\t")
de["BH_q"]=bh(de.P)
de.to_csv(OUT/"M03_RRBS_DIFFERENTIAL_RESULTS.tsv.gz",sep="\t",index=False,compression="gzip")

# Join the already completed Batch02 C02 extreme-group RNA model at gene level.
# This is group-level cross-assay comparison only; no patient pairing is inferred.
g=pd.read_csv(OUT/"M03_CROSSOMIC_GENE_EFFECTS.tsv.gz",sep="\t")
r=pd.read_csv(DE_RNA,sep="\t")[["gene","effect","moderated_t","P","q"]].rename(columns={"gene":"gene_symbol","effect":"rna_effect","moderated_t":"rna_stat","P":"rna_P","q":"rna_q"})
g=g.drop(columns=[c for c in ["rna_effect","rna_stat"] if c in g]).merge(r,on="gene_symbol",how="left",validate="many_to_one")
g["same_direction"]=np.where(g.rna_effect.notna(),np.sign(g.rrbs_effect)==np.sign(g.rna_effect),np.nan)
g["pairing"]="GROUP_LEVEL_ONLY_NO_PATIENT_KEY"
g.to_csv(OUT/"M03_CROSSOMIC_GENE_EFFECTS.tsv.gz",sep="\t",index=False,compression="gzip")

rows=[]
for annotation in ["promoter","gene_body"]:
    x=g[(g.annotation==annotation)&g.rna_effect.notna()&g.rrbs_effect.notna()].copy()
    rho=stats.spearmanr(x.rrbs_effect,x.rna_effect,nan_policy="omit").statistic
    rows.append({
      "annotation":annotation,"genes":len(x),"spearman_rho":rho,
      "correlation_P":"NOT_COMPUTED_GENES_NOT_INDEPENDENT",
      "same_direction_n":int(x.same_direction.sum()),"same_direction_fraction":float(x.same_direction.mean()),
      "inverse_direction_n":int((~x.same_direction.astype(bool)).sum()),"inverse_direction_fraction":float((~x.same_direction.astype(bool)).mean()),
      "directional_heuristic":"inverse expected only for promoter; descriptive, non-causal",
      "patient_pairing":"NOT_USED_NO_VERIFIED_KEY","P_pooling":"NONE","status":"DESCRIPTIVE_GROUP_LEVEL"})
concord=pd.DataFrame(rows)
concord.to_csv(OUT/"M03_CROSSOMIC_GENE_CONCORDANCE.tsv",sep="\t",index=False)

# Make the RNA program endpoint explicit instead of relying on row order.
ps=pd.read_csv(OUT/"M03_PROGRAM_CROSSOMIC_SUMMARY.tsv",sep="\t")
rr=pd.read_csv(OUT/"M03_REUSED_CORRECTED_RNA_RESULTS.tsv",sep="\t")
ordr=rr[rr.test_family.eq("GSE148059_dMFD_ordinal")][["module_id","effect","permutation_two_sided_bh_q","n_total","status"]].rename(columns={"effect":"RNA_ordinal_rho","permutation_two_sided_bh_q":"RNA_ordinal_q","n_total":"RNA_n","status":"RNA_ordinal_status"})
ext=rr[rr.test_family.eq("GSE148059_Degenerator_vs_Regenerator")][["module_id","effect","permutation_two_sided_bh_q","n_total","status"]].rename(columns={"effect":"RNA_extreme_Hedges_g","permutation_two_sided_bh_q":"RNA_extreme_q","n_total":"RNA_extreme_n","status":"RNA_extreme_status"})
endpoint_cols=["RNA_group_effect","RNA_q","RNA_ordinal_rho","RNA_ordinal_q","RNA_n","RNA_ordinal_status","RNA_extreme_Hedges_g","RNA_extreme_q","RNA_extreme_n","RNA_extreme_status"]
endpoint_cols += [c for c in ps.columns if c.endswith("_x") or c.endswith("_y")]
ps=ps.drop(columns=[c for c in endpoint_cols if c in ps]).merge(ordr,on="module_id",how="left").merge(ext,on="module_id",how="left")
ps["directions_compatible"]=np.where(ps.RNA_extreme_Hedges_g.notna()&ps.RRBS_promoter_median_effect.notna(),np.where(np.sign(ps.RNA_extreme_Hedges_g)==-np.sign(ps.RRBS_promoter_median_effect),"INVERSE_HEURISTIC_COMPATIBLE","NOT_INVERSE_HEURISTIC"),"NOT_EVALUABLE")
ps.to_csv(OUT/"M03_PROGRAM_CROSSOMIC_SUMMARY.tsv",sep="\t",index=False)

# Explicitly classify every tested locus; overlapping transcript annotations are
# retained in the detailed mapping while this table supplies one hierarchy label.
ann=pd.read_csv(OUT/"M03_RRBS_GENE_ANNOTATION.tsv.gz",sep="\t")
cls=ann.groupby("key").annotation.agg(lambda s:"promoter" if (s=="promoter").any() else "gene_body").rename("locus_class")
locus=de[["key","chr","position"]].drop_duplicates().merge(cls,left_on="key",right_index=True,how="left")
locus["locus_class"]=locus.locus_class.fillna("other_or_ambiguous")
locus.to_csv(OUT/"M03_RRBS_LOCUS_CLASSIFICATION.tsv.gz",sep="\t",index=False,compression="gzip")

# Required cross-omic and program evidence figures use the fixed extreme-group
# RNA contrast and the already calculated RRBS promoter summaries.
FIG=RUN/"03_GSE14806x/figures"
prom=g[(g.annotation=="promoter")&g.rna_effect.notna()&g.rrbs_effect.notna()].copy()
fig,ax=plt.subplots(figsize=(7,5));ax.hexbin(prom.rrbs_effect,prom.rna_effect,gridsize=55,mincnt=1,cmap="Blues");ax.axhline(0,color="black",lw=.6);ax.axvline(0,color="black",lw=.6);ax.set_xlabel("Promoter methylation M-value difference");ax.set_ylabel("RNA rlog-scale difference");fig.tight_layout();fig.savefig(FIG/"M03_F2_RNA_PROMOTER_EFFECTS.png",dpi=300);plt.close(fig)
hm=ps.set_index("program")[["RNA_extreme_Hedges_g","RRBS_promoter_median_effect"]].rename(columns={"RNA_extreme_Hedges_g":"RNA Hedges g","RRBS_promoter_median_effect":"Promoter methylation median M difference"})
fig,ax=plt.subplots(figsize=(6,6));sns.heatmap(hm,cmap="vlag",center=0,linewidths=.35,linecolor="white",ax=ax,cbar_kws={"label":"Effect (assay-specific scale)"});ax.set_xlabel("");ax.set_ylabel("");fig.tight_layout();fig.savefig(FIG/"M03_F3_PROGRAM_EVIDENCE_MATRIX.png",dpi=300);plt.close(fig)

tr=pd.read_csv(OUT/"M03_RRBS_ORDERED_TREND.tsv.gz",sep="\t")
summary={
 "BH_repair_scope":"finite eligible CpGs only; raw effects and P unchanged",
 "DE_total_rows":len(de),"DE_finite_P":int(de.P.notna().sum()),"DE_q_lt_0_05":int((de.BH_q<.05).sum()),
 "DE_min_q":float(de.BH_q.min()),"trend_total_rows":len(tr),"trend_q_lt_0_05":int((tr.BH_q<.05).sum()),
 "trend_positive_q_lt_0_05":int(((tr.BH_q<.05)&(tr.M_slope_per_pathology_step>0)).sum()),
 "trend_negative_q_lt_0_05":int(((tr.BH_q<.05)&(tr.M_slope_per_pathology_step<0)).sum()),
 "crossomic_promoter_genes":len(prom),"crossomic_promoter_spearman_rho":float(concord.loc[concord.annotation.eq('promoter'),'spearman_rho'].iloc[0]),"crossomic_promoter_inverse_direction_fraction":float(concord.loc[concord.annotation.eq('promoter'),'inverse_direction_fraction'].iloc[0]),
 "locus_classes":locus.locus_class.value_counts().to_dict()
}
(LOG/"M03_POSTPROCESS_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
