import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"

# M02-F3 uses comparable beta-scale effects from discovery and held-out cohorts.
disc=pd.read_csv(P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/M02_results/M02_PROPGER_DMR_ALL.tsv",sep="\t")
test=pd.read_csv(RUN/"02_GSE286347/results/M02_FIXED_REGION_ADJUSTED_RETEST.tsv",sep="\t")
disc[["chr","rest"]]=disc.coord.str.split(":",n=1,expand=True);disc[["start","end"]]=disc.rest.str.split("-",n=1,expand=True);disc.start=disc.start.astype(int);disc.end=disc.end.astype(int)
x=test.merge(disc[["chr","start","end","meandiff"]],on=["chr","start","end"],how="left",validate="one_to_one").rename(columns={"meandiff":"PROPGER_discovery_mean_beta_difference","effect_beta_adjusted":"PROPENG_sex_adjusted_beta_difference"})
x.to_csv(RUN/"02_GSE286347/results/M02_FIXED_REGION_EFFECT_COMPARISON.tsv",sep="\t",index=False)
fig,ax=plt.subplots(figsize=(5.5,5.2));ax.scatter(x.PROPGER_discovery_mean_beta_difference,x.PROPENG_sex_adjusted_beta_difference,s=18,alpha=.65,color="#486b8a");ax.axhline(0,color="black",lw=.6);ax.axvline(0,color="black",lw=.6);ax.set_xlabel("PROPGER discovery mean beta difference");ax.set_ylabel("PROPENG sex-adjusted beta difference");fig.tight_layout();fig.savefig(RUN/"02_GSE286347/figures/M02_F3_fixed_regions.png",dpi=300);plt.close(fig)

# M05-F3 displays every locked burden outcome rather than selecting one.
z=pd.read_csv(RUN/"05_GSE295206/results/M05_NAGEOTTE_BURDEN_PROGRAM_ASSOCIATIONS.tsv",sep="\t");hm=z.pivot(index="program",columns="outcome",values="rho")
fig,ax=plt.subplots(figsize=(7,6));sns.heatmap(hm,cmap="vlag",center=0,vmin=-1,vmax=1,annot=True,fmt=".2f",linewidths=.3,linecolor="white",ax=ax,cbar_kws={"label":"Spearman rho"});ax.set_xlabel("");ax.set_ylabel("");fig.tight_layout();fig.savefig(RUN/"05_GSE295206/figures/M05_F3_program_burden.png",dpi=300);plt.close(fig)
print(f"M02 regions={len(x)} missing discovery effects={x.PROPGER_discovery_mean_beta_difference.isna().sum()}; M05 rows={len(z)}")
