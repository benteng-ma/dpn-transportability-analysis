import os
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; OUT=RUN/"03_GSE14806x/results"; FIG=RUN/"03_GSE14806x/figures"; LOG=RUN/"03_GSE14806x/logs"
RAW=P/"analysis_v3/data/raw/GSE148060/extracted"
meta=pd.read_csv(OUT/"M03_RRBS_SAMPLE_AUDIT.tsv",sep="\t")
order={"Regenerator":0,"Intermediate":1,"Degenerator":2};rows=[]
for x in meta.itertuples():
    d=pd.read_csv(RAW/x.source_file,sep="\t",header=None,usecols=[3],names=["percent"],dtype={"percent":"float32"})
    v=d.percent.to_numpy(float);v=v[np.isfinite(v)]
    rows.append({"sample_id":x.sample_id,"group":x.group,"ordinal":order[x.group],"reported_loci":len(v),"median_percent_methylation":float(np.median(v)),"mean_percent_methylation":float(np.mean(v)),"fraction_zero":float(np.mean(v==0)),"fraction_hundred":float(np.mean(v==100))})
q=pd.DataFrame(rows);q.to_csv(OUT/"M03_RRBS_GLOBAL_SAMPLE_QC.tsv",sep="\t",index=False)
rho,p=stats.spearmanr(q.ordinal,q.median_percent_methylation);kr=stats.kruskal(*[q.loc[q.group.eq(g),"median_percent_methylation"] for g in order])
summary={"n":len(q),"median_methylation_ordinal_rho":float(rho),"median_methylation_spearman_P_descriptive":float(p),"median_methylation_kruskal_P_descriptive":float(kr.pvalue),"loci_min":int(q.reported_loci.min()),"loci_max":int(q.reported_loci.max()),"interpretation":"sample-level QC; not a program test and not covariate-adjusted"}
(LOG/"M03_GLOBAL_QC_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
fig,axs=plt.subplots(1,2,figsize=(9,4));sns.boxplot(data=q,x="group",y="median_percent_methylation",order=list(order),color="#9eb6c8",ax=axs[0]);sns.stripplot(data=q,x="group",y="median_percent_methylation",order=list(order),color="black",size=3,ax=axs[0]);axs[0].tick_params(axis="x",rotation=25);axs[0].set_xlabel("");axs[0].set_ylabel("Sample median methylation (%)");sns.boxplot(data=q,x="group",y="reported_loci",order=list(order),color="#c9b3a6",ax=axs[1]);sns.stripplot(data=q,x="group",y="reported_loci",order=list(order),color="black",size=3,ax=axs[1]);axs[1].tick_params(axis="x",rotation=25);axs[1].set_xlabel("");axs[1].set_ylabel("Reported CpG loci");fig.tight_layout();fig.savefig(FIG/"M03_F4_GLOBAL_SAMPLE_QC.png",dpi=300);plt.close(fig)
print(json.dumps(summary,indent=2))
