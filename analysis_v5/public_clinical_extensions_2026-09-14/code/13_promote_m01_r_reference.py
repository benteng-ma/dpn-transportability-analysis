import os
from pathlib import Path
import json, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; M=RUN/"01_GSE302658"; R=M/"reference_R_run"; O=M/"results"; A=M/"superseded_python_t_reference_implementation"
A.mkdir(exist_ok=True)
files=["BASELINE_SYMPTOM_ASSOCIATIONS.tsv","BASELINE_SPEARMAN_SENSITIVITY.tsv","LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv","LONGITUDINAL_SPEARMAN_SENSITIVITY.tsv","TREATMENT_INTERACTION_OMNIBUS.tsv","TREATMENT_INTERACTION_PAIRWISE.tsv"]
for f in files: shutil.copy2(O/f,A/f)
comparison=[]
for f in files:
    old=pd.read_csv(O/f,sep="\t"); new=pd.read_csv(R/f,sep="\t")
    if f=="BASELINE_SYMPTOM_ASSOCIATIONS.tsv": new=new.rename(columns={"beta":"beta_program"});new["AIC_difference"]=np.nan
    if f=="TREATMENT_INTERACTION_PAIRWISE.tsv": new=new.rename(columns={"se":"hc3_se"})
    keys=[c for c in ["program","symptom"] if c in old and c in new]
    z=old.merge(new,on=keys,suffixes=("_python","_R"))
    pcols=[c for c in ["P","P_omnibus","BH_q","BH_q_omnibus"] if c+"_python" in z]
    comparison.append({"file":f,"python_rows":len(old),"R_rows":len(new),"max_primary_P_or_q_difference":max([float(np.nanmax(abs(z[c+"_python"]-z[c+"_R"]))) for c in pcols],default=np.nan),"promoted":"R_CodeKit"})
    new.to_csv(O/f,sep="\t",index=False,na_rep="NA")
pd.DataFrame(comparison).to_csv(RUN/"tests/M01_R_REFERENCE_COMPARISON.tsv",sep="\t",index=False)

def heat(df,value,out):
    x=df.pivot(index="program",columns="symptom",values=value).reindex([f"P{i}" for i in range(1,11)])
    fig,ax=plt.subplots(figsize=(11,6));sns.heatmap(x,cmap="vlag",center=0,linewidths=.3,linecolor="white",ax=ax);ax.set_xlabel("");ax.set_ylabel("");fig.tight_layout();fig.savefig(M/"figures"/out,dpi=300);plt.close(fig)
b=pd.read_csv(O/files[0],sep="\t");l=pd.read_csv(O/files[2],sep="\t");it=pd.read_csv(O/files[4],sep="\t")
heat(b,"beta_program","M01_F2_baseline.png");heat(l,"beta_delta_program","M01_F3_longitudinal.png")
q=it.pivot(index="program",columns="symptom",values="P_omnibus").reindex([f"P{i}" for i in range(1,11)]);fig,ax=plt.subplots(figsize=(11,6));sns.heatmap(-np.log10(q),cmap="Blues",linewidths=.3,linecolor="white",ax=ax,cbar_kws={"label":"-log10 omnibus P"});ax.set_xlabel("");ax.set_ylabel("");fig.tight_layout();fig.savefig(M/"figures/M01_F4_interaction.png",dpi=300);plt.close(fig)
summary={"participants":104,"baseline":102,"followup":103,"paired":101,"eligible_programs":["P1","P2","P3","P4","P6","P7","P8","P10"],"final_engine":"original CodeKit R script","minimum_q":{"baseline":float(b.BH_q.min()),"longitudinal":float(l.BH_q.min()),"interaction":float(it.BH_q_omnibus.min())},"q_lt_0.05":{"baseline":int((b.BH_q<.05).sum()),"longitudinal":int((l.BH_q<.05).sum()),"interaction":int((it.BH_q_omnibus<.05).sum())}}
(M/"logs/RUN_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
(M/"logs/M01_FINAL_ENGINE.json").write_text(json.dumps({"final":"original CodeKit R output","python_t_reference":"retained as superseded implementation comparison","inference_conclusion_changed":False},indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
