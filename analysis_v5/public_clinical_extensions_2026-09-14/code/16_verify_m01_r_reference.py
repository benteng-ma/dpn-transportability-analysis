import os
from pathlib import Path
import json, re
import numpy as np
import pandas as pd

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; M=RUN/"01_GSE302658"; A=M/"superseded_python_t_reference_implementation"; R=M/"reference_R_run"
files=["BASELINE_SYMPTOM_ASSOCIATIONS.tsv","BASELINE_SPEARMAN_SENSITIVITY.tsv","LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv","LONGITUDINAL_SPEARMAN_SENSITIVITY.tsv","TREATMENT_INTERACTION_OMNIBUS.tsv","TREATMENT_INTERACTION_PAIRWISE.tsv"]
rows=[]
for f in files:
    a=pd.read_csv(A/f,sep="\t");b=pd.read_csv(R/f,sep="\t")
    keys=["program","symptom"]
    if "PAIRWISE" in f:
        a["dose"]=a.term.astype(str).map(lambda s:"150" if "150" in s else "20")
        b["dose"]=b.term.astype(str).map(lambda s:"150" if "150" in s else "20")
        keys.append("dose")
    z=a.merge(b,on=keys,suffixes=("_python","_R"),validate="one_to_one")
    pairs=[]
    for c in ["P","P_omnibus","BH_q","BH_q_omnibus","BH_q_pairwise"]:
        if c+"_python" in z and c+"_R" in z: pairs.append(float(np.nanmax(abs(z[c+"_python"]-z[c+"_R"]))))
    rows.append({"file":f,"rows_python":len(a),"rows_R":len(b),"matched_rows":len(z),"max_primary_P_or_q_difference":max(pairs,default=np.nan),"pass_1e_10":max(pairs,default=0)<1e-10,"final_engine":"original CodeKit R"})
out=pd.DataFrame(rows);out.to_csv(RUN/"tests/M01_R_REFERENCE_COMPARISON.tsv",sep="\t",index=False)
summary={"files":len(out),"all_pass":bool(out.pass_1e_10.all()),"max_P_or_q_difference":float(out.max_primary_P_or_q_difference.max()),"final_engine":"original CodeKit R"}
(RUN/"tests/M01_R_REFERENCE_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
