import os
from pathlib import Path
import csv, gzip, json, time
import numpy as np
import pandas as pd
from scipy import stats

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; OUT=RUN/"02_GSE286347/results"; LOG=RUN/"02_GSE286347/logs"
BETA=P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz"
QC=OUT/"M02_PUBLIC_METADATA_AUDIT.tsv"; OLD=P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/M02_results"
EPS=1e-6; CHUNK=2500

def bh(p):
    p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));z=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q[o]=np.minimum(z,1);return q
def ols_many(y,X,coef=1):
    inv=np.linalg.inv(X.T@X);b=(inv@X.T@y.T).T;res=y-b@X.T;df=X.shape[0]-X.shape[1];s2=np.sum(res*res,axis=1)/df;se=np.sqrt(np.maximum(s2*inv[coef,coef],0));t=b[:,coef]/se;return b[:,coef],se,t,2*stats.t.sf(abs(t),df),df

md=pd.read_csv(QC,sep="\t",dtype=str); md=md[md.cohort.isin(["PROPGER","PROPENG"])].copy();md["sex_male"]=md.sex.str.lower().eq("male").astype(int)
with gzip.open(BETA,"rt",newline="") as f: header=next(csv.reader(f))
beta={x[:-8]:x for x in header if x.endswith("_BetaVal")};md["beta_col"]=md.sample_key.map(beta)
started=time.time(); summaries=[]
for cohort in ["PROPGER","PROPENG"]:
    mm=md[md.cohort.eq(cohort)].copy(); batch=pd.get_dummies(mm.sentrix_id,prefix="slide",drop_first=True,dtype=float)
    X=np.column_stack([np.ones(len(mm)),mm.pain_group.eq("painful").astype(int),mm.sex_male.astype(int),batch.to_numpy()])
    if np.linalg.matrix_rank(X)<X.shape[1]: raise RuntimeError(f"rank-deficient batch design {cohort}")
    elig=pd.read_csv(OLD/f"M02_{cohort}_PROBE_ELIGIBILITY.tsv.gz",sep="\t",usecols=["probe_id","eligibility"]); keep=set(elig.loc[elig.eligibility.eq("eligible"),"probe_id"].astype(str))
    rows=[]; use=["ID_REF"]+list(mm.beta_col)
    for ch in pd.read_csv(BETA,usecols=use,chunksize=CHUNK):
        mask=ch.ID_REF.astype(str).isin(keep).to_numpy()
        if not mask.any(): continue
        ids=ch.loc[mask,"ID_REF"].astype(str).to_numpy();b=ch.loc[mask,list(mm.beta_col)].to_numpy(float);finite=np.all(np.isfinite(b),axis=1);ids=ids[finite];b=b[finite]
        y=np.log2(np.clip(b,EPS,1-EPS)/(1-np.clip(b,EPS,1-EPS))); eff,se,t,pv,df=ols_many(y,X,1);pain=mm.pain_group.eq("painful").to_numpy()
        rows.append(pd.DataFrame({"probe_id":ids,"M_value_difference":eff,"SE":se,"t":t,"P":pv,"df":df,"beta_difference":b[:,pain].mean(1)-b[:,~pain].mean(1),"n":len(mm),"n_painful":int(pain.sum()),"n_painless":int((~pain).sum()),"model":"sex_plus_sentrix_id"}))
    r=pd.concat(rows,ignore_index=True);r["q"]=bh(r.P);r=r.sort_values("P",kind="stable");r.to_csv(OUT/f"M02_{cohort}_EWAS_SEX_BATCH_ALL.tsv.gz",sep="\t",index=False,compression="gzip")
    summaries.append(dict(cohort=cohort,model="sex_plus_batch",status="COMPLETED_EXPLORATORY",eligible_CpGs=len(r),min_q=float(r.q.min()),q_lt_0_05=int((r.q<.05).sum()),design_rank=int(np.linalg.matrix_rank(X)),design_columns=X.shape[1],n=len(mm)))

# Fixed-region batch sensitivity in PROPENG using frozen region means.
rm=pd.read_csv(P/"analysis_v4/closeout_20260912T062015Z/03_M02_fixed_regions/PROPENG_retest_v2/REGION_MEAN_M.tsv",sep="\t").set_index("region_id")
reg=pd.read_csv(P/"analysis_v4/closeout_20260912T062015Z/03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv",sep="\t").set_index("region_id")
mm=md[md.cohort.eq("PROPENG")].set_index("sample_key");cols=[x for x in rm.columns if x in mm.index];m=mm.loc[cols];batch=pd.get_dummies(m.sentrix_id,prefix="slide",drop_first=True,dtype=float);X=np.column_stack([np.ones(len(m)),m.pain_group.eq("painful").astype(int),m.sex_male.astype(int),batch.to_numpy()]);y=rm[cols].to_numpy(float);eff,se,t,pv,df=ols_many(y,X,1)
fr=pd.DataFrame({"region_id":rm.index,"effect_M":eff,"se_M":se,"t":t,"P_two_sided":pv,"df":df,"model":"sex_plus_sentrix_id","n":len(m)});fr["BH_full_discovery_family"]=bh(fr.P_two_sided);fr=fr.join(reg[["chr","start","end","discovery_direction"]],on="region_id");fr["direction_matches_discovery"]=np.sign(fr.effect_M)==np.sign(fr.discovery_direction);fr.to_csv(OUT/"M02_FIXED_REGION_BATCH_ADJUSTED_RETEST.tsv",sep="\t",index=False)

comp=pd.read_csv(OUT/"M02_EWAS_MODEL_COMPARISON.tsv",sep="\t");comp=comp[~comp.model.eq("sex_plus_batch")];comp=pd.concat([comp,pd.DataFrame(summaries)],ignore_index=True);comp.to_csv(OUT/"M02_EWAS_MODEL_COMPARISON.tsv",sep="\t",index=False)
summary={"status":"COMPLETED_EXPLORATORY","seconds":time.time()-started,"cohorts":summaries,"fixed_regions_PROPENG":{"n":len(fr),"q_lt_0.05":int((fr.BH_full_discovery_family<.05).sum()),"min_q":float(fr.BH_full_discovery_family.min())}}
(LOG/"M02_BATCH_SENSITIVITY_RUN.json").write_text(json.dumps(summary,indent=2),encoding="utf-8");print(json.dumps(summary,indent=2))
