from __future__ import annotations
import csv, gzip, hashlib, json, math, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OLD=ROOT.parent/"all_extensions_2026-09-11"
BETA=ROOT/"07_resource_recovery/downloads/GSE286347_MatrixBetaVal.csv.gz"
ANNO=ROOT/"01_inputs/EPIC_V1_HG19_ANNOTATION.tsv.gz"
META=OLD/"00_admin/PUBLIC_SAMPLE_METADATA.tsv"
OUT=ROOT/"07_resource_recovery/M02_results"; OUT.mkdir(parents=True,exist_ok=True)
LOG=ROOT/"logs/05_m02_beta_ewas.json"
EPS=1e-6; CHUNK=2500

def bh(p):
    p=np.asarray(p,float); n=len(p); order=np.argsort(p); q=np.empty(n,float)
    ranked=p[order]*n/np.arange(1,n+1); ranked=np.minimum.accumulate(ranked[::-1])[::-1]
    q[order]=np.minimum(ranked,1.0); return q

def ols_many(y, X, coef=1):
    # y: probes x samples; complete finite rows only
    inv=np.linalg.inv(X.T@X); b=(inv@X.T@y.T).T
    resid=y-b@X.T; df=X.shape[0]-X.shape[1]
    s2=np.sum(resid*resid,axis=1)/df
    se=np.sqrt(np.maximum(s2*inv[coef,coef],0)); t=b[:,coef]/se
    p=2*stats.t.sf(np.abs(t),df)
    return b[:,coef],se,t,p,df

with gzip.open(BETA,"rt",newline="") as f:
    header=next(csv.reader(f))
beta_cols=[x for x in header if x.endswith("_BetaVal")]
p_cols=[x for x in header if x.endswith("_DetectionPVal")]
assert len(beta_cols)==315 and len(p_cols)==315
beta_key={x[:-8]:x for x in beta_cols}; p_key={x[:-14]:x for x in p_cols}

md=pd.read_csv(META,sep="\t",dtype=str)
md=md[md.resource.eq("GSE286347")].copy()
md["sample_key"]=md.donor_key
md["pain_group"]=md.analysis_group.replace({"painles diabetic neuropathy":"painless","painful diabetic neuropathy":"painful","control":"control"})
md["cohort"]=md.sample_key.str.extract(r"^(PROPGER|PROPENG)")[0]
md["sex_male"]=md.sex.str.lower().eq("male").astype(int)
md["beta_col"]=md.sample_key.map(beta_key); md["p_col"]=md.sample_key.map(p_key)
assert md.beta_col.notna().all() and md.p_col.notna().all() and md.sample_key.nunique()==315
dpn=md[md.pain_group.isin(["painful","painless"])].copy()
assert dpn.groupby("cohort").size().to_dict()=={"PROPENG":92,"PROPGER":139}

# Pass 1: detection-P sample QC over all assayed probes.
fail=np.zeros(len(dpn),dtype=np.int64); total=0
for ch in pd.read_csv(BETA,usecols=list(dpn.p_col),chunksize=CHUNK):
    a=ch[list(dpn.p_col)].to_numpy(float)
    fail += np.sum((~np.isfinite(a)) | (a>0.01),axis=0)
    total += a.shape[0]
dpn["probe_count"]=total; dpn["detection_fail_count"]=fail
dpn["detection_fail_fraction"]=fail/total
dpn["sample_qc_pass"]=dpn.detection_fail_fraction<=0.01
dpn.to_csv(OUT/"M02_SAMPLE_QC.tsv",sep="\t",index=False)
passed=dpn[dpn.sample_qc_pass].copy()

anno=pd.read_csv(ANNO,sep="\t",dtype={"probe_id":str,"chr":str})
anno=anno.drop_duplicates("probe_id",keep=False).set_index("probe_id")
sex_chr=set(anno.index[anno.chr.astype(str).str.replace("chr","",regex=False).isin(["X","Y"])]).intersection

results={}; exclusions=[]
usecols=["ID_REF"]+list(passed.beta_col)+list(passed.p_col)
for cohort in ["PROPGER","PROPENG"]:
    mm=passed[passed.cohort.eq(cohort)].copy()
    X=np.column_stack([np.ones(len(mm)),mm.pain_group.eq("painful").astype(int),mm.sex_male.astype(int)])
    rows=[]; ex=[]
    for ch in pd.read_csv(BETA,usecols=usecols,chunksize=CHUNK):
        ids=ch.ID_REF.astype(str).to_numpy()
        b=ch[list(mm.beta_col)].to_numpy(float)
        d=ch[list(mm.p_col)].to_numpy(float)
        failfrac=np.mean((~np.isfinite(d))|(d>0.01),axis=1)
        finite=np.all(np.isfinite(b),axis=1); bounds=np.all((b>=0)&(b<=1),axis=1)
        annotated=np.fromiter((x in anno.index for x in ids),bool,len(ids))
        autosome=np.fromiter((((x in anno.index) and str(anno.at[x,"chr"]).replace("chr","") not in ("X","Y")) for x in ids),bool,len(ids))
        keep=(failfrac<=0.05)&finite&bounds&annotated&autosome
        reason=np.full(len(ids),"eligible",object)
        reason[failfrac>0.05]="detection_failure"
        reason[(failfrac<=0.05)&~finite]="nonfinite_beta"
        reason[(failfrac<=0.05)&finite&~bounds]="beta_out_of_bounds"
        reason[(failfrac<=0.05)&finite&bounds&~annotated]="annotation_missing"
        reason[(failfrac<=0.05)&finite&bounds&annotated&~autosome]="sex_chromosome"
        ex.append(pd.DataFrame({"probe_id":ids,"cohort":cohort,"detection_fail_fraction":failfrac,"eligibility":reason}))
        if keep.any():
            bk=np.clip(b[keep],EPS,1-EPS); y=np.log2(bk/(1-bk))
            eff,se,t,p,df=ols_many(y,X,1)
            painful=mm.pain_group.eq("painful").to_numpy(); painless=~painful
            bd=b[keep][:,painful].mean(1)-b[keep][:,painless].mean(1)
            rows.append(pd.DataFrame({"probe_id":ids[keep],"M_value_difference":eff,"SE":se,"t":t,"P":p,"df":df,"beta_difference":bd,"n":len(mm),"n_painful":int(painful.sum()),"n_painless":int(painless.sum())}))
    rr=pd.concat(rows,ignore_index=True); rr["q"]=bh(rr.P.to_numpy())
    rr=rr.sort_values("P",kind="stable")
    rr.to_csv(OUT/f"M02_{cohort}_EWAS_ALL_ELIGIBLE.tsv.gz",sep="\t",index=False,compression="gzip")
    pd.concat(ex,ignore_index=True).to_csv(OUT/f"M02_{cohort}_PROBE_ELIGIBILITY.tsv.gz",sep="\t",index=False,compression="gzip")
    results[cohort]=rr

g=results["PROPGER"].set_index("probe_id"); e=results["PROPENG"].set_index("probe_id")
j=g[["M_value_difference","P","q","beta_difference"]].join(e[["M_value_difference","P","q","beta_difference"]],lsuffix="_PROPGER",rsuffix="_PROPENG",how="inner")
j["direction_concordant"]=np.sign(j.M_value_difference_PROPGER)==np.sign(j.M_value_difference_PROPENG)
j["PROPGER_q05"]=j.q_PROPGER<0.05; j["PROPENG_q05"]=j.q_PROPENG<0.05
j.to_csv(OUT/"M02_CROSS_COHORT_PROBE_SUPPORT.tsv.gz",sep="\t",compression="gzip")

summary=[]
for c,r in results.items():
    summary.append({"cohort":c,"samples":int(r.n.iloc[0]),"eligible_probes":len(r),"q_lt_005":int((r.q<.05).sum()),"q_lt_010":int((r.q<.10).sum()),"minimum_q":float(r.q.min())})
summary.append({"cohort":"cross_cohort","samples":len(passed),"eligible_probes":len(j),"q_lt_005":int(((j.q_PROPGER<.05)&(j.q_PROPENG<.05)&j.direction_concordant).sum()),"q_lt_010":int(((j.q_PROPGER<.10)&(j.q_PROPENG<.10)&j.direction_concordant).sum()),"minimum_q":float("nan")})
pd.DataFrame(summary).to_csv(OUT/"M02_EWAS_SUMMARY.tsv",sep="\t",index=False)
run={"status":"COMPLETED_EXPLORATORY_BETA_BRANCH","matrix_sha256":"6c3b6a7f4481463d5cb9c7f320734582c2c284917e10e5a54ec8e7eb7dd5f0f4","all_samples":315,"DPN_samples":len(dpn),"DPN_samples_pass":len(passed),"probe_rows":total,"cohort_summaries":summary,"primary_age_batch_adjusted_model":"BLOCKED_COVARIATES"}
LOG.write_text(json.dumps(run,indent=2),encoding="utf-8")
print(json.dumps(run,indent=2))
