import os
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
ROOT=P/"analysis_v5/public_clinical_extensions_2026-09-14/01_GSE302658";R=ROOT/"results";T=P/"analysis_v5/public_clinical_extensions_2026-09-14/tests"
SYM=["npsi_total","average_pain_nrs","worst_pain_nrs","burning_superficial","pressing_deep","paroxysmal","evoked","paresthesia_dysesthesia"]
def z(x):x=np.asarray(x,float);return (x-x.mean())/x.std(ddof=1)
def hc3(y,X,j):
 X=np.asarray(X,float);y=np.asarray(y,float);inv=np.linalg.inv(X.T@X);b=inv@X.T@y;e=y-X@b;h=np.sum((X@inv)*X,1);V=inv@(X.T@((e/(1-h))[:,None]**2*X))@inv;se=np.sqrt(V[j,j]);df=len(y)-X.shape[1];p=2*stats.t.sf(abs(b[j]/se),df);return b[j],se,p,V,b
def design_base(d,p,s):
 q=d[[p,s,'age','sex','bmi']].dropna().copy();X=np.column_stack([np.ones(len(q)),z(q[p]),z(q.age),(q.sex.str.upper()=='M').astype(int),z(q.bmi)]);return q,z(q[s]),X,1
def design_long(pair,p,s):
 q=pair[[p+'_b',p+'_f',s+'_b',s+'_f','treatment_b','age_b','sex_b','bmi_b']].dropna().copy();td=pd.get_dummies(q.treatment_b,drop_first=False,dtype=float).drop(columns=['Placebo']);X=np.column_stack([np.ones(len(q)),z(q[p+'_f']-q[p+'_b']),z(q[s+'_b']),td.to_numpy(),z(q.age_b),(q.sex_b.str.upper()=='M').astype(int),z(q.bmi_b)]);return q,z(q[s+'_f']-q[s+'_b']),X,1
def design_int(pair,p,s):
 q=pair[[p+'_b',s+'_b',s+'_f','treatment_b','age_b','sex_b','bmi_b']].dropna().copy();a=(q.treatment_b=='AZD2423 20 mg').astype(int);b=(q.treatment_b=='AZD2423 150 mg').astype(int);pz=z(q[p+'_b']);X=np.column_stack([np.ones(len(q)),z(q[s+'_b']),pz,a,b,a*pz,b*pz,z(q.age_b),(q.sex_b.str.upper()=='M').astype(int),z(q.bmi_b)]);return q,z(q[s+'_f']),X,[5,6]
d=pd.read_csv(R/"GSE302658_CLINICAL_ANALYSIS_WIDE.tsv",sep='\t');base=d[d.visit.eq('baseline')];fup=d[d.visit.eq('followup')];pair=base.merge(fup,on='participant_id',suffixes=('_b','_f'))
rng=np.random.default_rng(20260914);checks=[]
for fam,file,pcol,bcol,fun in [('baseline','BASELINE_SYMPTOM_ASSOCIATIONS.tsv','P','beta_program',design_base),('longitudinal','LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv','P','beta_delta_program',design_long)]:
 r=pd.read_csv(R/file,sep='\t');r=r[r[pcol].notna()&r.status.eq('COMPLETED')];sel=r.iloc[rng.choice(len(r),min(10,len(r)),replace=False)]
 for x in sel.itertuples():
  q,y,X,j=fun(base if fam=='baseline' else pair,x.program,x.symptom);bb,se,pv,V,b=hc3(y,X,j);checks.append(dict(family=fam,program=x.program,symptom=x.symptom,n=len(q),saved_beta=getattr(x,bcol),recomputed_beta=bb,abs_beta_diff=abs(getattr(x,bcol)-bb),saved_P=x.P,recomputed_P=pv,abs_P_diff=abs(x.P-pv),pass_tolerance=abs(x.P-pv)<1e-10))
ri=pd.read_csv(R/'TREATMENT_INTERACTION_OMNIBUS.tsv',sep='\t');ri=ri[ri.P_omnibus.notna()&ri.status.eq('COMPLETED')];sel=ri.iloc[rng.choice(len(ri),min(10,len(ri)),replace=False)]
for x in sel.itertuples():
 q,y,X,idx=design_int(pair,x.program,x.symptom);_,_,_,V,b=hc3(y,X,1);bb=b[idx];VV=V[np.ix_(idx,idx)];W=float(bb@np.linalg.inv(VV)@bb);pv=float(stats.chi2.sf(W,len(idx)));checks.append(dict(family='interaction',program=x.program,symptom=x.symptom,n=len(q),saved_beta=np.nan,recomputed_beta=np.nan,abs_beta_diff=np.nan,saved_P=x.P_omnibus,recomputed_P=pv,abs_P_diff=abs(x.P_omnibus-pv),pass_tolerance=abs(x.P_omnibus-pv)<1e-10))
checks=pd.DataFrame(checks);checks.to_csv(T/'M01_INDEPENDENT_RECOMPUTATIONS.tsv',sep='\t',index=False)

# Full LOO for all primary q<0.10 longitudinal rows.
lr=pd.read_csv(R/'LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv',sep='\t');loo=[];diag=[]
for x in lr[(lr.BH_q<.10)&lr.status.eq('COMPLETED')].itertuples():
 q,y,X,j=design_long(pair,x.program,x.symptom);vals=[]
 for i in range(len(q)):
  k=np.arange(len(q))!=i;vals.append(hc3(y[k],X[k],j)[0])
 fit=sm.OLS(y,X).fit();c=fit.get_influence().cooks_distance[0]
 loo.append(dict(family='longitudinal_HC3',program=x.program,symptom=x.symptom,n=len(q),primary_beta=x.beta_delta_program,sign_preservation_fraction=float(np.mean(np.sign(vals)==np.sign(x.beta_delta_program))),loo_beta_min=float(np.min(vals)),loo_beta_max=float(np.max(vals)),max_absolute_beta_change=float(np.max(abs(np.asarray(vals)-x.beta_delta_program))),status='COMPLETED_DIAGNOSTIC'))
 diag.append(dict(family='longitudinal_HC3',program=x.program,symptom=x.symptom,n=len(q),max_cooks_D=float(c.max()),n_cooks_gt_4_over_n=int((c>4/len(q)).sum())))
pd.DataFrame(loo).to_csv(R/'LOO_STABILITY_RESULTS.tsv',sep='\t',index=False);pd.DataFrame(diag).to_csv(R/'MODEL_DIAGNOSTICS.tsv',sep='\t',index=False)
summary={'independent_checks':len(checks),'all_pass':bool(checks.pass_tolerance.all()),'max_abs_P_diff':float(checks.abs_P_diff.max()),'LOO_results':len(loo),'all_LOO_sign_preservation_min':float(pd.DataFrame(loo).sign_preservation_fraction.min()) if loo else None};(T/'M01_VERIFICATION.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
