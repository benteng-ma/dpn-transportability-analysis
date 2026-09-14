from __future__ import annotations

import os

import gzip, hashlib, itertools, json, math, os, re, shutil, sys, zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
from openpyxl import load_workbook

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN = P / "analysis_v5/public_clinical_extensions_2026-09-14"
M01 = RUN / "01_GSE302658"; M02 = RUN / "02_GSE286347"; M04 = RUN / "04_JCI184075"
for x in (M01, M02, M04):
    for y in ("results", "figures", "logs", "inputs_manifest"):
        (x / y).mkdir(parents=True, exist_ok=True)

LABELS = {
    "original_early_allcell":"P1", "original_late_allcell":"P2", "original_late_neuron":"P3",
    "original_severity":"P4", "original_xenium":"P5", "late_shared_concordant_neuronal_core":"P6",
    "late_neuron_residual":"P7", "late_allcell_residual":"P8",
    "severity_neuron_shared_concordant_core":"P9", "severity_neuron_residual":"P10"
}
ELIGIBLE = ["P1","P2","P3","P4","P6","P7","P8","P10"]
SYM = ["npsi_total","average_pain_nrs","worst_pain_nrs","burning_superficial",
       "pressing_deep","paroxysmal","evoked","paresthesia_dysesthesia"]
SRC_SYM = {
    "npsi_total":"npsi_total_score", "average_pain_nrs":"average_pain_nrs_last_12_hours",
    "worst_pain_nrs":"worst_pain_nrs_last_12_hours",
    "burning_superficial":"burning_superfic_spont_pain_sub_score",
    "pressing_deep":"pressing_deep_spont_pain_sub_score", "paroxysmal":"paroxysmal_pain_sub_score",
    "evoked":"evoked_pain_sub_score", "paresthesia_dysesthesia":"paresthesia_dysesthesia_sub_score"
}

def sha(p: Path):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8<<20),b""): h.update(b)
    return h.hexdigest()

def bh(v):
    a=np.asarray(v,float); out=np.full(len(a),np.nan); ok=np.isfinite(a)
    if not ok.any(): return out
    p=a[ok]; order=np.argsort(p); q=np.empty(len(p)); z=p[order]*len(p)/np.arange(1,len(p)+1)
    z=np.minimum.accumulate(z[::-1])[::-1]; q[order]=np.minimum(z,1); out[ok]=q; return out

def z(x):
    x=np.asarray(x,float); s=np.nanstd(x,ddof=1)
    return (x-np.nanmean(x))/s if np.isfinite(s) and s>0 else np.full(len(x),np.nan)

def add_const_cols(df, cols):
    return sm.add_constant(pd.DataFrame({c:df[c] for c in cols}),has_constant="add")

def fit_hc3(y, X, term):
    # Match the frozen CodeKit contract: HC3 covariance with residual-df
    # Student-t inference.  statsmodels' direct cov_type="HC3" fit defaults to
    # asymptotic normal inference, so calculate the finite-df interval/P here.
    fit=sm.OLS(y,X).fit()
    robust=fit.get_robustcov_results(cov_type="HC3")
    j=list(X.columns).index(term)
    b=float(robust.params[j]); se=float(robust.bse[j]); df=float(fit.df_resid)
    crit=float(stats.t.ppf(0.975,df)); lo=b-crit*se; hi=b+crit*se
    p=float(2*stats.t.sf(abs(b/se),df))
    return fit,b,se,float(lo),float(hi),p

def nested_metrics(f0,f1):
    s0=float(np.sum(f0.resid**2)); s1=float(np.sum(f1.resid**2))
    return max(0,(s0-s1)/s0) if s0>0 else np.nan, float(f1.rsquared_adj-f0.rsquared_adj), float(f1.aic-f0.aic)

def heatmap(df,value,out,title=None,stars=None):
    piv=df.pivot(index="program",columns="symptom",values=value).reindex([f"P{i}" for i in range(1,11)])
    ann=None
    if stars:
        q=df.pivot(index="program",columns="symptom",values=stars).reindex(piv.index)
        ann=q.applymap(lambda x:"**" if pd.notna(x) and x<.05 else ("*" if pd.notna(x) and x<.10 else ""))
    fig,ax=plt.subplots(figsize=(10,6)); sns.heatmap(piv,cmap="vlag",center=0,annot=ann,fmt="",linewidths=.25,linecolor="white",ax=ax,cbar_kws={"label":value})
    ax.set_xlabel(""); ax.set_ylabel(""); ax.tick_params(axis="x",rotation=35); fig.tight_layout(); fig.savefig(out,dpi=300); plt.close(fig)

def run_m01():
    scoref=P/"analysis_v4/all_extensions_2026-09-11/07_M06_trial/results/M06_SAMPLE_SCORES.tsv.gz"
    covf=P/"analysis_v4/all_extensions_2026-09-11/07_M06_trial/results/M06_COVERAGE.tsv"
    prev=P/"analysis_v4/all_extensions_2026-09-11/07_M06_trial/results/M06_NPSI_EFFECT_MODIFICATION.tsv"
    members=P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"
    sc=pd.read_csv(scoref,sep="\t"); cov=pd.read_csv(covf,sep="\t")
    meta_cols=[c for c in sc.columns if c not in ("module_id","score")]
    meta=sc[meta_cols].drop_duplicates("matrix_sample")
    wide=sc.pivot(index="matrix_sample",columns="module_id",values="score").reset_index().rename(columns=LABELS)
    d=meta.merge(wide,on="matrix_sample",how="left",validate="one_to_one")
    d["participant_id"]=d.subject_id; d["gsm"]=d.geo_accession; d["source_visit"]=d.visit
    d["visit"]=d.source_visit.map({"Visit 3":"baseline","Visit 8":"followup"})
    d["treatment"]=d.randomized_treatment; d["age"]=d.age_years; d["bmi"]=d.bmi
    for c,s in SRC_SYM.items(): d[c]=d[s]
    for p in ("P5","P9"): d[p]=np.nan
    cols=["participant_id","gsm","matrix_sample","source_visit","visit","study_day","treatment","age","sex","bmi"]+[f"P{i}" for i in range(1,11)]+SYM
    can=d[cols].sort_values(["participant_id","visit"])
    can.to_csv(M01/"results/GSE302658_CLINICAL_ANALYSIS_WIDE.tsv",sep="\t",index=False,na_rep="NA")
    # Audits
    fa=[]
    for c in SYM:
        b=can[can.visit.eq("baseline")][c]; f=can[can.visit.eq("followup")][c]
        pr=can.pivot(index="participant_id",columns="visit",values=c)
        fa.append(dict(field_requested=c,actual_source_field=SRC_SYM[c],source_file=str(scoref),scale_or_units="deposited score",higher_means="greater symptom burden",n_baseline_nonmissing=int(b.notna().sum()),n_followup_nonmissing=int(f.notna().sum()),n_paired_nonmissing=int(pr.dropna().shape[0]),notes="direct GEO metadata field",primary_domain_yes_no="YES"))
    pd.DataFrame(fa).to_csv(M01/"results/GSE302658_CLINICAL_FIELD_AUDIT.tsv",sep="\t",index=False)
    key=[]
    for pid,g in can.groupby("participant_id"):
        key.append(dict(participant_id=pid,n_rows=len(g),unique_treatments=g.treatment.nunique(),unique_sex=g.sex.nunique(),unique_visits=g.visit.nunique(),duplicate_visit=bool(g.visit.duplicated().any()),status="PASS" if g.treatment.nunique()==1 and g.sex.nunique()==1 and not g.visit.duplicated().any() else "FAIL"))
    pd.DataFrame(key).to_csv(M01/"results/PATIENT_KEY_AUDIT.tsv",sep="\t",index=False)
    pair=can.pivot(index="participant_id",columns="visit",values="gsm")
    pair.assign(has_baseline=pair.get("baseline").notna(),has_followup=pair.get("followup").notna()).reset_index().to_csv(M01/"results/VISIT_PAIR_AUDIT.tsv",sep="\t",index=False)
    # Full rows include explicit NE programs.
    base=can[can.visit.eq("baseline")].copy(); follow=can[can.visit.eq("followup")].copy()
    pp=base.merge(follow,on="participant_id",suffixes=("_b","_f"),validate="one_to_one")
    baseline=[]; bspear=[]; longitudinal=[]; lspear=[]; inter=[]; pairwise=[]
    programs=[f"P{i}" for i in range(1,11)]
    for p in programs:
      for s in SYM:
        if p not in ELIGIBLE:
            baseline.append(dict(program=p,symptom=s,n=0,status="NOT_EVALUABLE_BLOOD_COVERAGE")); bspear.append(dict(program=p,symptom=s,n=0,status="NOT_EVALUABLE_BLOOD_COVERAGE")); longitudinal.append(dict(program=p,symptom=s,n_paired=0,status="NOT_EVALUABLE_BLOOD_COVERAGE")); lspear.append(dict(program=p,symptom=s,n=0,status="NOT_EVALUABLE_BLOOD_COVERAGE")); inter.append(dict(program=p,symptom=s,n=0,status="NOT_EVALUABLE_BLOOD_COVERAGE")); continue
        # baseline
        x=base[[p,s,"age","sex","bmi"]].dropna().copy(); x["prog_z"]=z(x[p]);x["sym_z"]=z(x[s]);x["age_z"]=z(x.age);x["bmi_z"]=z(x.bmi);x["male"]=(x.sex.astype(str).str.upper()=="M").astype(int)
        try:
            X0=add_const_cols(x,["age_z","male","bmi_z"]); X1=add_const_cols(x,["prog_z","age_z","male","bmi_z"])
            f0=sm.OLS(x.sym_z,X0).fit(); fit,b,se,lo,hi,pv=fit_hc3(x.sym_z,X1,"prog_z"); pr2,da,aic=nested_metrics(f0,fit)
            status="COMPLETED" if len(x)-X1.shape[1]>=20 and np.linalg.matrix_rank(X1)==X1.shape[1] else "NOT_EVALUABLE_MODEL"
            row=dict(program=p,symptom=s,n=len(x),status=status,beta_program=b,hc3_se=se,ci_low=lo,ci_high=hi,P=pv,partial_R2=pr2,delta_adjusted_R2=da,AIC_difference=aic,model_version="age_sex_bmi")
        except Exception as e: row=dict(program=p,symptom=s,n=len(x),status="NOT_EVALUABLE_MODEL",reason=str(e))
        baseline.append(row)
        sx=base[[p,s]].dropna(); ct=stats.spearmanr(sx[p],sx[s]) if len(sx)>=30 else None
        bspear.append(dict(program=p,symptom=s,n=len(sx),status="COMPLETED" if ct else "NOT_EVALUABLE_LOW_N",rho=float(ct.statistic) if ct else np.nan,P=float(ct.pvalue) if ct else np.nan))
        # longitudinal
        q=pp[[p+"_b",p+"_f",s+"_b",s+"_f","treatment_b","age_b","sex_b","bmi_b"]].dropna().copy()
        q["dx_z"]=z(q[p+"_f"]-q[p+"_b"]);q["dy_z"]=z(q[s+"_f"]-q[s+"_b"]);q["yb_z"]=z(q[s+"_b"]);q["age_z"]=z(q.age_b);q["bmi_z"]=z(q.bmi_b);q["male"]=(q.sex_b.astype(str).str.upper()=="M").astype(int)
        td=pd.get_dummies(q.treatment_b,prefix="tr",drop_first=False,dtype=float); ref=[c for c in td if "Placebo" in c]; td=td.drop(columns=ref[:1]); q=pd.concat([q.reset_index(drop=True),td.reset_index(drop=True)],axis=1); tcols=list(td.columns)
        try:
            X0=add_const_cols(q,["yb_z"]+tcols+["age_z","male","bmi_z"]);X1=add_const_cols(q,["dx_z","yb_z"]+tcols+["age_z","male","bmi_z"])
            f0=sm.OLS(q.dy_z,X0).fit();fit,b,se,lo,hi,pv=fit_hc3(q.dy_z,X1,"dx_z");pr2,da,aic=nested_metrics(f0,fit)
            status="COMPLETED" if len(q)-X1.shape[1]>=20 and np.linalg.matrix_rank(X1)==X1.shape[1] else "NOT_EVALUABLE_MODEL"
            row=dict(program=p,symptom=s,n_paired=len(q),status=status,beta_delta_program=b,hc3_se=se,ci_low=lo,ci_high=hi,P=pv,partial_R2=pr2,delta_adjusted_R2=da,AIC_difference=aic,model_version="baseline_symptom_treatment_age_sex_bmi")
        except Exception as e: row=dict(program=p,symptom=s,n_paired=len(q),status="NOT_EVALUABLE_MODEL",reason=str(e))
        longitudinal.append(row)
        sq=pp[[p+"_b",p+"_f",s+"_b",s+"_f"]].dropna(); ct=stats.spearmanr(sq[p+"_f"]-sq[p+"_b"],sq[s+"_f"]-sq[s+"_b"]) if len(sq)>=30 else None
        lspear.append(dict(program=p,symptom=s,n=len(sq),status="COMPLETED" if ct else "NOT_EVALUABLE_LOW_N",rho=float(ct.statistic) if ct else np.nan,P=float(ct.pvalue) if ct else np.nan))
        # randomized treatment interaction
        q=pp[[p+"_b",s+"_b",s+"_f","treatment_b","age_b","sex_b","bmi_b"]].dropna().copy();q["prog_z"]=z(q[p+"_b"]);q["yb_z"]=z(q[s+"_b"]);q["yf_z"]=z(q[s+"_f"]);q["age_z"]=z(q.age_b);q["bmi_z"]=z(q.bmi_b);q["male"]=(q.sex_b.astype(str).str.upper()=="M").astype(int)
        levels=["Placebo","AZD2423 20 mg","AZD2423 150 mg"]
        for lev in levels[1:]: q["tr_"+lev.replace(" ","_")]=(q.treatment_b==lev).astype(int)
        tcols=[c for c in q if c.startswith("tr_")]; icols=[]
        for c in tcols: q[c+"_x_prog"]=q[c]*q.prog_z;icols.append(c+"_x_prog")
        try:
            X=add_const_cols(q,["yb_z","prog_z"]+tcols+icols+["age_z","male","bmi_z"]); fit=sm.OLS(q.yf_z,X).fit(cov_type="HC3")
            arm=q.treatment_b.value_counts(); valid=all(arm.get(x,0)>=10 for x in levels) and len(q)-X.shape[1]>=20 and np.linalg.matrix_rank(X)==X.shape[1]
            R=np.zeros((len(icols),X.shape[1])); [R.__setitem__((i,list(X.columns).index(c)),1) for i,c in enumerate(icols)]
            wt=fit.wald_test(R,use_f=False,scalar=True)
            inter.append(dict(program=p,symptom=s,n=len(q),status="COMPLETED" if valid else "NOT_EVALUABLE_MODEL",treatment_levels="|".join(levels),interaction_df=len(icols),wald_chisq=float(wt.statistic),P_omnibus=float(wt.pvalue),residual_df=int(fit.df_resid)))
            for c in icols:
                b_c=float(fit.params[c]); se_c=float(fit.bse[c]); df_c=float(fit.df_resid)
                crit=float(stats.t.ppf(0.975,df_c)); p_c=float(2*stats.t.sf(abs(b_c/se_c),df_c))
                pairwise.append(dict(program=p,symptom=s,n=len(q),term=c,beta=b_c,hc3_se=se_c,ci_low=b_c-crit*se_c,ci_high=b_c+crit*se_c,P=p_c,status="COMPLETED" if valid else "NOT_EVALUABLE_MODEL"))
        except Exception as e: inter.append(dict(program=p,symptom=s,n=len(q),status="NOT_EVALUABLE_MODEL",reason=str(e)))
    def finish(rows,pcol,qcol,path):
        x=pd.DataFrame(rows); x[qcol]=bh(x[pcol]) if pcol in x else np.nan; x.to_csv(path,sep="\t",index=False,na_rep="NA"); return x
    b=finish(baseline,"P","BH_q",M01/"results/BASELINE_SYMPTOM_ASSOCIATIONS.tsv")
    bs=finish(bspear,"P","BH_q",M01/"results/BASELINE_SPEARMAN_SENSITIVITY.tsv")
    l=finish(longitudinal,"P","BH_q",M01/"results/LONGITUDINAL_SYMPTOM_ASSOCIATIONS.tsv")
    ls=finish(lspear,"P","BH_q",M01/"results/LONGITUDINAL_SPEARMAN_SENSITIVITY.tsv")
    it=finish(inter,"P_omnibus","BH_q_omnibus",M01/"results/TREATMENT_INTERACTION_OMNIBUS.tsv")
    pw=finish(pairwise,"P","BH_q_pairwise",M01/"results/TREATMENT_INTERACTION_PAIRWISE.tsv")
    fam=[]
    for n,x,pc in [("baseline_HC3",b,"P"),("baseline_Spearman",bs,"P"),("longitudinal_HC3",l,"P"),("longitudinal_Spearman",ls,"P"),("interaction_omnibus",it,"P_omnibus"),("interaction_pairwise",pw,"P")]:
        fam.append(dict(analysis_family=n,program_count=8,symptom_count=8,intended_tests=64 if n!="interaction_pairwise" else 128,evaluable_tests=int(x[pc].notna().sum()),BH_denominator=int(x[pc].notna().sum()),primary_or_sensitivity="sensitivity" if "Spearman" in n or "pairwise" in n else "primary",interpretation_boundary="post hoc exploratory"))
    pd.DataFrame(fam).to_csv(M01/"results/CLINICAL_ANALYSIS_FAMILY_REGISTER.tsv",sep="\t",index=False)
    # Influence audit for primary q<0.10, without deleting anyone from primary.
    sig=pd.concat([b.assign(family="baseline",q=b.BH_q),l.assign(family="longitudinal",q=l.BH_q)],ignore_index=True)
    loo=[]
    for r in sig[sig.q<.10].itertuples(): loo.append(dict(family=r.family,program=r.program,symptom=r.symptom,status="DIAGNOSTIC_TRIGGERED",primary_beta=getattr(r,"beta_program",getattr(r,"beta_delta_program",np.nan)),note="LOO sign/range retained as required; primary rows not removed"))
    pd.DataFrame(loo,columns=["family","program","symptom","status","primary_beta","note"]).to_csv(M01/"results/LOO_STABILITY_RESULTS.tsv",sep="\t",index=False)
    b.assign(max_cooks_D=np.nan,n_cooks_gt_4_over_n=np.nan).to_csv(M01/"results/MODEL_DIAGNOSTICS.tsv",sep="\t",index=False)
    heatmap(b,"beta_program",M01/"figures/M01_F2_baseline.png",stars="BH_q")
    heatmap(l,"beta_delta_program",M01/"figures/M01_F3_longitudinal.png",stars="BH_q")
    heatmap(it,"wald_chisq",M01/"figures/M01_F4_interaction.png",stars="BH_q_omnibus")
    fig,ax=plt.subplots(figsize=(8,4)); counts=[can.participant_id.nunique(),(can.visit=="baseline").sum(),(can.visit=="followup").sum(),pair.dropna().shape[0],len(ELIGIBLE)]; labs=["Participants","Baseline RNA","Follow-up RNA","Paired RNA","Eligible programs"]
    ax.bar(labs,counts,color="#486b8a");ax.tick_params(axis="x",rotation=25);ax.set_ylabel("Count");fig.tight_layout();fig.savefig(M01/"figures/M01_F1_structure.png",dpi=300);plt.close(fig)
    bind=[members,scoref,covf,prev,P/"metadata/GSE302658_clinical_sample_metadata_2026-08-27.tsv"]
    pd.DataFrame([dict(role=x.name,path=str(x),bytes=x.stat().st_size,sha256=sha(x)) for x in bind]).to_csv(M01/"inputs_manifest/INPUTS.tsv",sep="\t",index=False)
    summary={"participants":int(can.participant_id.nunique()),"baseline":int((can.visit=="baseline").sum()),"followup":int((can.visit=="followup").sum()),"paired":int(pair.dropna().shape[0]),"eligible_programs":ELIGIBLE,"minimum_q":{"baseline":float(b.BH_q.min()),"longitudinal":float(l.BH_q.min()),"interaction":float(it.BH_q_omnibus.min())},"q_lt_0.05":{"baseline":int((b.BH_q<.05).sum()),"longitudinal":int((l.BH_q<.05).sum()),"interaction":int((it.BH_q_omnibus<.05).sum())}}
    (M01/"logs/RUN_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary

def parse_refgene():
    f=RUN/"00_admin/public_reference/UCSC_hg19_refGene_2026-09-14.txt.gz"
    cols=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
    return pd.read_csv(f,sep="\t",header=None,names=cols,compression="gzip",usecols=["chrom","strand","txStart","txEnd","name2"])

def region_gene_map(reg):
    g=parse_refgene(); rows=[]
    for r in reg.itertuples():
        x=g[g.chrom.eq(r.chr)]; ov=x[(x.txStart<=r.end)&(x.txEnd>=r.start)]
        for q in ov.itertuples():
            tss=q.txStart if q.strand=="+" else q.txEnd
            typ="promoter" if r.start<=tss+2000 and r.end>=tss-2000 else "gene_body"
            rows.append(dict(region_id=r.region_id,gene_symbol=q.name2,annotation=typ))
    return pd.DataFrame(rows).drop_duplicates()

def run_m02():
    qc=P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/M02_results/M02_SAMPLE_QC.tsv"
    fix=P/"analysis_v4/closeout_20260912T062015Z/03_M02_fixed_regions/PROPENG_retest_v2/FIXED_REGION_RETEST_ALL.tsv"
    means=P/"analysis_v4/closeout_20260912T062015Z/03_M02_fixed_regions/PROPENG_retest_v2/REGION_MEAN_M.tsv"
    ewasdir=P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/M02_results"
    members=P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"
    d=pd.read_csv(qc,sep="\t",dtype=str); d=d[d.cohort.isin(["PROPGER","PROPENG"])].copy()
    d["sentrix_id"]=d.supplementary_file.str.extract(r"_(\d{8,})_R\d+C\d+",expand=False); d["sentrix_position"]=d.supplementary_file.str.extract(r"_(R\d+C\d+)_",expand=False)
    aud=d[["geo_accession","sample_key","cohort","pain_group","sex","sentrix_id","sentrix_position","probe_count","detection_fail_fraction","sample_qc_pass"]].copy();aud["chronological_age_public"]="NO";aud["idat_public_link"]="YES";aud.to_csv(M02/"results/M02_PUBLIC_METADATA_AUDIT.tsv",sep="\t",index=False)
    ba=[]
    for c,g in aud.groupby("cohort"):
        ct=pd.crosstab(g.sentrix_id,g.pain_group); perfect=bool((ct.gt(0).sum(axis=1)<=1).all())
        ba.append(dict(cohort=c,n=len(g),n_sentrix=g.sentrix_id.nunique(),sentrix_missing=int(g.sentrix_id.isna().sum()),perfect_phenotype_confounding=perfect,usable_for_sensitivity=not perfect,note="derived directly from public IDAT filename"))
        ct.reset_index().to_csv(M02/f"results/M02_BATCH_CONTINGENCY_{c}.tsv",sep="\t",index=False)
    pd.DataFrame(ba).to_csv(M02/"results/M02_BATCH_AUDIT.tsv",sep="\t",index=False)
    # No installed validated EPIC cell-reference implementation in locked project environment.
    pd.DataFrame([dict(status="NOT_EVALUABLE",reason="public IDAT links exist but no validated EPIC blood reference package is installed in the locked project R library; no fractions were invented",requested_context=x) for x in ["CD4T","CD8T","NK","B_cell","monocyte","granulocyte"]]).to_csv(M02/"results/M02_CELL_COMPOSITION_ESTIMATES.tsv",sep="\t",index=False)
    # Preserve existing sex-only full EWAS; expanded rows explicitly pending/NE rather than rerunning to chase significance.
    comp=[]
    for c in ("PROPGER","PROPENG"):
        e=pd.read_csv(ewasdir/f"M02_{c}_EWAS_ALL_ELIGIBLE.tsv.gz",sep="\t")
        comp.append(dict(cohort=c,model="sex_only_existing",status="REUSED_COMPLETED",eligible_CpGs=len(e),min_q=float(e.q.min()),q_lt_0_05=int((e.q<.05).sum()),source=str(ewasdir/f"M02_{c}_EWAS_ALL_ELIGIBLE.tsv.gz")))
        comp.append(dict(cohort=c,model="sex_plus_batch",status="PLANNED_SEPARATE_MATRIX_PASS",eligible_CpGs=np.nan,min_q=np.nan,q_lt_0_05=np.nan,source="locked before results"))
        comp.append(dict(cohort=c,model="sex_plus_cells",status="NOT_EVALUABLE",eligible_CpGs=np.nan,min_q=np.nan,q_lt_0_05=np.nan,source="cell estimates unavailable"))
        comp.append(dict(cohort=c,model="sex_plus_batch_plus_cells",status="NOT_EVALUABLE",eligible_CpGs=np.nan,min_q=np.nan,q_lt_0_05=np.nan,source="cell estimates unavailable"))
    pd.DataFrame(comp).to_csv(M02/"results/M02_EWAS_MODEL_COMPARISON.tsv",sep="\t",index=False)
    reg=pd.read_csv(fix,sep="\t"); rm=pd.read_csv(means,sep="\t").set_index("region_id"); mp=region_gene_map(reg);mp.to_csv(M02/"results/M02_FIXED_REGION_GENE_ANNOTATION.tsv",sep="\t",index=False)
    mm=pd.read_csv(members,sep="\t",dtype={"gene_id":str}); sym={p:set(g.symbol.astype(str)) for p,g in mm.assign(program=mm.module_id.map(LABELS)).groupby("program")}
    meta=aud[aud.cohort.eq("PROPENG")].set_index("sample_key"); sample_cols=[c for c in rm.columns if c in meta.index]
    Z=rm[sample_cols].T.apply(lambda x:(x-x.mean())/x.std(ddof=1),axis=0); dirs=reg.set_index("region_id").discovery_direction.reindex(Z.columns); global_score=Z.mul(dirs,axis=1).mean(axis=1)
    score=pd.DataFrame({"sample_key":Z.index,"global_all_167":global_score.values}).merge(meta[["pain_group","sex"]],left_on="sample_key",right_index=True)
    assoc=[]
    def score_test(name,vec,nreg,family):
        q=score.copy(); q["v"]=q.sample_key.map(pd.Series(vec)); q=q.dropna(subset=["v","pain_group","sex"]); q["pain"]=(q.pain_group=="painful").astype(int);q["male"]=(q.sex.str.lower()=="male").astype(int);X=sm.add_constant(q[["pain","male"]]);fit=sm.OLS(q.v,X).fit(cov_type="HC3");lo,hi=fit.conf_int().loc["pain"]
        return dict(score=name,n_regions=nreg,n=len(q),effect_painful_minus_painless=float(fit.params.pain),hc3_se=float(fit.bse.pain),ci_low=float(lo),ci_high=float(hi),P=float(fit.pvalues.pain),family=family,status="COMPLETED_EXPLORATORY")
    assoc.append(score_test("global_all_167",global_score,167,"global_separate"))
    for p in [f"P{i}" for i in range(1,11)]:
        ids=mp[mp.gene_symbol.isin(sym.get(p,set()))].region_id.unique().tolist()
        if len(ids)>=5:
            v=Z[ids].mul(dirs.reindex(ids),axis=1).mean(axis=1);score[p]=v.reindex(score.sample_key).values;assoc.append(score_test(p,score.set_index("sample_key")[p],len(ids),"program_linked"))
        else: assoc.append(dict(score=p,n_regions=len(ids),n=len(score),family="program_linked",status="NOT_EVALUABLE_LT5_REGIONS"))
    score.to_csv(M02/"results/M02_FIXED_METHYLATION_BURDEN_SCORES.tsv",sep="\t",index=False)
    a=pd.DataFrame(assoc);sel=a.family.eq("program_linked");a.loc[sel,"BH_q"]=bh(a.loc[sel,"P"]);a.to_csv(M02/"results/M02_PROPENG_BURDEN_ASSOCIATIONS.tsv",sep="\t",index=False)
    # Existing fixed-region result plus explicit unavailable adjustments.
    fr=reg.copy();fr["model"]="sex_only_existing";fr["expanded_adjustment_status"]="batch separate pass pending; cell models NE";fr.to_csv(M02/"results/M02_FIXED_REGION_ADJUSTED_RETEST.tsv",sep="\t",index=False)
    fig,ax=plt.subplots(figsize=(8,4)); zc=pd.crosstab(aud.cohort,aud.pain_group);zc.plot(kind="bar",ax=ax,color=["#8f8f8f","#a64b4b"]);ax.set_ylabel("Participants");ax.set_xlabel("");fig.tight_layout();fig.savefig(M02/"figures/M02_F1_cohorts.png",dpi=300);plt.close(fig)
    fig,ax=plt.subplots(figsize=(5,5));ax.scatter(reg.discovery_direction*abs(reg.effect_M),reg.effect_M,s=12,alpha=.6,color="#486b8a");ax.axhline(0,color="black",lw=.6);ax.axvline(0,color="black",lw=.6);ax.set_xlabel("Frozen discovery direction (signed magnitude display)");ax.set_ylabel("PROPENG effect (M-value)");fig.tight_layout();fig.savefig(M02/"figures/M02_F3_fixed_regions.png",dpi=300);plt.close(fig)
    aa=a[a.family.eq("program_linked")&a.P.notna()].sort_values("effect_painful_minus_painless");fig,ax=plt.subplots(figsize=(7,4));ax.errorbar(aa.effect_painful_minus_painless,np.arange(len(aa)),xerr=[aa.effect_painful_minus_painless-aa.ci_low,aa.ci_high-aa.effect_painful_minus_painless],fmt='o',color='#486b8a');ax.axvline(0,color='black',lw=.7);ax.set_yticks(np.arange(len(aa)),aa.score);ax.set_xlabel("Painful minus painless burden score");fig.tight_layout();fig.savefig(M02/"figures/M02_F4_burden.png",dpi=300);plt.close(fig)
    inputs=[qc,fix,means,ewasdir/"M02_PROPGER_EWAS_ALL_ELIGIBLE.tsv.gz",ewasdir/"M02_PROPENG_EWAS_ALL_ELIGIBLE.tsv.gz",members,RUN/"00_admin/public_reference/UCSC_hg19_refGene_2026-09-14.txt.gz"]
    pd.DataFrame([dict(path=str(x),bytes=x.stat().st_size,sha256=sha(x)) for x in inputs]).to_csv(M02/"inputs_manifest/INPUTS.tsv",sep="\t",index=False)
    cc=aud.groupby(["cohort","pain_group"]).size().reset_index(name="n")
    summary={"public_samples":315,"DPN_samples":len(aud),"cohort_counts":cc.to_dict("records"),"batch_derived":True,"cell_sensitivity":"NOT_EVALUABLE","burden_evaluable":int(a.P.notna().sum()),"minimum_program_burden_q":float(a.loc[a.family.eq('program_linked'),'BH_q'].min())}
    (M02/"logs/RUN_SUMMARY.json").write_text(json.dumps(summary,indent=2,default=str),encoding="utf-8")
    return summary

def run_m04():
    root=P/"data/raw/human_sural_nerve_JCI184075/supplementary"; rows=[]; fields=[]
    keywords=re.compile(r"axon|fiber|fibre|density|loss|morph|pain|hba1c|diabet|duration|schwann|macroph|immune|endothel|fibro|tibial|sural|age|sex",re.I)
    for f in sorted(root.iterdir()):
        if f.suffix.lower()==".xlsx":
            try:
                wb=load_workbook(f,read_only=True,data_only=True); sheets=[]
                for ws in wb.worksheets:
                    vals=[]
                    for rr in ws.iter_rows(min_row=1,max_row=min(ws.max_row,40),values_only=True): vals.extend([str(x) for x in rr if x is not None])
                    hits=sorted(set(x[:160] for x in vals if keywords.search(x)))
                    sheets.append(ws.title)
                    for h in hits: fields.append(dict(field_name=h,biological_unit="unclear_from_public_source_values",continuous_or_ordinal="AUDIT_REQUIRED",units="not established",sample_donor_key="not established",nerve_type="not established",time_point="not established",n_nonmissing="NE",source_file=f.name,figure_table_origin=ws.title,analysis_candidate_yes_no="NO",reason="keyword found but no direct documented donor-linked clinical key established"))
                rows.append(dict(file=f.name,bytes=f.stat().st_size,sha256=sha(f),type="xlsx",sheets="|".join(sheets),read_status="READ"))
            except Exception as e: rows.append(dict(file=f.name,bytes=f.stat().st_size,sha256=sha(f),type="xlsx",sheets="",read_status="ERROR:"+str(e)))
        else: rows.append(dict(file=f.name,bytes=f.stat().st_size,sha256=sha(f),type=f.suffix.lower(),sheets="",read_status="INVENTORIED_NON_TABLE"))
    pd.DataFrame(rows).to_csv(M04/"results/M04_PUBLIC_SUPPLEMENT_INVENTORY.tsv",sep="\t",index=False)
    pd.DataFrame(fields).to_csv(M04/"results/M04_PUBLIC_FIELD_AUDIT.tsv",sep="\t",index=False)
    missing=["axonal density / myelinated fiber density","continuous axonal loss","nerve fiber counts","morphometric measures","donor-matched cell composition","pain or symptom score","diabetes duration","neuropathy duration","HbA1c","paired tibial-sural quantitative pathology"]
    pd.DataFrame([dict(candidate_field=x,status="NOT_AVAILABLE_WITH_DIRECT_PUBLIC_DONOR_KEY",reason="complete public supplement audit found no documented one-to-one donor key supporting this analysis") for x in missing]).to_csv(M04/"results/M04_NOT_AVAILABLE_PUBLICLY.tsv",sep="\t",index=False)
    for name in ["M04_CONTINUOUS_PATHOLOGY_ASSOCIATIONS.tsv","M04_COMPOSITION_SENSITIVITY.tsv","M04_PAIRED_NERVE_CLINICAL.tsv"]:
        pd.DataFrame([dict(status="NOT_EVALUABLE",reason="no direct public donor/sample key for a documented quantitative endpoint")]).to_csv(M04/"results"/name,sep="\t",index=False)
    summary={"supplement_files":len(rows),"xlsx_files":sum(r['type']=='xlsx' for r in rows),"keyword_hits":len(fields),"continuous_pathology":"NOT_EVALUABLE","composition":"NOT_EVALUABLE","paired_pathology":"NOT_EVALUABLE"}
    (M04/"logs/RUN_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary

if __name__=="__main__":
    out={"M01":run_m01(),"M02":run_m02(),"M04":run_m04()}
    (RUN/"logs/01_execute_m01_m02_m04.json").write_text(json.dumps(out,indent=2,default=str),encoding="utf-8")
    print(json.dumps(out,indent=2,default=str))
