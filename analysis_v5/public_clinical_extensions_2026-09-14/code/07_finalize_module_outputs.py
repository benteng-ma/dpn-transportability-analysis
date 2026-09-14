import os
from pathlib import Path
import csv, gzip, hashlib, json, platform, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from openpyxl import load_workbook

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"
PY = Path(sys.executable)
def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);f=np.isfinite(p);pf=p[f];o=np.argsort(pf);z=np.minimum.accumulate((pf[o]*len(pf)/np.arange(1,len(pf)+1))[::-1])[::-1];qf=np.empty(len(pf));qf[o]=np.minimum(z,1);q[f]=qf;return q
def sha(f):
    h=hashlib.sha256()
    with open(f,"rb") as z:
        for b in iter(lambda:z.read(8<<20),b""):h.update(b)
    return h.hexdigest()

# M01 implementation correction receipt.
m01=RUN/"01_GSE302658"
receipt={"superseded_scope":"initial HC3 inference used asymptotic-normal reference","reason":"frozen CodeKit requires residual-df Student-t reference","unchanged":"designs, coefficients, HC3 standard errors and analysis families","replacement":"canonical results now use HC3 covariance with residual-df t inference","verification":json.loads((RUN/"tests/M01_VERIFICATION.json").read_text())}
(m01/"logs/M01_IMPLEMENTATION_CORRECTION.json").write_text(json.dumps(receipt,indent=2),encoding="utf-8")

# M02-F2: compare sex-only and public Sentrix-adjusted effects for all eligible CpGs.
m02=RUN/"02_GSE286347"; old=P/"analysis_v4/targeted_repair_2026-09-12/07_resource_recovery/M02_results"
# The analysis QC table contains the 231 DPN cases.  Build the required 315-row
# public inventory directly from the official GEO series matrix and merge the
# richer DPN QC fields without using controls in the painful/painless models.
series=P/"analysis_v4/all_extensions_2026-09-11/01_inputs/public_downloads/GSE286347/GSE286347_series_matrix.txt.gz"
raw={}; chars=[]; supp=[]
with gzip.open(series,"rt",encoding="utf-8",errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"): break
        if not line.startswith("!Sample_"): continue
        z=next(csv.reader([line],delimiter="\t")); key=z[0][8:]; vals=z[1:]
        if key=="characteristics_ch1": chars.append(vals)
        elif key=="supplementary_file": supp.append(vals)
        elif key in ["title","geo_accession","platform_id","source_name_ch1"]: raw[key]=vals
n=len(raw["geo_accession"]); full=pd.DataFrame({k:v for k,v in raw.items()})
for vals in chars:
    parsed=[x.split(":",1) for x in vals]
    if all(len(x)==2 for x in parsed): full[parsed[0][0].strip().lower()]=[x[1].strip() for x in parsed]
if supp: full["green_idat_url"]=supp[0]
full["public_series_n"]=n; full["analysis_role"]=np.where(full.get("phenotype","").str.lower().eq("control"),"CONTROL_NOT_PRIMARY","DPN_PAIN_ANALYSIS")
dpn=pd.read_csv(m02/"results/M02_PUBLIC_METADATA_AUDIT.tsv",sep="\t");full=full.merge(dpn,left_on="geo_accession",right_on="geo_accession",how="left",suffixes=("_geo","_analysis"),validate="one_to_one")
full.to_csv(m02/"results/M02_PUBLIC_METADATA_AUDIT.tsv",sep="\t",index=False)
fig,axs=plt.subplots(1,2,figsize=(10,4.5),sharex=False,sharey=False)
for ax,cohort in zip(axs,["PROPGER","PROPENG"]):
    a=pd.read_csv(old/f"M02_{cohort}_EWAS_ALL_ELIGIBLE.tsv.gz",sep="\t",usecols=["probe_id","M_value_difference"]).rename(columns={"M_value_difference":"sex_only"})
    b=pd.read_csv(m02/f"results/M02_{cohort}_EWAS_SEX_BATCH_ALL.tsv.gz",sep="\t",usecols=["probe_id","M_value_difference"]).rename(columns={"M_value_difference":"sex_batch"})
    x=a.merge(b,on="probe_id",validate="one_to_one");ax.hexbin(x.sex_only,x.sex_batch,gridsize=70,mincnt=1,cmap="Blues");lo=min(x.sex_only.min(),x.sex_batch.min());hi=max(x.sex_only.max(),x.sex_batch.max());ax.plot([lo,hi],[lo,hi],color="black",lw=.7);ax.set_title(cohort);ax.set_xlabel("Sex-only M-value effect");ax.set_ylabel("Sex + Sentrix effect")
fig.tight_layout();fig.savefig(m02/"figures/M02_F2_MODEL_EFFECT_COMPARISON.png",dpi=300);plt.close(fig)

# M04 explicit public-key audit. Categorical morphology already used historically
# is retained; it is not transformed into a continuous endpoint.
m04=RUN/"04_JCI184075"; supp=P/"data/raw/human_sural_nerve_JCI184075/supplementary"
s171=supp/"jci-135-184075-s171.xlsx"; d=pd.read_excel(s171,sheet_name="DPN_patient_sample_info")
mor="Sural nerve morphology assessment (Toluidine blue)"; cols={
 "age":"Age","sex":"Sex","diabetes_status":"Diabetes/control","morphology_category":mor,
 "sural_tibial_sample_labels":"Bulk RNA-seq Nextseq500 sural versus tibial","severity_sample_label":"Bulk RNA-seq Novaseq moderate versus severe axonal loss",
 "visium_sample_label":"Visium","pain_score":"ABSENT","HbA1c":"ABSENT","diabetes_duration":"ABSENT","neuropathy_duration":"ABSENT",
 "continuous_axonal_density":"ABSENT","quantitative_cell_composition":"ABSENT","paired_quantitative_pathology":"ABSENT"}
rows=[]
for field,c in cols.items():
    if c=="ABSENT": rows.append({"field":field,"source_file":s171.name,"public":False,"n_nonmissing":0,"direct_key":False,"scale":"absent","new_analysis":"NOT_EVALUABLE"})
    else:
        v=d[c].dropna();v=v[~v.astype(str).str.lower().isin(["not performed","na","n/a"])]
        scale="continuous" if field=="age" else ("ordinal_category" if field=="morphology_category" else "categorical_or_sample_label")
        rows.append({"field":field,"source_file":s171.name,"public":True,"n_nonmissing":len(v),"direct_key":True,"scale":scale,"new_analysis":"REUSE_HISTORICAL_ONLY" if field=="morphology_category" else "AUDIT_ONLY"})
audit=pd.DataFrame(rows);audit.to_csv(m04/"results/M04_DIRECT_KEY_AUDIT.tsv",sep="\t",index=False)
mc=d[mor].dropna().astype(str);mc=mc[~mc.str.lower().isin(["not performed","na","n/a"])]
mc.value_counts().rename_axis("deposited_morphology_category").reset_index(name="n_records").to_csv(m04/"results/M04_MORPHOLOGY_CATEGORY_COUNTS.tsv",sep="\t",index=False)
status=pd.DataFrame([
 {"analysis":"continuous pathology association","status":"NOT_EVALUABLE","reason":"No public donor-linked continuous pathology measure; only categorical axonal-loss grades already analyzed."},
 {"analysis":"composition sensitivity","status":"NOT_EVALUABLE","reason":"No public donor-linked quantitative cell fractions/counts for the fixed sural program-score cohort."},
 {"analysis":"paired tibial-sural quantitative pathology","status":"NOT_EVALUABLE","reason":"Direct S/T RNA sample labels exist, but no quantitative pathology measure in both nerves."},
 {"analysis":"categorical moderate-versus-severe axonal loss","status":"REUSED_HISTORICAL","reason":"Existing corrected analysis retained; not rerun to seek a different P value."},
 {"analysis":"S174 pain history","status":"AUDIT_ONLY_NOT_LINKED","reason":"Six DRG/sciatic proteomic donors are a different assay cohort and are not treated as a sural-score clinical crosswalk."}
]);status.to_csv(m04/"results/M04_ANALYSIS_STATUS.tsv",sep="\t",index=False)
fig,ax=plt.subplots(figsize=(8,4.8));tmp=audit.assign(label=np.where(audit.new_analysis.eq("NOT_EVALUABLE"),"Absent / NE",np.where(audit.new_analysis.eq("REUSE_HISTORICAL_ONLY"),"Existing categorical result","Public audit field"))).label.value_counts();ax.barh(tmp.index,tmp.values,color=["#8a8a8a","#486b8a","#9d6b53"][:len(tmp)]);ax.set_xlabel("Audited field count");ax.set_ylabel("");fig.tight_layout();fig.savefig(m04/"figures/M04_F1_PUBLIC_FIELD_AVAILABILITY.png",dpi=300);plt.close(fig)

# M05 missing required displays from existing computed outputs.
m05=RUN/"05_GSE295206"; res=m05/"results"; figd=m05/"figures"
arm=pd.read_csv(res/"M05_NAGEOTTE_BURDEN_ARM_ASSOCIATIONS.tsv",sep="\t");hm=arm.pivot(index="program",columns="arm",values="rho")
fig,ax=plt.subplots(figsize=(6,6));sns.heatmap(hm,cmap="vlag",center=0,vmin=-1,vmax=1,annot=True,fmt=".2f",ax=ax,cbar_kws={"label":"Spearman rho"});ax.set_xlabel("");ax.set_ylabel("");fig.tight_layout();fig.savefig(figd/"M05_F4_ARM_SPECIFIC_BURDEN.png",dpi=300);plt.close(fig)
mk=pd.read_csv(res/"M05_SPATIAL_MARKER_PROXY_SUMMARY.tsv",sep="\t");hm=mk.pivot(index="program",columns="context",values="rho_proxy_delta_vs_program_contrast")
fig,ax=plt.subplots(figsize=(8,6));sns.heatmap(hm,cmap="vlag",center=0,vmin=-1,vmax=1,ax=ax,cbar_kws={"label":"Spearman rho"});ax.set_xlabel("Fixed external marker-expression proxy");ax.set_ylabel("");fig.tight_layout();fig.savefig(figd/"M05_F5_MARKER_PROXY_ROI_DIFFERENCES.png",dpi=300);plt.close(fig)
ring=pd.read_csv(res/"M05_TOPOLOGICAL_RING_ANALYSIS.tsv",sep="\t");dr=ring.groupby(["donor_id","program","ring"],as_index=False).mean(numeric_only=True);avg=dr.groupby(["program","ring"],as_index=False).agg(mean_score=("mean_score","mean"),n_donors=("donor_id","nunique"));order=["0","1","2","3+"];avg["ring"]=pd.Categorical(avg.ring.astype(str),order,ordered=True)
fig,axs=plt.subplots(3,3,figsize=(10,8),sharex=True); programs=sorted(avg.program.unique(),key=lambda x:int(x[1:]));
for ax,p in zip(axs.flat,programs):
    z=avg[avg.program.eq(p)].sort_values("ring");ax.plot(z.ring.astype(str),z.mean_score,marker="o",color="#486b8a");ax.axhline(0,color="black",lw=.5);ax.set_title(f"{p}, donors={int(z.n_donors.min()) if len(z) else 0}")
for ax in axs[-1,:]:ax.set_xlabel("Spot graph steps")
for ax in axs[:,0]:ax.set_ylabel("Mean fixed score")
fig.tight_layout();fig.savefig(figd/"M05_F6_TOPOLOGICAL_RING_TRENDS.png",dpi=300);plt.close(fig)

# Semantic tests target scientific invariants rather than file counts alone.
tests=[]
def add(name,passed,observed,expected):tests.append({"test":name,"passed":bool(passed),"observed":str(observed),"expected":str(expected)})
v=json.loads((RUN/"tests/M01_VERIFICATION.json").read_text());add("M01 independent HC3 t-reference replay",v["all_pass"],v["max_abs_P_diff"],"<=1e-10")
wide=pd.read_csv(m01/"results/GSE302658_CLINICAL_ANALYSIS_WIDE.tsv",sep="\t");add("M01 participant count",wide.participant_id.nunique()==104,wide.participant_id.nunique(),104)
md=pd.read_csv(m02/"results/M02_PUBLIC_METADATA_AUDIT.tsv",sep="\t");add("M02 all public samples",len(md)==315,len(md),315);add("M02 DPN samples",md.cohort.isin(["PROPGER","PROPENG"]).sum()==231,md.cohort.isin(["PROPGER","PROPENG"]).sum(),231)
for cohort in ["PROPGER","PROPENG"]:
    x=pd.read_csv(m02/f"results/M02_{cohort}_EWAS_SEX_BATCH_ALL.tsv.gz",sep="\t");q=bh(x.P);add(f"M02 {cohort} full-family BH",np.nanmax(np.abs(q-x.q))<1e-12,np.nanmax(np.abs(q-x.q)),"<1e-12")
de=pd.read_csv(RUN/"03_GSE14806x/results/M03_RRBS_DIFFERENTIAL_RESULTS.tsv.gz",sep="\t");q=bh(de.P);add("M03 finite-family BH",np.nanmax(np.abs(q-de.BH_q))<1e-12,np.nanmax(np.abs(q-de.BH_q)),"<1e-12")
sa=pd.read_csv(RUN/"03_GSE14806x/results/M03_RRBS_SAMPLE_AUDIT.tsv",sep="\t");add("M03 no fabricated patient pairing",sa.pairing_status.eq("BLOCKED").all(),sa.pairing_status.value_counts().to_dict(),"all BLOCKED")
bur=pd.read_csv(res/"M05_NAGEOTTE_BURDEN_BY_DONOR.tsv",sep="\t");add("M05 donor inference count",bur.donor_id.nunique()==6,bur.donor_id.nunique(),6);add("M05 P9 remains NE/absent",not arm.program.eq("P9").any(),arm.program.unique().tolist(),"P9 absent")
pd.DataFrame(tests).to_csv(RUN/"tests/MODULE_SEMANTIC_TESTS.tsv",sep="\t",index=False)

# Environment and actual source hashes used by this finalizer.
env=RUN/"00_admin/environment";env.mkdir(parents=True,exist_ok=True)
(env/"PYTHON_ENVIRONMENT.txt").write_text(f"executable={sys.executable}\nversion={sys.version}\nplatform={platform.platform()}\n",encoding="utf-8")
sources=[s171,DE_RNA if 'DE_RNA' in globals() else P/"analysis_v3/batch02/results/de/DE_ALL_GENES_C02.tsv.gz"]
pd.DataFrame([{"path":str(f),"bytes":f.stat().st_size,"sha256":sha(f)} for f in sources]).to_csv(RUN/"00_admin/FINALIZER_INPUTS.tsv",sep="\t",index=False)
print(json.dumps({"semantic_tests":len(tests),"passed":sum(x['passed'] for x in tests)},indent=2))
