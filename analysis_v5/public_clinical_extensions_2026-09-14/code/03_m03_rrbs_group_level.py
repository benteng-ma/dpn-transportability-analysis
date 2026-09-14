import os
from pathlib import Path
import gzip, hashlib, json, math, time
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14"; OUT=RUN/"03_GSE14806x/results"; FIG=RUN/"03_GSE14806x/figures"; LOG=RUN/"03_GSE14806x/logs"; MAN=RUN/"03_GSE14806x/inputs_manifest"
for x in (OUT,FIG,LOG,MAN): x.mkdir(parents=True,exist_ok=True)
RAW=P/"analysis_v3/data/raw/GSE148060/extracted"; META=P/"analysis_v3/data/processed/GSE148060_sample_metadata.tsv"; RNA_META=P/"analysis_v3/batch02/inputs/GSE148059_metadata.tsv"; RNA_RES=P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/results/BATCH01_STRICT_GSE148059.tsv"; MEMBERS=P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"; REF=RUN/"00_admin/public_reference/UCSC_hg19_refGene_2026-09-14.txt.gz"
CHR={str(i):i for i in range(1,23)}|{"X":23,"Y":24,"M":25,"MT":25}
def enc(ch,pos):
    c=np.array([CHR.get(str(x).replace("chr",""),0) for x in ch],dtype=np.uint64);return (c<<np.uint64(32))|np.asarray(pos,dtype=np.uint64)
def dec(keys): return np.array(["chr"+str(int(x>>np.uint64(32))) if int(x>>np.uint64(32))<=22 else {23:"chrX",24:"chrY",25:"chrM"}.get(int(x>>np.uint64(32)),"chr0") for x in keys]),(keys & np.uint64(0xffffffff)).astype(np.int64)
def bh(p):
    p=np.asarray(p,float);q=np.full(len(p),np.nan);finite=np.isfinite(p)
    pf=p[finite];o=np.argsort(pf);z=np.minimum.accumulate((pf[o]*len(pf)/np.arange(1,len(pf)+1))[::-1])[::-1]
    qf=np.empty(len(pf));qf[o]=np.minimum(z,1);q[finite]=qf;return q
def load(f):
    d=pd.read_csv(f,sep="\t",header=None,names=["chr","start","end","percent"],na_values="NA",dtype={"chr":str,"start":np.int64,"end":np.int64,"percent":np.float32});k=enc(d.chr,d.start);o=np.argsort(k);return k[o],d.percent.to_numpy(np.float32)[o]
def update(union,counts,k,group_index):
    nu=np.union1d(union,k); new=np.zeros((3,len(nu)),dtype=np.uint16);ix=np.searchsorted(nu,union);new[:,ix]=counts;new[group_index,np.searchsorted(nu,k)]+=1;return nu,new
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

started=time.time();meta=pd.read_csv(META,sep="\t"); meta["group"]=meta.dmfd_class; order={"Regenerator":0,"Intermediate":1,"Degenerator":2}; files=[]
for r in meta.itertuples():
    f=RAW/r.source_file; files.append((r.sample_id,r.group,order[r.group],f))
meta.assign(processed_scale="CpG percent methylation",genome_build="hg19",direct_RNA_RRBS_key="NO_VERIFIED_CROSS_ASSAY_KEY",pairing_status="BLOCKED").to_csv(OUT/"M03_RRBS_SAMPLE_AUDIT.tsv",sep="\t",index=False)
rm=pd.read_csv(RNA_META,sep="\t");rm.assign(processed_scale="rlog processed expression",direct_RNA_RRBS_key="NO_VERIFIED_CROSS_ASSAY_KEY",pairing_status="BLOCKED").to_csv(OUT/"M03_RNA_SAMPLE_AUDIT.tsv",sep="\t",index=False)

union=np.array([],dtype=np.uint64);counts=np.zeros((3,0),dtype=np.uint16)
for i,(_,g,gi,f) in enumerate(files):
    k,_=load(f);union,counts=update(union,counts,k,gi);print("pass1",i+1,len(files),len(union),flush=True)
nby=np.bincount([x[2] for x in files],minlength=3); contrast=(counts[0]>=math.ceil(.7*nby[0]))&(counts[2]>=math.ceil(.7*nby[2]));trend=counts.sum(0)>=math.ceil(.7*len(files));cand=union[contrast|trend];del union,counts
mat=np.full((len(cand),len(files)),np.nan,np.float32)
for j,(_,g,gi,f) in enumerate(files):
    k,v=load(f);ix=np.searchsorted(cand,k);ok=(ix<len(cand))&(cand[np.minimum(ix,len(cand)-1)]==k);mat[ix[ok],j]=v[ok];print("pass2",j+1,len(files),flush=True)
beta=np.clip(mat/100,.001,.999);M=np.log2(beta/(1-beta));groups=np.array([x[2] for x in files]);R=groups==0;D=groups==2
def welch(a,b):
    na=np.sum(np.isfinite(a),1);nb=np.sum(np.isfinite(b),1);ma=np.nanmean(a,1);mb=np.nanmean(b,1);va=np.nanvar(a,1,ddof=1);vb=np.nanvar(b,1,ddof=1);se=np.sqrt(va/na+vb/nb);t=(ma-mb)/se;df=(va/na+vb/nb)**2/((va/na)**2/(na-1)+(vb/nb)**2/(nb-1));return ma-mb,t,df,2*stats.t.sf(abs(t),df),na,nb
eff,t,df,pv,nD,nR=welch(M[:,D],M[:,R]);be=np.nanmean(beta[:,D],1)-np.nanmean(beta[:,R],1);eligible=(nD>=math.ceil(.7*D.sum()))&(nR>=math.ceil(.7*R.sum()))
ch,pos=dec(cand);de=pd.DataFrame({"key":cand[eligible],"chr":ch[eligible],"position":pos[eligible],"n_deg":nD[eligible],"n_reg":nR[eligible],"M_difference_Deg_minus_Reg":eff[eligible],"beta_difference_Deg_minus_Reg":be[eligible],"t":t[eligible],"df":df[eligible],"P":pv[eligible]});de["BH_q"]=bh(de.P);de.to_csv(OUT/"M03_RRBS_DIFFERENTIAL_RESULTS.tsv.gz",sep="\t",index=False,compression="gzip")
# ordered trend using complete finite observations and simple slope.
X=np.tile(groups,(len(cand),1)).astype(float);mask=np.isfinite(M);n=mask.sum(1);sx=(X*mask).sum(1);sy=np.nansum(M,1);sxx=(X*X*mask).sum(1);sxy=np.nansum(M*X,1);den=sxx-sx*sx/n;slope=(sxy-sx*sy/n)/den;inter=(sy-slope*sx)/n;res=np.where(mask,M-(inter[:,None]+slope[:,None]*X),0);s2=(res*res).sum(1)/(n-2);se=np.sqrt(s2/den);tt=slope/se;pp=2*stats.t.sf(abs(tt),n-2);ok=(n>=math.ceil(.7*len(files)))&(den>0)&np.isfinite(pp)
tr=pd.DataFrame({"key":cand[ok],"chr":ch[ok],"position":pos[ok],"n":n[ok],"M_slope_per_pathology_step":slope[ok],"t":tt[ok],"P":pp[ok]});tr["BH_q"]=bh(tr.P);tr.to_csv(OUT/"M03_RRBS_ORDERED_TREND.tsv.gz",sep="\t",index=False,compression="gzip")

# hg19 refGene annotation by 100-kb bins.
cols=["bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames"]
genes=pd.read_csv(REF,sep="\t",header=None,names=cols,usecols=["chrom","strand","txStart","txEnd","name2"],compression="gzip").drop_duplicates();bins={}
for q in genes.itertuples():
    for b in range(max(0,(q.txStart-2000)//100000),(q.txEnd+2000)//100000+1): bins.setdefault((q.chrom,b),[]).append(q)
ann=[]
for r in de.itertuples():
    for q in bins.get((r.chr,r.position//100000),[]):
        tss=q.txStart if q.strand=='+' else q.txEnd
        if tss-2000<=r.position<=tss+2000: ann.append((r.key,r.chr,r.position,q.name2,"promoter"))
        elif q.txStart<=r.position<=q.txEnd: ann.append((r.key,r.chr,r.position,q.name2,"gene_body"))
ad=pd.DataFrame(ann,columns=["key","chr","position","gene_symbol","annotation"]).drop_duplicates();ad.to_csv(OUT/"M03_RRBS_GENE_ANNOTATION.tsv.gz",sep="\t",index=False,compression="gzip")
dd=de.merge(ad,on=["key","chr","position"],how="inner");gene=dd.groupby(["gene_symbol","annotation"]).agg(rrbs_effect=("M_difference_Deg_minus_Reg","median"),rrbs_stat=("t","median"),loci=("key","nunique")).reset_index()
# Reuse corrected RNA program-level results; no patient-level RNA-RRBS link asserted.
rna=pd.read_csv(RNA_RES,sep="\t");rna.to_csv(OUT/"M03_REUSED_CORRECTED_RNA_RESULTS.tsv",sep="\t",index=False)
gene.assign(rna_effect=np.nan,rna_stat=np.nan,pairing="GROUP_LEVEL_ONLY_NO_PATIENT_KEY").to_csv(OUT/"M03_CROSSOMIC_GENE_EFFECTS.tsv.gz",sep="\t",index=False,compression="gzip")
mem=pd.read_csv(MEMBERS,sep="\t");labels={"original_early_allcell":"P1","original_late_allcell":"P2","original_late_neuron":"P3","original_severity":"P4","original_xenium":"P5","late_shared_concordant_neuronal_core":"P6","late_neuron_residual":"P7","late_allcell_residual":"P8","severity_neuron_shared_concordant_core":"P9","severity_neuron_residual":"P10"};summary=[]
for m,g in mem.groupby("module_id"):
    rr=gene[(gene.annotation=="promoter")&gene.gene_symbol.isin(set(g.symbol.astype(str)))]; rx=rna[rna.module_id.eq(m)]
    summary.append(dict(program=labels[m],module_id=m,RNA_status="REUSED_CORRECTED",RNA_group_effect=rx.effect.iloc[0] if len(rx) else np.nan,RNA_q=rx.permutation_two_sided_bh_q.iloc[0] if len(rx) else np.nan,RRBS_promoter_loci=int(rr.loci.sum()) if len(rr) else 0,RRBS_promoter_median_effect=float(rr.rrbs_effect.median()) if len(rr) else np.nan,RRBS_enrichment_status="NOT_EVALUABLE_UNEQUAL_LOCUS_REPRESENTATION",RRBS_q=np.nan,directions_compatible="DESCRIPTIVE_ONLY",participants_independent="UNKNOWN_NOT_ASSUMED"))
ps=pd.DataFrame(summary);ps.to_csv(OUT/"M03_PROGRAM_CROSSOMIC_SUMMARY.tsv",sep="\t",index=False)
fig,ax=plt.subplots(figsize=(7,4));meta.group.value_counts().reindex(["Regenerator","Intermediate","Degenerator"]).plot.bar(ax=ax,color="#486b8a");ax.set_ylabel("RRBS samples");ax.set_xlabel("");fig.tight_layout();fig.savefig(FIG/"M03_F1_groups.png",dpi=300);plt.close(fig)
fig,ax=plt.subplots(figsize=(7,5));ax.hexbin(de.beta_difference_Deg_minus_Reg,de.M_difference_Deg_minus_Reg,gridsize=50,mincnt=1,cmap="Blues");ax.axhline(0,color='black',lw=.5);ax.axvline(0,color='black',lw=.5);ax.set_xlabel("Beta difference");ax.set_ylabel("M-value difference");fig.tight_layout();fig.savefig(FIG/"M03_F2_RRBS_effects.png",dpi=300);plt.close(fig)
inputs=[META,RNA_META,RNA_RES,MEMBERS,REF]+[x[3] for x in files];pd.DataFrame([dict(path=str(x),bytes=x.stat().st_size,sha256=sha(x)) for x in inputs]).to_csv(MAN/"INPUTS.tsv",sep="\t",index=False)
summaryj={"status":"COMPLETED_GROUP_LEVEL","RRBS_samples":len(files),"groups":meta.group.value_counts().to_dict(),"candidate_CpGs":len(cand),"DE_eligible":len(de),"DE_q_lt_0.05":int((de.BH_q<.05).sum()),"trend_eligible":len(tr),"trend_q_lt_0.05":int((tr.BH_q<.05).sum()),"annotated_gene_locus_rows":len(ad),"patient_pairing":"BLOCKED_NOT_USED","seconds":time.time()-started}
(LOG/"RUN_SUMMARY.json").write_text(json.dumps(summaryj,indent=2),encoding="utf-8");print(json.dumps(summaryj,indent=2))
