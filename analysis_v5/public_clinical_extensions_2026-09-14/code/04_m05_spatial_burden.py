import os
from pathlib import Path
from collections import defaultdict, deque
import csv,gzip,hashlib,itertools,json,re,time
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats
import matplotlib.pyplot as plt

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
RUN=P/"analysis_v5/public_clinical_extensions_2026-09-14";OUT=RUN/"05_GSE295206/results";FIG=RUN/"05_GSE295206/figures";LOG=RUN/"05_GSE295206/logs";MAN=RUN/"05_GSE295206/inputs_manifest"
for x in (OUT,FIG,LOG,MAN):x.mkdir(parents=True,exist_ok=True)
ROOT=P/"analysis_v4/targeted_repair_2026-09-12";INP=ROOT/"07_resource_recovery/M03_selected_inputs";OLDRES=ROOT/"07_resource_recovery/M03_results";CLOSE=P/"analysis_v4/closeout_20260912T062015Z/04_M03_spatial";MEM=pd.read_csv(P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv",sep="\t",dtype={"gene_id":str});MARK=pd.read_csv(P/"analysis_v3/batch03/results/REFERENCE_MARKER_PANELS.tsv",sep="\t",dtype={"canonical_ID":str});INFO=P/"data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz"
LAB={"original_early_allcell":"P1","original_late_allcell":"P2","original_late_neuron":"P3","original_severity":"P4","original_xenium":"P5","late_shared_concordant_neuronal_core":"P6","late_neuron_residual":"P7","late_allcell_residual":"P8","severity_neuron_shared_concordant_core":"P9","severity_neuron_residual":"P10"}
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));z=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q[o]=np.minimum(z,1);return q
def dec(x):return x.decode() if isinstance(x,(bytes,np.bytes_)) else str(x)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def exact_spear(x,y):
 obs=stats.spearmanr(x,y).statistic;vals=[stats.spearmanr(x,np.asarray(y)[list(p)]).statistic for p in itertools.permutations(range(len(y)))];return float(obs),float(np.mean(np.abs(vals)>=abs(obs)-1e-15)),len(vals)
def loo_range(x,y):
 v=[]
 for i in range(len(x)):
  k=np.arange(len(x))!=i;v.append(stats.spearmanr(np.asarray(x)[k],np.asarray(y)[k]).statistic)
 return float(np.nanmin(v)),float(np.nanmax(v))

# Deterministic symbol resolver used by the existing spatial branch.
gi=pd.read_csv(INFO,sep="\t",dtype=str,keep_default_na=False)[["GeneID","Symbol","Synonyms"]];maps=[defaultdict(set) for _ in range(3)]
for r in gi.itertuples(index=False):
 maps[0][r.Symbol].add(r.GeneID);maps[1][r.Symbol.upper()].add(r.GeneID)
 if r.Synonyms and r.Synonyms!='-':
  for s in r.Synonyms.split('|'):maps[2][s.upper()].add(r.GeneID)
exact={k:next(iter(v)) for k,v in maps[0].items() if len(v)==1};fold={k:next(iter(v)) for k,v in maps[1].items() if len(v)==1};syn={k:next(iter(v)) for k,v in maps[2].items() if len(v)==1 and k not in fold}
def resolve(s):return exact.get(s) or fold.get(s.upper()) or syn.get(s.upper())

section=pd.read_csv(OLDRES/"M03_SECTION_ROI_EFFECTS.tsv",sep="\t");donfx=pd.read_csv(OLDRES/"M03_DONOR_EFFECTS.tsv",sep="\t");coords=pd.read_csv(CLOSE/"M03_REAL_COORDINATES_WITH_AUTHOR_ROI.tsv.gz",sep="\t")
qcs=pd.read_csv(OLDRES/"M03_SECTION_QC.tsv",sep="\t");gsmmap=qcs.set_index('geo_accession')[['sample_id','donor_id','DRG_ID']]
sec_b=coords.groupby(["geo_accession","roi"]).size().unstack(fill_value=0).reset_index();sec_b=sec_b.merge(qcs[["geo_accession","sample_id","donor_id","DRG_ID"]],on="geo_accession");sec_b["all_author_labelled_eligible_spots"]=sec_b.get("Nageotte",0)+sec_b.get("Neuronal",0);sec_b["mixed_excluded_spots"]=sec_b.get("mixed_excluded",0);sec_b["nageotte_spot_fraction"]=sec_b.get("Nageotte",0)/sec_b.all_author_labelled_eligible_spots
sec_b.to_csv(OUT/"M05_NAGEOTTE_BURDEN_BY_SECTION.tsv",sep="\t",index=False)
drg_b=sec_b.groupby(["donor_id","DRG_ID"]).nageotte_spot_fraction.mean().reset_index();don_b=drg_b.groupby("donor_id").nageotte_spot_fraction.mean().reset_index();don_b.to_csv(OUT/"M05_NAGEOTTE_BURDEN_BY_DONOR.tsv",sep="\t",index=False)

whole=[];proxy=[];ring=[];start=time.time()
context={"Schwann":["Schwann"],"immune_macrophage":["Macrophage","Lymphoid"],"endothelial_vascular":["Vascular_endothelial","Mural"],"stromal_perineurial_fibroblast":["Endoneurial","Epineurial","Perineurial"]}
panel={c:set(MARK[MARK.celltype.isin(types)].canonical_ID.astype(str)) for c,types in context.items()}
for hp in sorted(INP.glob('*filtered_feature_bc_matrix.h5')):
 gsm=hp.name.split('_')[0];sm=gsmmap.loc[gsm];cf=next(INP.glob(gsm+'*tissue_positions*.csv.gz'))
 with h5py.File(hp,'r') as h:
  m=h['matrix'];bars=[dec(x) for x in m['barcodes'][:]];names=[dec(x) for x in m['features/name'][:]];mat=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
 gids=[resolve(x) for x in names];ok=[i for i,g in enumerate(gids) if g];ug=sorted(set(gids[i] for i in ok),key=int);ix={g:i for i,g in enumerate(ug)};agg=sparse.csr_matrix((np.ones(len(ok)),([ix[gids[i]] for i in ok],ok)),shape=(len(ug),len(names)))@mat;nz=np.asarray(agg.sum(1)).ravel()>0;agg=agg[nz];ug=np.asarray(ug)[nz];gix={g:i for i,g in enumerate(ug)}
 dense=agg.toarray().astype(np.float32);ranks=stats.rankdata(dense,axis=0,method='average')/dense.shape[0]
 for mod,g in MEM.groupby('module_id'):
  up=[gix[x] for x in g[g.direction.eq('up')].gene_id if x in gix];dn=[gix[x] for x in g[g.direction.eq('down')].gene_id if x in gix]
  if len(up)>=10 and len(dn)>=10:whole.append(dict(geo_accession=gsm,sample_id=sm.sample_id,donor_id=sm.donor_id,DRG_ID=sm.DRG_ID,module_id=mod,program=LAB[mod],whole_section_score=float(np.mean(ranks[up].mean(0)-ranks[dn].mean(0))),measured_up=len(up),measured_down=len(dn)))
 totals=np.asarray(agg.sum(0)).ravel();barix={b:i for i,b in enumerate(bars)};roi=coords[coords.geo_accession.eq(gsm)].set_index('barcode').roi.to_dict()
 for c,gs in panel.items():
  mids=[gix[x] for x in gs if x in gix];v=np.log1p(1e4*np.asarray(agg[mids].sum(0)).ravel()/np.maximum(totals,1)) if mids else np.full(len(bars),np.nan)
  for rr in ('Nageotte','Neuronal'):
   ids=[barix[b] for b,r in roi.items() if r==rr and b in barix];proxy.append(dict(geo_accession=gsm,sample_id=sm.sample_id,donor_id=sm.donor_id,DRG_ID=sm.DRG_ID,context=c,roi=rr,n_spots=len(ids),measured_markers=len(mids),proxy_mean=float(np.nanmean(v[ids])) if ids else np.nan))
 proxy.append(dict(geo_accession=gsm,sample_id=sm.sample_id,donor_id=sm.donor_id,DRG_ID=sm.DRG_ID,context='neuronal',roi='NOT_EVALUABLE',n_spots=0,measured_markers=0,proxy_mean=np.nan))
 # graph steps: nearest observed array-lattice connections inferred from minimal coordinate distance.
 with gzip.open(cf,'rt',newline='') as f: cr=list(csv.reader(f));
 if cr and cr[0][0].lower()=='barcode':cr=cr[1:]
 pos={r[0]:(int(r[2]),int(r[3])) for r in cr if len(r)>=4 and r[0] in barix};nodes=list(pos);arr=np.array([pos[b] for b in nodes]);nage=[i for i,b in enumerate(nodes) if roi.get(b)=='Nageotte']
 if nage and len(nodes)>1:
  # connect each node to its six closest coordinate neighbours; unit is explicitly graph step.
  from scipy.spatial import cKDTree
  tree=cKDTree(arr);_,nei=tree.query(arr,k=min(7,len(nodes)));adj=[set() for _ in nodes]
  for i,js in enumerate(nei[:,1:]):
   for j in np.atleast_1d(js):adj[i].add(int(j));adj[int(j)].add(i)
  dist=np.full(len(nodes),999,int);dq=deque(nage)
  for i in nage:dist[i]=0
  while dq:
   i=dq.popleft()
   for j in adj[i]:
    if dist[j]>dist[i]+1:dist[j]=dist[i]+1;dq.append(j)
  for mod,g in MEM.groupby('module_id'):
   up=[gix[x] for x in g[g.direction.eq('up')].gene_id if x in gix];dn=[gix[x] for x in g[g.direction.eq('down')].gene_id if x in gix]
   if len(up)<10 or len(dn)<10:continue
   score=ranks[up].mean(0)-ranks[dn].mean(0)
   for rr,mask in [('0',dist==0),('1',dist==1),('2',dist==2),('3+',dist>=3)]:
    ids=[barix[nodes[i]] for i in np.where(mask)[0]];ring.append(dict(geo_accession=gsm,donor_id=sm.donor_id,program=LAB[mod],ring=rr,unit='spot_graph_steps',n_spots=len(ids),mean_score=float(np.mean(score[ids])) if ids else np.nan))

whole=pd.DataFrame(whole);whole.to_csv(OUT/"M05_WHOLE_SECTION_SCORES.tsv",sep='\t',index=False);wd=whole.groupby(['donor_id','program']).whole_section_score.mean().reset_index()
# Adjacent-neuronal and ROI contrast from corrected existing section-level scores, section->DRG->donor.
adj=section.groupby(['donor_id','DRG_ID','program']).agg(adjacent=('S_mean_Neuronal','mean'),contrast=('S_delta','mean'),up_contrast=('U_delta','mean'),down_contrast=('D_delta','mean')).reset_index().groupby(['donor_id','program']).mean(numeric_only=True).reset_index();allv=wd.merge(adj,on=['donor_id','program'],how='outer').merge(don_b,on='donor_id')
assoc=[]
for outcome in ['whole_section_score','adjacent','contrast']:
 for prog,g in allv.groupby('program'):
  q=g[[outcome,'nageotte_spot_fraction']].dropna();
  if len(q)==6:
   rho,pv,nperm=exact_spear(q[outcome],q.nageotte_spot_fraction);lo,hi=loo_range(q[outcome],q.nageotte_spot_fraction);assoc.append(dict(program=prog,outcome=outcome,n_donors=6,rho=rho,P_exact=pv,permutations=nperm,LOO_rho_min=lo,LOO_rho_max=hi,status='COMPLETED_EXPLORATORY'))
  else:assoc.append(dict(program=prog,outcome=outcome,n_donors=len(q),status='NOT_EVALUABLE'))
a=pd.DataFrame(assoc);a['BH_q']=np.nan
for o,g in a.groupby('outcome'):a.loc[g.index,'BH_q']=bh(g.P_exact)
a.to_csv(OUT/"M05_NAGEOTTE_BURDEN_PROGRAM_ASSOCIATIONS.tsv",sep='\t',index=False)
arms=[]
for arm in ['up_contrast','down_contrast']:
 for prog,g in allv.groupby('program'):
  q=g[[arm,'nageotte_spot_fraction']].dropna();
  if len(q)==6:
   rho,pv,nperm=exact_spear(q[arm],q.nageotte_spot_fraction);lo,hi=loo_range(q[arm],q.nageotte_spot_fraction);arms.append(dict(program=prog,arm=arm,n_donors=6,rho=rho,P_exact=pv,permutations=nperm,LOO_rho_min=lo,LOO_rho_max=hi,status='COMPLETED_EXPLORATORY'))
  else:arms.append(dict(program=prog,arm=arm,n_donors=len(q),status='NOT_EVALUABLE'))
arms=pd.DataFrame(arms);arms['BH_q']=bh(arms.P_exact);arms.to_csv(OUT/"M05_NAGEOTTE_BURDEN_ARM_ASSOCIATIONS.tsv",sep='\t',index=False)
pr=pd.DataFrame(proxy);pv=pr[pr.roi.isin(['Nageotte','Neuronal'])].pivot_table(index=['donor_id','DRG_ID','context'],columns='roi',values='proxy_mean').reset_index();pv['proxy_delta']=pv.Nageotte-pv.Neuronal;pd.DataFrame(proxy).to_csv(OUT/"M05_SPATIAL_MARKER_PROXY_BY_SECTION.tsv",sep='\t',index=False);pvd=pv.groupby(['donor_id','context']).proxy_delta.mean().reset_index();corr=[]
for c,g in pvd.groupby('context'):
 for prog,h in adj.groupby('program'):
  q=g.merge(h[['donor_id','contrast']],on='donor_id').dropna();rho,pv0,nperm=exact_spear(q.proxy_delta,q.contrast);corr.append(dict(context=c,program=prog,n_donors=len(q),rho_proxy_delta_vs_program_contrast=rho,P_exact=pv0,status='DESCRIPTIVE_SENSITIVITY'))
pd.DataFrame(corr).to_csv(OUT/"M05_SPATIAL_MARKER_PROXY_SUMMARY.tsv",sep='\t',index=False)
ring=pd.DataFrame(ring);ring.to_csv(OUT/"M05_TOPOLOGICAL_RING_ANALYSIS.tsv",sep='\t',index=False)
fig,ax=plt.subplots(figsize=(7,4));ax.bar(don_b.donor_id,don_b.nageotte_spot_fraction,color='#486b8a');ax.tick_params(axis='x',rotation=45);ax.set_ylabel('Nageotte / author-labelled spots');fig.tight_layout();fig.savefig(FIG/"M05_F2_burden.png",dpi=300);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4));aa=a[a.outcome.eq('contrast')].sort_values('rho');ax.scatter(aa.rho,np.arange(len(aa)),color='#486b8a');ax.axvline(0,color='black',lw=.7);ax.set_yticks(np.arange(len(aa)),aa.program);ax.set_xlabel('Spearman rho: burden vs ROI contrast');fig.tight_layout();fig.savefig(FIG/"M05_F3_program_burden.png",dpi=300);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5));x=coords[coords.geo_accession.isin(coords.geo_accession.unique()[:4])];colors={'Nageotte':'#a64b4b','Neuronal':'#486b8a','mixed_excluded':'#999999'}
for i,(gsm,g) in enumerate(x.groupby('geo_accession')):ax.scatter(g.pixel_col,g.pixel_row,s=5,c=g.roi.map(colors),alpha=.6,label=gsm if i<4 else None)
ax.invert_yaxis();ax.set_xlabel('Public pixel column');ax.set_ylabel('Public pixel row');fig.tight_layout();fig.savefig(FIG/"M05_F1_real_coordinates.png",dpi=300);plt.close(fig)
inputs=list(INP.glob('*filtered_feature_bc_matrix.h5'))+[CLOSE/"M03_REAL_COORDINATES_WITH_AUTHOR_ROI.tsv.gz",OLDRES/"M03_SECTION_ROI_EFFECTS.tsv",P/"analysis_v3/batch03/results/REFERENCE_MARKER_PANELS.tsv",P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv"]
pd.DataFrame([dict(path=str(x),bytes=x.stat().st_size,sha256=sha(x)) for x in inputs]).to_csv(MAN/"INPUTS.tsv",sep='\t',index=False)
summary={'status':'COMPLETED_EXPLORATORY','sections':16,'DRGs':7,'donors':6,'evaluable_programs':int(a.program.nunique()),'program_burden_q_lt_0.05':int((a.BH_q<.05).sum()),'arm_q_lt_0.05':int((arms.BH_q<.05).sum()),'topology_unit':'spot_graph_steps','marker_proxies':'expression proxies, not cell fractions','seconds':time.time()-start};(LOG/"RUN_SUMMARY.json").write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
