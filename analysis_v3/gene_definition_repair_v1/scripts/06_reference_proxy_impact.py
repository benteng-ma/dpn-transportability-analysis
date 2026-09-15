from common import *
import pandas as pd,numpy as np,shutil
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
def imported(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.isfinite(p);out[ok]=multipletests(p[ok],method='fdr_bh')[1] if ok.any() else [];return out
def main():
 b3=P/'analysis_v3/batch03';refdir=B/'reference_impact';(refdir/'results').mkdir(parents=True,exist_ok=True)
 ra=imported(b3/'scripts/reference_analysis.py','reference_analysis');ra.B=b3;sys.modules['reference_analysis']=ra
 members=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str});m12=load('12_validate_hdrg_components_in_human_sural_nerve.py');look=m12.build_symbol_lookup(m12.load_gene_info(m12.NCBI/'Homo_sapiens.gene_info.gz'))
 cov=[];maps=[]
 for acc in ['GSE168243','GSE285983']:
  genes,k,cpm,oldranks,det=ra.reference(acc)
  ann=m12.add_resolution(pd.DataFrame({'gene':genes}),'gene',look);ann['row']=np.arange(len(genes));maps.append(ann.assign(reference=acc))
  keep=ann.human_gene_id.notna();expr=pd.DataFrame(np.log2(cpm[keep]+1),index=ann.loc[keep,'human_gene_id']).groupby(level=0,sort=True).median();rank=rankdata(expr,axis=0,method='average')/len(expr)-.5
  idx={g:i for i,g in enumerate(expr.index)};rows=[]
  for mod,f in members.groupby('module_id'):
   u=[idx[g] for g in f.loc[f.direction=='up','gene_id'] if g in idx];d=[idx[g] for g in f.loc[f.direction=='down','gene_id'] if g in idx];eligible=len(u)>=10 and len(d)>=10
   score=rank[u].mean(axis=0)-rank[d].mean(axis=0) if eligible else np.full(len(k),np.nan)
   cov.append(dict(reference=acc,object=mod,up=len(u),down=len(d),measured=len(u)+len(d),source_members=len(f),background=len(expr),eligible=eligible))
   for j,key in k.iterrows():rows.append(dict(reference=acc,**key.to_dict(),object=mod,score=score[j] if key.nuclei>=30 else np.nan,score_20_sensitivity=score[j] if key.nuclei>=20 else np.nan,coverage_eligible=eligible))
  new=pd.DataFrame(rows);new.to_csv(refdir/'results'/f'{acc}_DONOR_LOCALIZATION.tsv.gz',sep='\t',index=False)
  old=pd.read_csv(b3/'results'/f'{acc}_DONOR_LOCALIZATION.tsv.gz',sep='\t');old=old[~old.object.str.startswith('module:')]
  put(new.merge(old,on=['reference','object','donor','level','celltype'],suffixes=('_new','_old'),validate='one_to_one'),f'results/{acc}_LOCALIZATION_OLD_NEW.tsv.gz')
 put(pd.DataFrame(cov),'results/REFERENCE_COVERAGE.tsv');put(pd.concat(maps),'results/REFERENCE_TARGET_CANONICAL_MAPPING.tsv.gz')
 shutil.copyfile(b3/'results/GSE285983_DONOR_AUDIT.tsv',refdir/'results/GSE285983_DONOR_AUDIT.tsv')
 assoc=imported(b3/'scripts/13_reference_associations.py','reference_assoc_repair');assoc.B=refdir;assoc.put=lambda f,p:f.to_csv(refdir/p,sep='\t',index=False);assoc.assert_locked=lambda:None
 assert (B/'config/LOCK.json').exists();assoc.main()
 for fn in ['ALL_REFERENCE_PROGRAM_ASSOCIATIONS.tsv','ALL_REFERENCE_COVARIATE_SENSITIVITY.tsv']:
  new=pd.read_csv(refdir/'results'/fn,sep='\t');old=pd.read_csv(b3/'results'/fn,sep='\t');keys=['program','celltype','contrast','minimum_nuclei'];put(new.merge(old,on=keys,suffixes=('_new','_old'),validate='one_to_one'),'results/'+fn.replace('.tsv','_OLD_NEW.tsv'))
 # Fixed marker panels and PCs; no new marker selection or module computations.
 proxy=imported(b3/'scripts/14_marker_proxies.py','proxy_functions');markers=pd.read_csv(b3/'results/REFERENCE_MARKER_PANELS.tsv',sep='\t',dtype={'canonical_ID':str});assert not set(markers.canonical_ID)&set(members.gene_id)
 ys=pd.read_csv(B/'results/DONOR_SCORES.tsv.gz',sep='\t');ps=pd.read_csv(b3/'results/BULK_PROXY_SCORES.tsv',sep='\t');effects=[];corr=[]
 for acc in ['GSE24290','GSE148059']:
  eid=acc+('_C01' if acc=='GSE24290' else '_C02');meta=pd.read_csv(B/'inputs/matrices'/f'{eid}_metadata.tsv',sep='\t').set_index('sample_id')
  pc=pd.read_csv(b3/'results'/f'{acc}_PROXY_PC_SCORES.tsv',sep='\t').set_index('sample_id').loc[meta.index].to_numpy()
  for version in ['repaired_inherited','repaired_statfree']:
   score=ys[(ys.endpoint==eid)&(ys.version==version)].pivot(index='sample_id',columns='module_id',values='score').loc[meta.index]
   for mod in score:
    for model,pcs in [('unadjusted',None),('proxy_adjusted',pc)]:effects.append(dict(source=acc,program=mod,version=version,model=model,**proxy.fit_model(score[mod],meta,pcs)))
    for cell,f in ps[ps.source==acc].groupby('celltype'):
     vv=f.set_index('sample_id').loc[meta.index].proxy;valid=score[mod].notna().all();rho,p=spearmanr(vv,score[mod]) if valid else (np.nan,np.nan);corr.append(dict(source=acc,program=mod,version=version,celltype=cell,rho=rho,p=p,status='COMPLETED' if valid else 'NOT_EVALUABLE'))
 er=pd.DataFrame(effects);cr=pd.DataFrame(corr)
 for df,keys in [(er,['version','model']),(cr,['version'])]:
  df['q']=np.nan
  for _,ix in df.groupby(keys).groups.items():df.loc[ix,'q']=bh(df.loc[ix,'p'])
 old=pd.read_csv(b3/'results/BULK_COMPOSITION_SENSITIVITY_ALL.tsv',sep='\t');put(er.merge(old,on=['source','program','model'],suffixes=('_new','_old'),validate='many_to_one'),'results/PROXY_ADJUSTMENT_OLD_NEW.tsv')
 old=pd.read_csv(b3/'results/PROXY_PROGRAM_MODULE_CORRELATIONS.tsv',sep='\t');old=old[old.family=='proxy_program'].rename(columns={'object':'program'});put(cr.merge(old,on=['source','program','celltype'],suffixes=('_new','_old'),validate='many_to_one'),'results/PROXY_CORRELATIONS_OLD_NEW.tsv')
 dump(dict(reference_disease_interpretation='CIDP/CIAP associations, not DPN validation; VN center confounding and DPN n=2 remain unchanged',fixed_markers_overlap_repaired_programs=0,module_only_analyses_rerun=False,ADRA1B_and_candidate_expression_recomputed=False,reference_dedup='canonical median log2 CPM before rank; not counts reconstruction'),'reports/REFERENCE_PROXY_ADAPTER_AUDIT.json')
 print('Reference and fixed-proxy impact completed',flush=True)
if __name__=='__main__':run(main)
