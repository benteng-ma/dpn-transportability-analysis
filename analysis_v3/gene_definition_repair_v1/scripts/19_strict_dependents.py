from common import *
import pandas as pd,numpy as np,shutil
from scipy.stats import rankdata,spearmanr,t
from statsmodels.stats.multitest import multipletests
S=B/'strict_background_addendum'
def imp(path,name):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def put_s(x,p):
 dst=S/p;dst.parent.mkdir(parents=True,exist_ok=True);x.to_csv(dst,sep='\t',index=False)
def bh(p):
 p=np.array(p,float);q=np.full(len(p),np.nan);v=np.isfinite(p);q[v]=multipletests(p[v],method='fdr_bh')[1] if v.any() else [];return q
def main():
 assert (S/'results/ASSOCIATION_IMPACT_ALL.tsv').exists()
 m06=load('06_validate_hdrg_signatures_in_independent_human_bulk.py');m12=load('12_validate_hdrg_components_in_human_sural_nerve.py')
 ss=pd.read_csv(S/'results/DONOR_SCORES.tsv.gz',sep='\t');old=pd.read_csv(B/'results/HISTORICAL_DONOR_TEST_IMPACT.tsv',sep='\t');out=[]
 md,_=m06.load_metadata(ss[ss.endpoint=='original_hDRG_DPN_control'].sample_id.unique().tolist())
 for z in old[(old.version=='repaired_statfree')].to_dict('records'):
  f=ss[(ss.endpoint==z['endpoint'])&(ss.module_id==z['module_id'])].set_index('sample_id');score=f.score;g=f.group
  if score.notna().all():
   if z['endpoint'].startswith('original_hDRG'):
    pv,n,_=m06.sex_stratified_exact_permutation(score,md);pos=g=='DPN';neg=g=='Control'
    for age in [False,True]:
     fit=m06.ols_condition_effect(score,md,age);name='sex_age' if age else 'sex';z[name+'_beta']=fit['coefficient'];z[name+'_p']=fit['two_sided_p'];z[name+'_ci_low']=fit['coefficient']-t.ppf(.975,fit['df'])*fit['standard_error'];z[name+'_ci_high']=fit['coefficient']+t.ppf(.975,fit['df'])*fit['standard_error']
   else:
    pos=g==('DPN' if z['endpoint'].endswith('DPN_control') else 'Severe');neg=~pos;pv,n=m12.exact_label_p(np.r_[score[pos],score[neg]],int(pos.sum()))
   loo=[score[pos&(score.index!=i)].mean()-score[neg&(score.index!=i)].mean() for i in score.index]
   z.update(p=pv,mean_difference=score[pos].mean()-score[neg].mean(),allocations=n,loo_all_positive=all(v>0 for v in loo),loo_min=min(loo),loo_max=max(loo),old_recomputed_p_difference=pv-z['archived_p'])
  out.append(z)
 new=pd.DataFrame(out);new['q']=new.groupby(['endpoint','version','family'],group_keys=False).p.apply(m12.bh_adjust)
 put_s(new,'results/HISTORICAL_STRICT.tsv');put_s(pd.concat([old[old.version!='repaired_statfree'],new]),'consolidated_results/HISTORICAL_DONOR_TEST_IMPACT.tsv')
 # Only affected GSE148059 stat-free pathology/ordinal estimates are refitted.
 m=imp(P/'analysis_v3/scripts/04_run_batch01_frozen_score_associations.py','strict_b1');cfg=json.loads(m.CONFIG_PATH.read_text());eid='GSE148059_C02';f=ss[ss.endpoint==eid].copy();f['cohort']='GSE148059';f['mapping_set']='primary_unambiguous';f['module_display']=f.module_id.map(m.DISPLAY)
 meta=pd.read_csv(S/'inputs/matrices'/f'{eid}_metadata.tsv',sep='\t');meta['progression_group']=meta.group;meta['dmfd_class']=meta.group;meta['pathology_order']=meta.group.map({'Regenerator':0,'Intermediate':1,'Degenerator':2})
 r=m.run_group_analyses(f[['cohort','mapping_set','module_id','module_display','sample_id','patient_id','score']],meta[['sample_id','patient_id','progression_group','dmfd_class','pathology_order']],cfg);r['version']='repaired_statfree';put_s(r,'results/BATCH01_STRICT_GSE148059.tsv')
 orig=pd.read_csv(P/'analysis_v3/results/PATHOLOGY_ASSOCIATION_RESULTS.tsv',sep='\t').query("mapping_set=='primary_unambiguous'");r=r.merge(orig,on=['cohort','module_id','mapping_set','test_family'],suffixes=('_new','_old'),validate='one_to_one');old=pd.read_csv(B/'results/BATCH01_PATHOLOGY_OLD_NEW.tsv',sep='\t');keep=~((old.version=='repaired_statfree')&(old.cohort=='GSE148059'));put_s(pd.concat([old[keep],r]),'consolidated_results/BATCH01_PATHOLOGY_OLD_NEW.tsv')
 print('strict historical and pathology completed',flush=True)
 # Recompute the two program-only references with truly unique official GeneIDs.
 b3=P/'analysis_v3/batch03';refdir=S/'reference_impact';(refdir/'results').mkdir(parents=True,exist_ok=True)
 ra=imp(b3/'scripts/reference_analysis.py','reference_analysis');ra.B=b3;sys.modules['reference_analysis']=ra
 info=m12.load_gene_info(m12.NCBI/'Homo_sapiens.gene_info.gz');look=m12.build_symbol_lookup(info);dup=set(info.loc[info.Symbol.duplicated(keep=False),'Symbol']);look=({k:v for k,v in look[0].items() if k not in dup},look[1],{k:v for k,v in look[2].items() if k not in {d.upper() for d in dup}},look[3])
 members=pd.read_csv(S/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str});cov=[]
 for acc in ['GSE168243','GSE285983']:
  genes,keys,cpm,_,_=ra.reference(acc);ann=m12.add_resolution(pd.DataFrame({'gene':genes}),'gene',look);keep=ann.human_gene_id.notna();expr=pd.DataFrame(np.log2(cpm[keep]+1),index=ann.loc[keep,'human_gene_id']).groupby(level=0,sort=True).median();rank=rankdata(expr,axis=0,method='average')/len(expr)-.5;idx={g:i for i,g in enumerate(expr.index)};rows=[]
  put_s(ann,'inputs/'+acc+'_STRICT_MAPPING.tsv')
  for mod,f in members.groupby('module_id'):
   u=[idx[g] for g in f.loc[f.direction=='up','gene_id'] if g in idx];d=[idx[g] for g in f.loc[f.direction=='down','gene_id'] if g in idx];ok=min(len(u),len(d))>=10;sc=rank[u].mean(axis=0)-rank[d].mean(axis=0) if ok else np.full(len(keys),np.nan)
   cov.append(dict(reference=acc,object=mod,up=len(u),down=len(d),measured=len(u)+len(d),source_members=len(f),background=len(expr),eligible=ok))
   for j,key in keys.iterrows():rows.append(dict(reference=acc,**key.to_dict(),object=mod,score=sc[j] if key.nuclei>=30 else np.nan,score_20_sensitivity=sc[j] if key.nuclei>=20 else np.nan,coverage_eligible=ok))
  new=pd.DataFrame(rows);new.to_csv(refdir/'results'/f'{acc}_DONOR_LOCALIZATION.tsv.gz',sep='\t',index=False);old=pd.read_csv(b3/'results'/f'{acc}_DONOR_LOCALIZATION.tsv.gz',sep='\t');old=old[~old.object.str.startswith('module:')];put_s(new.merge(old,on=['reference','object','donor','level','celltype'],suffixes=('_new','_old'),validate='one_to_one'),f'consolidated_results/{acc}_LOCALIZATION_OLD_NEW.tsv.gz')
 put_s(pd.DataFrame(cov),'consolidated_results/REFERENCE_COVERAGE.tsv');shutil.copyfile(b3/'results/GSE285983_DONOR_AUDIT.tsv',refdir/'results/GSE285983_DONOR_AUDIT.tsv')
 a=imp(b3/'scripts/13_reference_associations.py','strict_ref_assoc');a.B=refdir;a.put=lambda f,p:f.to_csv(refdir/p,sep='\t',index=False);a.assert_locked=lambda:None;a.main()
 for fn in ['ALL_REFERENCE_PROGRAM_ASSOCIATIONS.tsv','ALL_REFERENCE_COVARIATE_SENSITIVITY.tsv']:
  new=pd.read_csv(refdir/'results'/fn,sep='\t');old=pd.read_csv(b3/'results'/fn,sep='\t');put_s(new.merge(old,on=['program','celltype','contrast','minimum_nuclei'],suffixes=('_new','_old'),validate='one_to_one'),'consolidated_results/'+fn.replace('.tsv','_OLD_NEW.tsv'))
 # Preserve fixed nuisance proxies. Only response changed; no marker or PC reselection.
 proxy=imp(b3/'scripts/14_marker_proxies.py','strict_proxy');meta=pd.read_csv(S/'inputs/matrices/GSE148059_C02_metadata.tsv',sep='\t').set_index('sample_id');pc=pd.read_csv(b3/'results/GSE148059_PROXY_PC_SCORES.tsv',sep='\t').set_index('sample_id').loc[meta.index].to_numpy();sc=ss[ss.endpoint=='GSE148059_C02'].pivot(index='sample_id',columns='module_id',values='score').loc[meta.index];ps=pd.read_csv(b3/'results/BULK_PROXY_SCORES.tsv',sep='\t');effects=[];corr=[]
 for mod in sc:
  for model,pcs in [('unadjusted',None),('proxy_adjusted',pc)]:effects.append(dict(source='GSE148059',program=mod,version='repaired_statfree',model=model,**proxy.fit_model(sc[mod],meta,pcs)))
  for cell,f in ps[ps.source=='GSE148059'].groupby('celltype'):
   rho,p=spearmanr(f.set_index('sample_id').loc[meta.index].proxy,sc[mod]) if sc[mod].notna().all() else (np.nan,np.nan);corr.append(dict(source='GSE148059',program=mod,version='repaired_statfree',celltype=cell,rho=rho,p=p,status='COMPLETED' if np.isfinite(p) else 'NOT_EVALUABLE'))
 for fn,new,keys in [('PROXY_ADJUSTMENT_OLD_NEW.tsv',pd.DataFrame(effects),['source','program','version','model']),('PROXY_CORRELATIONS_OLD_NEW.tsv',pd.DataFrame(corr),['source','program','version','celltype'])]:
  old=pd.read_csv(B/'results'/fn,sep='\t').set_index(keys);new=new.set_index(keys)
  for col in new:
   target=col+'_new' if col+'_new' in old else col
   if target in old:old.loc[new.index,target]=new[col]
  old=old.reset_index();qkeys=['version','model'] if 'model' in old else ['version'];pcol='p_new';qcol='q_new'
  for _,ix in old.groupby(qkeys).groups.items():old.loc[ix,qcol]=bh(old.loc[ix,pcol])
  put_s(old,'consolidated_results/'+fn)
 print('strict references and proxy completed',flush=True)
if __name__=='__main__':run(main)
