from common import *
import pandas as pd,numpy as np
from scipy import stats

def module_frames():
 m=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str});c=pd.read_csv(B/'inputs/CANONICAL_SOURCE_COMPARISONS.tsv',sep='\t',dtype={'gene_id':str})
 m=m.merge(c[['contrast','gene_id','median_effect']],left_on=['source_contrast','gene_id'],right_on=['contrast','gene_id'],validate='many_to_one')
 m=m.rename(columns={'symbol':'gene','median_effect':'source_log2fc'});m['audit_only']=False
 return {k:f.copy() for k,f in m.groupby('module_id')}

def main():
 mods=module_frames();m12=load('12_validate_hdrg_components_in_human_sural_nerve.py');m06=load('06_validate_hdrg_signatures_in_independent_human_bulk.py');m13=load('13_build_cross_target_component_transportability_atlas.py')
 look=m12.build_symbol_lookup(m12.load_gene_info(m12.NCBI/'Homo_sapiens.gene_info.gz'))
 scores=pd.read_csv(B/'results/DONOR_SCORES.tsv.gz',sep='\t');coverage=pd.read_csv(B/'results/COVERAGE.tsv',sep='\t')
 historical=[];metadata,_=m06.load_metadata(scores.query("endpoint=='original_hDRG_DPN_control'").sample_id.unique().tolist())
 old5=pd.read_csv(P/'results/tables/independent_human_DPN_bulk_signature_tests_2026-08-27.tsv',sep='\t');old3=pd.read_csv(P/'results/tables/cross_target_hDRG_component_transportability_tests_2026-08-27.tsv',sep='\t')
 contrast={f.source_contrast.iloc[0]:mod for mod,f in mods.items() if mod.startswith('original_')};family={contrast[r.contrast_id]:r.test_family for r in old5.itertuples()}
 oldj=pd.read_csv(P/'results/tables/JCI184075_hDRG_component_sample_tests_2026-08-27.tsv',sep='\t')
 for (eid,ver,mod),f in scores[scores.endpoint.str.startswith('original_')].groupby(['endpoint','version','module_id']):
  if eid.startswith('original_hDRG'):
   if mod not in family and mod not in m13.MODULE_IDS:continue
   fam=family.get(mod,'three_late_components');oldrow=old5[old5.contrast_id.map(contrast)==mod] if mod in family else old3[(old3.target_id=='independent_human_hDRG')&(old3.module_id==mod)]
   oldp=oldrow.iloc[0]['donor_sex_stratified_exact_one_sided_p' if mod in family else 'exact_p'];oldq=oldrow.iloc[0]['donor_exact_bh_q' if mod in family else 'exact_q']
  else:
   target='DPN_vs_control' if eid.endswith('DPN_control') else 'Severe_vs_moderate'
   allowed=m12.DPN_MODULES if target=='DPN_vs_control' else m12.SEVERITY_MODULES
   if mod not in allowed:continue
   fam=m12.family_for(target,mod);oldrow=oldj[(oldj.module_id==mod)&(oldj.expression_source=='DESeq2_normalized_counts')&(oldj.target_id==target)]
   # Actual archive's severity label may differ; match the sole non-DPN target, never patient keys.
   if oldrow.empty:oldrow=oldj[(oldj.module_id==mod)&(oldj.expression_source=='DESeq2_normalized_counts')&(oldj.target_id!='DPN_vs_control')]
   assert len(oldrow)==1
   oldp=oldrow.iloc[0].exact_one_sided_p;oldq=oldrow.iloc[0].exact_bh_q
  row=dict(endpoint=eid,version=ver,module_id=mod,family=fam,archived_p=oldp,archived_q=oldq,status='NOT_EVALUABLE',p=np.nan,mean_difference=np.nan,loo_all_positive=False)
  if f.score.notna().all():
   ss=f.set_index('sample_id').score;g=f.set_index('sample_id').group
   if eid.startswith('original_hDRG'):
    pv,n,null=m06.sex_stratified_exact_permutation(ss,metadata);pos=g=='DPN';neg=g=='Control'
    for age in [False,True]:
     fit=m06.ols_condition_effect(ss,metadata,age);name='sex_age' if age else 'sex';row[name+'_beta']=fit['coefficient'];row[name+'_p']=fit['two_sided_p'];row[name+'_ci_low']=fit['coefficient']-stats.t.ppf(.975,fit['df'])*fit['standard_error'];row[name+'_ci_high']=fit['coefficient']+stats.t.ppf(.975,fit['df'])*fit['standard_error']
   else:
    pos=g==('DPN' if eid.endswith('DPN_control') else 'Severe');neg=~pos;pv,n=m12.exact_label_p(np.r_[ss[pos],ss[neg]],int(pos.sum()))
   loo=[ss[pos & (ss.index!=i)].mean()-ss[neg & (ss.index!=i)].mean() for i in ss.index]
   row.update(status='COMPLETED',p=pv,mean_difference=ss[pos].mean()-ss[neg].mean(),loo_all_positive=all(d>0 for d in loo),loo_min=min(loo),loo_max=max(loo),allocations=n)
  historical.append(row)
 h=pd.DataFrame(historical);h['q']=h.groupby(['endpoint','version','family'],group_keys=False).p.apply(m12.bh_adjust);h['old_recomputed_p_difference']=h.p-h.archived_p
 put(h,'results/HISTORICAL_DONOR_TEST_IMPACT.tsv')
 # Three-component human/rodent/ocular tests: use unchanged historical functions and families.
 for name,fn in [('human_components',m13.analyze_independent_human),('rat_DRG_components',m13.analyze_rat_drg),('ocular_components',m13.analyze_ocular)]:
  out=fn(mods,look)
  for j,value in enumerate(out):
   if isinstance(value,pd.DataFrame):put(value,f'results/historical/{name}_{j}.tsv.gz')
   else:dump(value,f'reports/{name}_audit.json')
  print('historical adapter',name,flush=True)
 # Ten/four JCI gene-concordance tasks and deposited TPM score sensitivity.
 for eid,target in [('original_JCI_sural_DPN_control','DPN_vs_control'),('original_JCI_sural_severe_moderate','Severe_vs_moderate')]:
  meta=pd.read_csv(B/'inputs/matrices'/f'{eid}_metadata.tsv',sep='\t');groups=dict(zip(meta.sample_id,meta.group));case=target=='DPN_vs_control';path=m12.FILES['dpn_control' if case else 'axon_severity'];sheet='deseq2_results' if case else 'Deseq_mod_vs_severe_axonLoss';ids=m12.DPN_MODULES if case else m12.SEVERITY_MODULES
  u,x,mm,au=m12.load_target(target,path,sheet,groups,1 if case else -1,look);mapped,ma,cv=m12.map_modules(mods,u,target,ids,look)
  tests=m12.gene_tests(target,mapped,u,ids);put(tests,'results/historical/'+target+'_gene_tests.tsv');put(ma,'results/historical/'+target+'_gene_mapping.tsv.gz');put(cv,'results/historical/'+target+'_coverage.tsv')
  # TPM input adapter signature is inspected separately before invoking.
 print('historical neural and cross-tissue adapters completed',flush=True)
if __name__=='__main__':run(main)
