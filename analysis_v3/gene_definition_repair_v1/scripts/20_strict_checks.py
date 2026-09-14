from common import *
import pandas as pd,numpy as np,shutil
from statsmodels.stats.multitest import multipletests
S=B/'strict_background_addendum'
def main():
 checks=[]
 def ck(k,v):checks.append(dict(check=k,status='PASS' if v else 'FAIL'))
 def read(f):return pd.read_csv(S/f,sep='\t')
 be=read('consolidated_results/REAL_BENCHMARK_ALL.tsv')
 # Partial execution cannot shrink the final family. Restore 40 eligible tests per background.
 for bg,ix in be.groupby('background').groups.items():
  ix=[i for i in ix if pd.notna(be.loc[i,'p'])];ck('40_benchmark_family:'+bg,len(ix)==40);be.loc[ix,'holm_p']=multipletests(be.loc[ix,'p'],method='holm')[1]
 be.to_csv(S/'consolidated_results/REAL_BENCHMARK_ALL.tsv',sep='\t',index=False)
 mem=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str});ss=read('results/DONOR_SCORES.tsv.gz');checkscores=[]
 ck('source_members_unchanged',sha(S/'inputs/REPAIRED_MEMBERS.tsv')==sha(B/'inputs/REPAIRED_MEMBERS.tsv'))
 ck('ambiguous_ids_not_members',not set(mem.gene_id)&{'100124696','100187828','100505381'})
 for r in read('config/LOCKED_FILES.tsv').itertuples():ck('strict_lock:'+r.file,sha(S/r.file)==r.sha256)
 for r in read('inputs/MATRIX_INVENTORY.tsv').itertuples():
  x=pd.read_csv(S/r.matrix,sep='\t',index_col=0);x.index=x.index.astype(str);rank=x.rank(axis=0,pct=True,method='average')
  ck('strict_background_excludes_ambiguous:'+r.endpoint,not set(x.index)&{'100124696','100187828','100505381'})
  for mod,f in mem.groupby('module_id'):
   u=list(set(f.loc[f.direction=='up','gene_id'])&set(x.index));d=list(set(f.loc[f.direction=='down','gene_id'])&set(x.index));expected=rank.loc[u].mean()-rank.loc[d].mean() if min(len(u),len(d))>=10 else pd.Series(np.nan,index=x.columns);got=ss[(ss.endpoint==r.endpoint)&(ss.module_id==mod)].set_index('sample_id').score.loc[x.columns]
   ck('strict_score:'+r.endpoint+':'+mod,np.allclose(expected,got,equal_nan=True,atol=1e-12))
  old=pd.read_csv(B/'results/DONOR_SCORES.tsv.gz',sep='\t');old=old[(old.endpoint==r.endpoint)&(old.version=='repaired_statfree')];new=ss[ss.endpoint==r.endpoint];j=new.merge(old,on=['endpoint','module_id','sample_id'],suffixes=('_strict','_initial'),validate='one_to_one');j['score_change']=j.score_strict-j.score_initial;checkscores.append(j)
 pd.concat(checkscores).to_csv(S/'reports/SCORE_IMPLEMENTATION_IMPACT.tsv.gz',sep='\t',index=False)
 old=pd.read_csv(B/'results/ASSOCIATION_IMPACT_ALL.tsv',sep='\t').query("version=='repaired_statfree'");new=read('results/ASSOCIATION_IMPACT_ALL.tsv');j=new.merge(old,on=['endpoint','module_id'],suffixes=('_strict','_initial'),validate='one_to_one');j.to_csv(S/'reports/ASSOCIATION_IMPLEMENTATION_IMPACT.tsv',sep='\t',index=False)
 oldhist=pd.read_csv(B/'results/HISTORICAL_DONOR_TEST_IMPACT.tsv',sep='\t');newhist=read('results/HISTORICAL_STRICT.tsv');h=newhist.merge(oldhist[oldhist.version=='repaired_statfree'],on=['endpoint','module_id','version'],suffixes=('_strict','_initial'),validate='one_to_one');h.to_csv(S/'reports/HISTORICAL_IMPLEMENTATION_IMPACT.tsv',sep='\t',index=False)
 for fn in ['PROXY_ADJUSTMENT_OLD_NEW.tsv','ALL_REFERENCE_PROGRAM_ASSOCIATIONS_OLD_NEW.tsv','ALL_REFERENCE_COVARIATE_SENSITIVITY_OLD_NEW.tsv']:
  f=read('consolidated_results/'+fn);ck('complete_rows:'+fn,len(f)==(80 if fn.startswith('PROXY') else 720))
 a=read('consolidated_results/ASSOCIATION_IMPACT_ALL.tsv');ck('all_150_associations_retained',len(a)==150);ck('no_source_rule_changes',sha(B/'inputs/REPAIRED_MEMBERS.tsv')==json.loads((B/'config/LOCK.json').read_text())['members_sha256'])
 # Recheck stored strict random partition membership and count strata.
 sp=importlib.util.spec_from_file_location('mathengine',B/'scripts/02_associations_benchmarks.py');e=importlib.util.module_from_spec(sp);sp.loader.exec_module(e)
 parts=pd.read_csv(B/'inputs/REPAIRED_PARTITIONS.tsv',sep='\t',dtype={'gene_id':str});rr=read('results/MATCHED_RANDOM_SPLITS.tsv.gz')
 for (eid,fam),f in rr.groupby(['endpoint','family']):
  x=pd.read_csv(S/'inputs/matrices'/f'{eid}__statfree.tsv.gz',sep='\t',index_col=0);x.index=x.index.astype(str);u=parts[(parts.family==fam)&parts.bin.isin(['core','residual'])&parts.gene_id.isin(x.index)].sort_values('gene_id');groups,_,_,_=e.bm.matching_strata(x.loc[u.gene_id].to_numpy(),np.where(u.direction=='up',1,-1),u.bin.eq('core').to_numpy())
  for row in f.itertuples():
   ids=set(row.core_gene_ids.split('|'));mask=u.gene_id.isin(ids).to_numpy();ck(f'strict_partition:{eid}:{fam}:{row.draw}',ids<=set(u.gene_id) and all(mask[g].sum()==k for g,k in groups) and hashlib.sha256(mask.tobytes()).hexdigest()==row.partition_sha256)
 put(pd.DataFrame(checks),'tests/STRICT_ADDENDUM_TESTS.tsv');assert all(x['status']=='PASS' for x in checks)
 summary={'strict_checks_passed':len(checks),'earlier_resolver_failures':4,'earlier_failed_results_retained':True,'max_absolute_score_change':float(pd.concat(checkscores).score_change.abs().max()),'common_P_changed_rows':int((abs(j.p_strict-j.p_initial)>1e-12).sum()),'historical_P_changed_rows':int((abs(h.p_strict-h.p_initial)>1e-12).sum()),'no_Holm_below_05':bool((be.holm_p>=.05).all() if be.holm_p.notna().all() else not(be.holm_p<.05).any()),'minimum_Holm':float(be.holm_p.min()),'strict_lock':json.loads((S/'config/LOCK.json').read_text())}
 for fn in ['PROXY_ADJUSTMENT_OLD_NEW.tsv','ALL_REFERENCE_PROGRAM_ASSOCIATIONS_OLD_NEW.tsv','ALL_REFERENCE_COVARIATE_SENSITIVITY_OLD_NEW.tsv']:
  f=read('consolidated_results/'+fn);summary[fn]={'q_below_05':int((f.q_new<.05).sum()),'q_below_10':int((f.q_new<.1).sum()),'minimum_q':float(f.q_new.min())}
 dump(summary,'reports/STRICT_ADDENDUM_SUMMARY.json');print(json.dumps(summary,indent=2),flush=True)
 # Archive initial presentation and failed handoff file before final display refresh.
 arc=B/'initial_presentation_before_strict_addendum';arc.mkdir(exist_ok=True)
 for folder,pattern in [('figures','Figure_R*'),('reports','REPAIR_REVIEW_REPORT.pdf'),('tests','HANDOFF_TESTS.tsv')]:
  for file in (B/folder).glob(pattern):
   dst=arc/folder/file.name;dst.parent.mkdir(parents=True,exist_ok=True)
   if not dst.exists():shutil.copyfile(file,dst)
if __name__=='__main__':run(main)
