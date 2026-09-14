from common import *
import pandas as pd,numpy as np,shutil,itertools,math
from scipy import stats
from statsmodels.stats.multitest import multipletests
# Use independently tested numerical functions without executing old scripts or writing to old paths.
s=importlib.util.spec_from_file_location('fixed_math',B4/'scripts/benchmark_math.py');bm=importlib.util.module_from_spec(s);s.loader.exec_module(bm)
FAMILIES=[('late_allcell','original_late_allcell','late_shared_concordant_neuronal_core','late_allcell_residual'),('late_neuron','original_late_neuron','late_shared_concordant_neuronal_core','late_neuron_residual'),('severity','original_severity','severity_neuron_shared_concordant_core','severity_neuron_residual')]

def bootstrap(scores,y,sex,rng):
 strata=[np.flatnonzero((y==g)&(sex==ss)) for g in sorted(set(y)) for ss in sorted(set(sex)) if np.any((y==g)&(sex==ss))]
 idx=np.concatenate([rng.choice(ix,size=(2000,len(ix)),replace=True) for ix in strata],axis=1)
 by=np.concatenate([y[ix] for ix in strata]);v=scores[idx]
 d=v[:,by==1].mean(axis=1)-v[:,by==0].mean(axis=1)
 return np.quantile(d,[.025,.975],axis=0)

def main():
 for z in pd.read_csv(B/'config/LOCKED_FILES.tsv',sep='\t').itertuples():assert sha(B/z.file)==z.sha256,z.file
 repaired=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str})
 part=pd.read_csv(B/'inputs/REPAIRED_PARTITIONS.tsv',sep='\t',dtype={'gene_id':str})
 old=pd.read_csv(B4/'inputs/SOURCE_CANONICAL_MAPPING.tsv',sep='\t',dtype={'human_gene_id':str}).rename(columns={'human_gene_id':'gene_id'})
 old=old[old.gene_id.notna()]
 allscores=[];coverage=[];assoc=[];bench=[];randoms=[];randomsummary=[];algebra=[];checks=[];permrecords=[]
 for eid,rows in pd.read_csv(B/'inputs/MATRIX_INVENTORY.tsv',sep='\t').groupby('endpoint',sort=False):
  meta=pd.read_csv(B/'inputs/matrices'/f'{eid}_metadata.tsv',sep='\t');y=meta.group_code.to_numpy();sex=meta.sex.to_numpy() if eid.startswith('original_hDRG') else np.repeat('all',len(y))
  alloc,exact,space=bm.label_allocations(y,np.random.default_rng(seed(eid,'permutation')),count=9999,exact_max=100000,strata=sex)
  permrecords.append(dict(endpoint=eid,exact=exact,allocation_space=str(space),draws=len(alloc),seed=str(seed(eid,'permutation')),allocation_sha256=hashlib.sha256(alloc.tobytes()).hexdigest()))
  models={}
  for row in rows.itertuples():
   x=pd.read_csv(B/row.matrix,sep='\t',index_col=0);x.index=x.index.astype(str);assert list(x.columns)==meta.sample_id.tolist()
   rank=bm.doubled_ranks(x.to_numpy());idx={g:i for i,g in enumerate(x.index)}
   for version,source in [('legacy',old),('repaired',repaired)]:
    # Do not add an unrequested fourth branch: historical legacy background only.
    if version=='legacy' and row.background!='inherited':continue
    label=version+'_'+row.background;sc={}
    for mod,f in source.groupby('module_id',sort=True):
     if version=='legacy':
      pair=eid.startswith('GSE') or (eid.startswith('original_hDRG') and mod.startswith('original_'))
      f=f.drop_duplicates(['gene_id','direction'] if pair else ['gene_id'])
     fm=f[f.gene_id.isin(idx)];up=[idx[g] for g in fm.loc[fm.direction=='up','gene_id']];down=[idx[g] for g in fm.loc[fm.direction=='down','gene_id']]
     v=bm.signed_score(rank,up,down);sc[mod]=v
     coverage.append(dict(endpoint=eid,version=label,module_id=mod,source_members=len(f),measured=len(fm),n_up=len(up),n_down=len(down),coverage=len(fm)/len(f),background_genes=len(x),eligible=bool(np.isfinite(v).all())))
     for j,sample in enumerate(meta.sample_id):allscores.append(dict(endpoint=eid,version=label,module_id=mod,sample_id=sample,patient_id=meta.patient_id.iloc[j],group=meta.group.iloc[j],score=v[j]))
     if not np.isfinite(v).all():assoc.append(dict(endpoint=eid,version=label,module_id=mod,status='NOT_EVALUABLE',reason='fewer than 10 measured genes per signed arm',p=np.nan));continue
     d,t=bm.all_stats(v[:,None],y,(y==1)[None,:].astype(float));_,null=bm.all_stats(v[:,None],y,alloc)
     p=float(bm.perm_p(t[0],null,exact)[0]);ci=bootstrap(v[:,None],y,sex,np.random.default_rng(seed(eid,mod,'bootstrap')))
     pos=v[y==1];neg=v[y==0];df=len(pos)+len(neg)-2;pooled=((len(pos)-1)*pos.var(ddof=1)+(len(neg)-1)*neg.var(ddof=1))/df
     auc=stats.mannwhitneyu(pos,neg,alternative='two-sided').statistic/(len(pos)*len(neg))
     assoc.append(dict(endpoint=eid,version=label,module_id=mod,status='COMPLETED',n_total=len(y),n_positive=len(pos),n_negative=len(neg),mean_difference=d[0,0],ci_low=ci[0,0],ci_high=ci[1,0],t=t[0,0],p=p,hedges_g=d[0,0]/np.sqrt(pooled)*(1-3/(4*df-1)),fixed_direction_auc=auc,missing=0,exact=exact,permutations=len(alloc)))
    models[label]=sc
   if row.background in ['inherited','statfree']:
    for fam,parent,core,residual in FAMILIES:
     f=part[(part.family==fam)&part.gene_id.isin(idx)].copy();sets={n:set(f.loc[f.bin.eq(n),'gene_id']) for n in ['core','residual']};e=set(f.gene_id)-sets['core']-sets['residual']
     assert not sets['core']&sets['residual'] and len(set(f.gene_id))==len(f)
     coremem=repaired[repaired.module_id==core].set_index('gene_id').direction.to_dict();resmem=repaired[repaired.module_id==residual].set_index('gene_id').direction.to_dict()
     assert sets['core']==set(coremem)&set(idx) and sets['residual']==set(resmem)&set(idx)
     assert all(coremem[g]==d for g,d in zip(f.loc[f.bin=='core','gene_id'],f.loc[f.bin=='core','direction']))
     assert all(resmem[g]==d for g,d in zip(f.loc[f.bin=='residual','gene_id'],f.loc[f.bin=='residual','direction']))
     full=models['repaired_'+row.background][parent];c=models['repaired_'+row.background][core];r=models['repaired_'+row.background][residual]
     union=f[f.bin.isin(['core','residual'])];u=[idx[g] for g in union.loc[union.direction=='up','gene_id']];dn=[idx[g] for g in union.loc[union.direction=='down','gene_id']];restricted=bm.signed_score(rank,u,dn)
     pu=f.direction.eq('up').sum();pdn=f.direction.eq('down').sum();terms=[]
     for ss in [sets['core'],sets['residual'],e]:
      sf=f[f.gene_id.isin(ss)];terms.append(rank[[idx[g] for g in sf.loc[sf.direction=='up','gene_id']]].sum(axis=0)/(2*len(x)*pu)-rank[[idx[g] for g in sf.loc[sf.direction=='down','gene_id']]].sum(axis=0)/(2*len(x)*pdn))
     err=float(np.max(abs(sum(terms)-full)));assert err<1e-12
     checks.append(dict(endpoint=eid,background=row.background,family=fam,gate='PASS',max_algebra_error=err))
     for j,sample in enumerate(meta.sample_id):algebra.append(dict(endpoint=eid,background=row.background,family=fam,sample_id=sample,parent_score=full[j],core_contribution=terms[0][j],residual_contribution=terms[1][j],excluded_contribution=terms[2][j],restricted_score=restricted[j]))
     vv=np.c_[full,restricted,c,r,c-r]
     if not np.isfinite(vv).all():
      for test in ['parent_full','parent_restricted','split_joint_maxT','core_minus_residual']:bench.append(dict(endpoint=eid,background=row.background,family=fam,test=test,status='NOT_EVALUABLE',p=np.nan,reason='arm coverage'))
      continue
     effect,obs=bm.all_stats(vv,y,(y==1)[None,:].astype(float));_,null=bm.all_stats(vv,y,alloc)
     cis=bootstrap(vv,y,sex,np.random.default_rng(seed(eid,fam,'benchmark_bootstrap')))
     for test,col in [('parent_full',0),('parent_restricted',1),('split_joint_maxT',None),('core_minus_residual',4)]:
      ot=np.max(abs(obs[0,2:4])) if col is None else abs(obs[0,col]);nt=np.max(abs(null[:,2:4]),axis=1) if col is None else abs(null[:,col]);pp=float(bm.perm_p(ot,nt,exact))
      bench.append(dict(endpoint=eid,background=row.background,family=fam,test=test,status='COMPLETED',statistic=ot,p=pp,effect=np.nan if col is None else effect[0,col],ci_low=np.nan if col is None else cis[0,col],ci_high=np.nan if col is None else cis[1,col],joint_core_effect=effect[0,2] if col is None else np.nan,joint_residual_effect=effect[0,3] if col is None else np.nan,n=len(y)))
     union=union.sort_values('gene_id');pos=np.array([idx[g] for g in union.gene_id]);sign=np.where(union.direction=='up',1,-1);cm=union.bin.eq('core').to_numpy();groups,strata,mean,mad=bm.matching_strata(x.iloc[pos].to_numpy(),sign,cm)
     space=math.prod(math.comb(len(g),k) for g,k in groups);enumerated=space<=999;rng=np.random.default_rng(seed(eid,fam,row.background,'matched_split'));seen=set();values=[]
     if enumerated:
      options=[list(itertools.combinations(g,k)) for g,k in groups]
      def gen():
       for choices in itertools.product(*options):
        a=np.zeros(len(union),bool);a[[v for c0 in choices for v in c0]]=True;yield a
      draws=gen()
     else:draws=(bm.random_split(groups,rng,len(union)) for _ in range(999))
     for j,mask in enumerate(draws):
      assert all(mask[g].sum()==k for g,k in groups)
      a=bm.signed_score(rank,pos[mask&(sign>0)],pos[mask&(sign<0)]);b=bm.signed_score(rank,pos[~mask&(sign>0)],pos[~mask&(sign<0)])
      _,t0=bm.all_stats(np.c_[a,b],y,(y==1)[None,:].astype(float));val=float(np.max(abs(t0)));values.append(val);sig=hashlib.sha256(mask.tobytes()).hexdigest();seen.add(sig)
      randoms.append(dict(endpoint=eid,background=row.background,family=fam,draw=j,max_abs_t=val,partition_sha256=sig,core_gene_ids='|'.join(union.gene_id.to_numpy()[mask]),mismatch_count=0))
     target=float(np.max(abs(obs[0,2:4])));values=np.array(values)
     randomsummary.append(dict(endpoint=eid,background=row.background,family=fam,source_max_abs_t=target,draws=len(values),unique_draws=len(seen),space=str(space),enumerated=enumerated,conditional_midrank_percentile=float((np.sum(values<target-1e-12)+.5*np.sum(abs(values-target)<=1e-12))/len(values)),not_validation_p=True))
  print('completed',eid,flush=True)
 a=pd.DataFrame(assoc);a['q_BH_common_method']=np.nan
 for _,ix in a.groupby(['endpoint','version']).groups.items():
  ix=[i for i in ix if np.isfinite(a.loc[i,'p'])];a.loc[ix,'q_BH_common_method']=multipletests(a.loc[ix,'p'],method='fdr_bh')[1] if ix else []
 be=pd.DataFrame(bench);be['holm_p']=np.nan
 for _,ix in be.groupby('background').groups.items():
  ix=[i for i in ix if np.isfinite(be.loc[i,'p'])];be.loc[ix,'holm_p']=multipletests(be.loc[ix,'p'],method='holm')[1] if ix else []
 for df,name in [(pd.DataFrame(allscores),'DONOR_SCORES.tsv.gz'),(pd.DataFrame(coverage),'COVERAGE.tsv'),(a,'ASSOCIATION_IMPACT_ALL.tsv'),(be,'REAL_BENCHMARK_ALL.tsv'),(pd.DataFrame(randoms),'MATCHED_RANDOM_SPLITS.tsv.gz'),(pd.DataFrame(randomsummary),'MATCHED_RANDOM_REFERENCE.tsv'),(pd.DataFrame(algebra),'ALGEBRA_RECONSTRUCTION.tsv'),(pd.DataFrame(checks),'GATES.tsv'),(pd.DataFrame(permrecords),'PERMUTATION_RECORDS.tsv')]:put(df,'results/'+name)
 # Wide old/new impact tables preserve like-for-like methods, not substituted historical P/q.
 a.set_index(['endpoint','module_id','version']).unstack('version').to_csv(B/'results/ASSOCIATION_IMPACT_WIDE.tsv',sep='\t')
 dump(dict(association_rows=len(a),benchmark_rows=len(be),random_splits=len(randoms),gates=len(checks),legacy_b4_scores_sha256=sha(B4/'results/DONOR_PROGRAM_SCORES_AUDIT.tsv.gz'),note='Old/new association P uses the same new retrospective test. Historical P/q remains in protected outputs and separate dependency register.'),'reports/ANALYSIS_RUN_SUMMARY.json')
if __name__=='__main__':run(main)
