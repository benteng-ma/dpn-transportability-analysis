from common import *
import pandas as pd,numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def main():
 checks=[]
 def check(name,value,detail=''):checks.append(dict(check=name,status='PASS' if value else 'FAIL',detail=str(detail)))
 for z in pd.read_csv(B/'config/LOCKED_FILES.tsv',sep='\t').itertuples():check('locked:'+z.file,sha(B/z.file)==z.sha256)
 for z in pd.read_csv(B/'inputs/PROTECTED_FILES.tsv',sep='\t').itertuples():check('protected:'+z.file,sha(P/z.file)==z.sha256)
 src=pd.read_csv(B/'inputs/SOURCE_ALL_ROWS.tsv',sep='\t',dtype={'gene_id':str}).fillna({'gene_id':''});can=pd.read_csv(B/'inputs/CANONICAL_SOURCE_COMPARISONS.tsv',sep='\t',dtype={'gene_id':str});members=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str});parts=pd.read_csv(B/'inputs/REPAIRED_PARTITIONS.tsv',sep='\t',dtype={'gene_id':str})
 check('unique_standard_gene_per_program',not members.duplicated(['module_id','gene_id']).any())
 check('no_signed_conflicts_in_repaired_members',set(members.direction)=={'up','down'})
 for row in can.itertuples():
  f=src[(src.contrast==row.contrast)&(src.gene_id==row.gene_id)];signs=np.sign(f.effect.to_numpy());ok=np.isfinite(f.effect).all() and len(set(signs))==1 and 0 not in signs and ((f.q<.05)&(abs(f.effect)>.585)).all()
  check('source_rule:'+row.contrast+':'+row.gene_id,(row.status=='eligible')==bool(ok))
 for fam,f in parts.groupby('family'):
  par=members[members.module_id==f.parent.iloc[0]].set_index('gene_id').direction
  check('partition_parent:'+fam,par.sort_index().equals(f.set_index('gene_id').direction.sort_index()))
  check('partition_disjoint:'+fam,f.gene_id.is_unique)
  cross=f[f.bin=='excluded_cross_source_opposed'];check('cross_conflicts_retained:'+fam,(cross.direction!=cross.partner_direction).all())
 scores=pd.read_csv(B/'results/DONOR_SCORES.tsv.gz',sep='\t');legacy=pd.read_csv(B4/'results/DONOR_PROGRAM_SCORES_AUDIT.tsv.gz',sep='\t')
 j=scores[scores.version=='legacy_inherited'].merge(legacy,on=['endpoint','sample_id'],suffixes=('_new','_old'));j=j[j.module_id==j.program]
 check('1530_legacy_score_parity',len(j)==1530 and np.allclose(j.score_new,j.score_old,atol=1e-12,equal_nan=True));put(j[['endpoint','sample_id','module_id','score_new','score_old']],'tests/LEGACY_SCORE_PARITY.tsv')
 h=pd.read_csv(B/'results/HISTORICAL_DONOR_TEST_IMPACT.tsv',sep='\t');hh=h[h.version=='legacy_inherited'];check('historical_original_P_reproduced',np.allclose(hh.p,hh.archived_p,atol=1e-12,equal_nan=True));check('historical_original_Q_reproduced',np.allclose(hh.q,hh.archived_q,atol=1e-12,equal_nan=True))
 numeric=[]
 for row in pd.read_csv(B/'inputs/MATRIX_INVENTORY.tsv',sep='\t').itertuples():
  x=pd.read_csv(B/row.matrix,sep='\t',index_col=0);x.index=x.index.astype(str);r=x.rank(axis=0,method='average',pct=True)
  for mod,f in members.groupby('module_id'):
   u=list(set(f.loc[f.direction=='up','gene_id'])&set(x.index));d=list(set(f.loc[f.direction=='down','gene_id'])&set(x.index));want=r.loc[u].mean()-r.loc[d].mean() if min(len(u),len(d))>=10 else pd.Series(np.nan,index=x.columns)
   got=scores[(scores.endpoint==row.endpoint)&(scores.version=='repaired_'+row.background)&(scores.module_id==mod)].set_index('sample_id').loc[x.columns].score
   check('independent_rank:'+row.endpoint+':'+row.background+':'+mod,np.allclose(got,want,equal_nan=True,atol=1e-12))
   numeric.append(dict(endpoint=row.endpoint,background=row.background,module_id=mod,max_abs_difference=float(abs(got-want).max()) if got.notna().any() else np.nan))
  old=pd.read_csv(B/'inputs/matrices'/f'{row.endpoint}__inherited.tsv.gz',sep='\t',index_col=0);old.index=old.index.astype(str);common=old.index.intersection(x.index)
  numeric.append(dict(endpoint=row.endpoint,background=row.background,module_id='BACKGROUND_VALUE_COMPARISON',max_abs_difference=float(abs(old.loc[common]-x.loc[common]).to_numpy().max()),genes_added=len(set(x.index)-set(old.index)),genes_removed=len(set(old.index)-set(x.index))))
 put(pd.DataFrame(numeric),'tests/INDEPENDENT_RANK_AND_BACKGROUND_CHECK.tsv')
 g=pd.read_csv(B/'results/GATES.tsv',sep='\t');check('all_30_measured_partition_checks',len(g)==30 and (g.max_algebra_error<1e-12).all())
 b=pd.read_csv(B/'results/REAL_BENCHMARK_ALL.tsv',sep='\t');check('120_tasks_retained',len(b)==120);check('severity_not_forced_past_coverage',b[b.family=='severity'].status.eq('NOT_EVALUABLE').all());check('80_late_tests_completed',len(b[b.status=='COMPLETED'])==80)
 for bg,f in b[b.status=='COMPLETED'].groupby('background'):check('holm:'+bg,np.allclose(f.holm_p,multipletests(f.p,method='holm')[1]))
 rs=pd.read_csv(B/'results/MATCHED_RANDOM_SPLITS.tsv.gz',sep='\t');check('matched_counts',len(rs)==19980 and rs.mismatch_count.eq(0).all())
 # Side-by-side complete historical files: retain missing rather than force a match.
 pairs=[('results/historical/PBMC/human_PBMC_stage_projection_tests_2026-08-27.tsv','results/tables/human_PBMC_stage_projection_tests_2026-08-27.tsv',['cohort','contrast_id','transition'],'PBMC'),('results/historical/GSE302658/GSE302658_clinical_signature_tests_2026-08-27.tsv','results/tables/GSE302658_clinical_signature_tests_2026-08-27.tsv',['contrast_id','test_family','window','endpoint'],'CLINICAL')]
 for new,old,keys,name in pairs:
  nn=pd.read_csv(B/new,sep='\t');oo=pd.read_csv(P/old,sep='\t');out=nn.merge(oo,on=keys,how='outer',suffixes=('_new','_old'),validate='one_to_one',indicator=True);check(name+'_all_rows_joined',out._merge.eq('both').all());put(out,'results/'+name+'_OLD_NEW.tsv')
 oldatlas=pd.read_csv(P/'results/tables/cross_target_hDRG_component_transportability_tests_2026-08-27.tsv',sep='\t')
 newatlas=pd.concat([pd.read_csv(B/f'results/historical/{name}_0.tsv.gz',sep='\t') for name in ['human_components','rat_DRG_components','ocular_components']])
 put(newatlas.merge(oldatlas,on=['target_id','module_id'],suffixes=('_new','_old'),validate='one_to_one'),'results/CROSS_TISSUE_COMPONENT_OLD_NEW.tsv')
 # Source-derived candidates are annotations, not new gene screening: retain all 405 original records.
 old=pd.read_csv(P/'analysis_v3/batch02/results/concordance/V2_405_CANDIDATES_PRESERVED.tsv',sep='\t',dtype={'human_gene_id':str});old['canonical_membership_still_present']=[(r.module_id,str(r.human_gene_id)) in set(zip(members.module_id,members.gene_id)) for r in old.itertuples()];put(old,'results/ORIGINAL_405_MEMBERSHIP_IMPACT.tsv');check('405_records_preserved',len(old)==405)
 # Do not misreport target-dependent aliases as new biological replication.
 check('source_conflicts_excluded_from_signed_parents',all(not set(f[f.status=='within_source_direction_conflict'].gene_id)&set(members[(members.module_id.str.startswith('original_'))&(members.source_contrast==k)].gene_id) for k,f in can.groupby('contrast')))
 put(pd.DataFrame(checks),'tests/ALL_TESTS.tsv');summary=pd.Series([c['status'] for c in checks]).value_counts().to_dict();dump(summary,'tests/TEST_SUMMARY.json');print(summary,flush=True)
 assert summary.get('FAIL',0)==0
if __name__=='__main__':run(main)
