from common import *
import pandas as pd,numpy as np

def main():
 checks=[]
 def ck(name,v,detail=''):checks.append(dict(check=name,status='PASS' if v else 'FAIL',detail=str(detail)))
 members=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str})
 # Independently check every stored random partition's signed stratum counts.
 sp=importlib.util.spec_from_file_location('assoc',B/'scripts/02_associations_benchmarks.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a)
 parts=pd.read_csv(B/'inputs/REPAIRED_PARTITIONS.tsv',sep='\t',dtype={'gene_id':str});r=pd.read_csv(B/'results/MATCHED_RANDOM_SPLITS.tsv.gz',sep='\t')
 for (eid,bg,fam),v in r.groupby(['endpoint','background','family']):
  x=pd.read_csv(B/'inputs/matrices'/f'{eid}__{bg}.tsv.gz',sep='\t',index_col=0);x.index=x.index.astype(str)
  u=parts[(parts.family==fam)&parts.bin.isin(['core','residual'])&parts.gene_id.isin(x.index)].sort_values('gene_id')
  sign=np.where(u.direction=='up',1,-1);cm=u.bin.eq('core').to_numpy();groups,_,_,_=a.bm.matching_strata(x.loc[u.gene_id].to_numpy(),sign,cm)
  for row in v.itertuples():
   ids=set(row.core_gene_ids.split('|'));mask=u.gene_id.isin(ids).to_numpy()
   ck(f'random_partition:{eid}:{bg}:{fam}:{row.draw}',ids<=set(u.gene_id) and all(mask[g].sum()==k for g,k in groups) and hashlib.sha256(mask.tobytes()).hexdigest()==row.partition_sha256)
 src=pd.read_csv(B/'inputs/SOURCE_ALL_ROWS.tsv',sep='\t')
 duplicates=int(src.duplicated(['contrast','raw_gene','effect','q']).sum());ck('exact_source_rows_duplicate_audit',True,duplicates)
 # Preserve the old spatial gate, now with corrected GeneIDs and the archived panel.
 info=pd.read_csv(P/'data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz',sep='\t',dtype=str,keep_default_na=False)
 symbols=dict(zip(info.Symbol,info.GeneID));panel=pd.read_csv(P/'analysis_v3/batch03/inputs/SPATIAL_99_GENE_PANEL.tsv',sep='\t');ids={symbols[g] for g in panel.gene if g in symbols}
 rows=[]
 for name,f in members.groupby('module_id'):
  matched=f[f.gene_id.isin(ids)];nu=matched.direction.eq('up').sum();nd=matched.direction.eq('down').sum()
  rows.append(dict(object=name,source_members=len(f),panel_members=len(matched),up=nu,down=nd,coverage=len(matched)/len(f),whole_set_eligible=nu>=10 and nd>=10 and len(matched)/len(f)>=.2))
 put(pd.DataFrame(rows),'results/SPATIAL_REPAIRED_PROGRAM_COVERAGE.tsv')
 ck('spatial_all_99_symbols_resolved',len(ids)==99)
 # Inspect all secondary result families without selecting a favorable subset.
 summaries={}
 for file,col in [('REGULATOR_EXCLUSION_OLD_NEW.tsv','q_new'),('PROXY_ADJUSTMENT_OLD_NEW.tsv','q_new'),('ALL_REFERENCE_PROGRAM_ASSOCIATIONS_OLD_NEW.tsv','q_new'),('ALL_REFERENCE_COVARIATE_SENSITIVITY_OLD_NEW.tsv','q_new')]:
  f=pd.read_csv(B/'results'/file,sep='\t');summaries[file]=dict(rows=len(f),finite_q=int(f[col].notna().sum()),q_below_05=int((f[col]<.05).sum()),q_below_10=int((f[col]<.1).sum()),minimum_q=float(f[col].min()))
 ck('canonical_lock_still_intact',sha(B/'inputs/REPAIRED_MEMBERS.tsv')==json.loads((B/'config/LOCK.json').read_text())['members_sha256'])
 for row in pd.read_csv(B/'inputs/PROTECTED_FILES.tsv',sep='\t').itertuples():ck('final_protected:'+row.file,sha(P/row.file)==row.sha256)
 put(pd.DataFrame(checks),'tests/ADDITIONAL_TESTS.tsv');dump({'counts':pd.Series([z['status'] for z in checks]).value_counts().to_dict(),'source_exact_duplicates':duplicates,'secondary_result_families':summaries,'spatial_evaluable_programs':sum(r['whole_set_eligible'] for r in rows)},'reports/FINAL_AUDIT.json')
 print(json.dumps(summaries,indent=2));print('additional checks',len(checks),flush=True)
 assert all(z['status']=='PASS' for z in checks)
if __name__=='__main__':run(main)
