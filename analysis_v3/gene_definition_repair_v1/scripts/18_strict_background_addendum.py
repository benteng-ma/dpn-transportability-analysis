from common import *
import pandas as pd,numpy as np,shutil
S=B/'strict_background_addendum'
def imported(path,name):
 sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def put_s(x,p):
 dst=S/p;dst.parent.mkdir(parents=True,exist_ok=True);x.to_csv(dst,sep='\t',index=False)
def dump_s(x,p):
 dst=S/p;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_text(json.dumps(x,indent=2,default=str),encoding='utf-8')
def main():
 assert not (S/'config/LOCK.json').exists(),'Strict branch already locked; do not overwrite'
 for d in ['inputs/matrices','config','results','reports']:(S/d).mkdir(parents=True,exist_ok=True)
 m12=load('12_validate_hdrg_components_in_human_sural_nerve.py');m06=load('06_validate_hdrg_signatures_in_independent_human_bulk.py')
 info=m12.load_gene_info(m12.NCBI/'Homo_sapiens.gene_info.gz');look=m12.build_symbol_lookup(info)
 dup=set(info.loc[info.Symbol.duplicated(keep=False),'Symbol']);look=( {k:v for k,v in look[0].items() if k not in dup},look[1],{k:v for k,v in look[2].items() if k not in {d.upper() for d in dup}},look[3])
 put_s(info[info.Symbol.isin(dup)],'inputs/AMBIGUOUS_OFFICIAL_SYMBOLS.tsv')
 inventories=[];delta=[]
 for z in pd.read_csv(B/'inputs/MATRIX_INVENTORY.tsv',sep='\t').query("background=='statfree' and endpoint!='GSE24290_C01'").itertuples():
  eid=z.endpoint;mf=B/'inputs/matrices'/f'{eid}_metadata.tsv';meta=pd.read_csv(mf,sep='\t');samples=meta.sample_id.tolist();shutil.copyfile(mf,S/'inputs/matrices'/mf.name)
  if eid=='original_hDRG_DPN_control':raw=pd.read_excel(m06.RAW/'41598_2022_8100_MOESM6_ESM.xlsx',sheet_name='S5 Normalized_data D_vs_C')
  elif eid.startswith('original_JCI'):
   case=eid.endswith('DPN_control');raw=pd.read_excel(m12.FILES['dpn_control' if case else 'axon_severity'],sheet_name='deseq2_results' if case else 'Deseq_mod_vs_severe_axonLoss');raw.columns=raw.columns.astype(str).str.strip()
  else:raw=pd.read_csv(P/'analysis_v3/data/processed/GSE148059_rlog_expression.tsv.gz',sep='\t')
  raw=raw.rename(columns={raw.columns[0]:'raw_gene'});mapped=m12.add_resolution(raw[['raw_gene',*samples]],'raw_gene',look)
  for col in samples:mapped[col]=pd.to_numeric(mapped[col],errors='coerce')
  mapped['include_row']=mapped.human_gene_id.notna()&np.isfinite(mapped[samples]).all(axis=1)
  put_s(mapped.drop(columns=samples),f'inputs/{eid}_STRICT_ROW_MAPPING.tsv')
  x=mapped[mapped.include_row].groupby('human_gene_id',sort=True)[samples].median()
  if eid!='GSE148059_C02':x=x.loc[(x!=0).any(axis=1)]
  x.index=x.index.astype(str);out=S/'inputs/matrices'/f'{eid}__statfree.tsv.gz';x.to_csv(out,sep='\t',index_label='human_gene_id')
  old=pd.read_csv(B/z.matrix,sep='\t',index_col=0);old.index=old.index.astype(str)
  common=old.index.intersection(x.index);delta.append(dict(endpoint=eid,old_genes=len(old),strict_genes=len(x),removed_ids='|'.join(sorted(set(old.index)-set(x.index))),added_ids='|'.join(sorted(set(x.index)-set(old.index))),common_value_max_difference=float(abs(x.loc[common]-old.loc[common]).to_numpy().max())))
  inventories.append(dict(endpoint=eid,background='statfree',n_genes=len(x),n_samples=len(samples),matrix=out.relative_to(S).as_posix(),sha256=sha(out)))
 for name in ['REPAIRED_MEMBERS.tsv','REPAIRED_PARTITIONS.tsv']:shutil.copyfile(B/'inputs'/name,S/'inputs'/name)
 put_s(pd.DataFrame(inventories),'inputs/MATRIX_INVENTORY.tsv');put_s(pd.DataFrame(delta),'reports/BACKGROUND_CHANGE.tsv')
 files=[f for f in (S/'inputs').rglob('*') if f.is_file()];put_s(pd.DataFrame([dict(file=f.relative_to(S).as_posix(),sha256=sha(f)) for f in files]),'config/LOCKED_FILES.tsv')
 dump_s({'locked_utc':now(),'nature':'Post-result implementation correction; not preregistration','rule':'Remove ambiguous official symbols rather than use the legacy last-row dictionary; no target-P/effect decision. Reaggregate actual finite expression rows, median per GeneID. Source members unchanged. Only four affected stat-free backgrounds recomputed; GSE24290 unchanged.','source_members_sha256':sha(B/'inputs/REPAIRED_MEMBERS.tsv'),'input_manifest_sha256':sha(S/'config/LOCKED_FILES.tsv'),'implementation_sha256':sha(__file__),'strict_associations_viewed':False},'config/LOCK.json')
 engine=imported(B/'scripts/02_associations_benchmarks.py','strict_engine');engine.B=S;engine.put=put_s;engine.dump=dump_s;engine.main()
 # Preserve engine-only outputs and build a complete view without overwriting any initial outputs.
 for fn in ['ASSOCIATION_IMPACT_ALL.tsv','COVERAGE.tsv','DONOR_SCORES.tsv.gz','REAL_BENCHMARK_ALL.tsv','MATCHED_RANDOM_REFERENCE.tsv','MATCHED_RANDOM_SPLITS.tsv.gz','GATES.tsv','ALGEBRA_RECONSTRUCTION.tsv']:
  old=pd.read_csv(B/'results'/fn,sep='\t');new=pd.read_csv(S/'results'/fn,sep='\t');col='version' if 'version' in old else 'background';val='repaired_statfree' if col=='version' else 'statfree'
  keep=~(old[col].eq(val)&old.endpoint.isin(new.endpoint));merged=pd.concat([old[keep],new],ignore_index=True);put_s(merged,'consolidated_results/'+fn)
 print('strict background main and benchmark completed',flush=True)
if __name__=='__main__':run(main)
