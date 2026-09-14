from common import *
import pandas as pd,numpy as np,shutil

def main():
 assert not (B/'config/LOCK.json').exists()
 m12=load('12_validate_hdrg_components_in_human_sural_nerve.py');m06=load('06_validate_hdrg_signatures_in_independent_human_bulk.py')
 look=m12.build_symbol_lookup(m12.load_gene_info(m12.NCBI/'Homo_sapiens.gene_info.gz'))
 inventory=[];audits=[];src=[]
 for z in pd.read_csv(B4/'results/INPUT_AVAILABILITY.tsv',sep='\t').itertuples():
  eid=z.endpoint;folder=B/'inputs/matrices';folder.mkdir(exist_ok=True)
  mf=B4/'inputs/matrices'/f'{eid}_metadata.tsv';meta=pd.read_csv(mf,sep='\t');samples=meta.sample_id.tolist()
  shutil.copyfile(mf,folder/mf.name)
  old=B4/'inputs/matrices'/f'{eid}.tsv.gz';shutil.copyfile(old,folder/f'{eid}__inherited.tsv.gz')
  src.extend([mf,old]); raw=None;sheet=''
  if eid=='original_hDRG_DPN_control':
   path=m06.RAW/'41598_2022_8100_MOESM6_ESM.xlsx';sheet='S5 Normalized_data D_vs_C';raw=pd.read_excel(path,sheet_name=sheet)
  elif eid.startswith('original_JCI'):
   path=m12.FILES['dpn_control' if eid.endswith('DPN_control') else 'axon_severity'];sheet='deseq2_results' if eid.endswith('DPN_control') else 'Deseq_mod_vs_severe_axonLoss';raw=pd.read_excel(path,sheet_name=sheet);raw.columns=raw.columns.astype(str).str.strip()
  elif eid=='GSE148059_C02':
   path=P/'analysis_v3/data/processed/GSE148059_rlog_expression.tsv.gz';raw=pd.read_csv(path,sep='\t')
  else:
   path=old;x=pd.read_csv(old,sep='\t',index_col=0);x.index=x.index.astype(str)
  if raw is not None:
   src.append(path);raw=raw.rename(columns={raw.columns[0]:'raw_gene'});raw['source_row']=np.arange(len(raw))+2
   mapped=m12.add_resolution(raw[['raw_gene','source_row',*samples]],'raw_gene',look)
   for s in samples:mapped[s]=pd.to_numeric(mapped[s],errors='coerce')
   mapped['finite_all']=np.isfinite(mapped[samples]).all(axis=1)
   mapped['include_row']=mapped.human_gene_id.notna()&mapped.finite_all
   put(mapped.drop(columns=samples),'inputs/'+eid+'_BACKGROUND_ROW_MAPPING.tsv')
   x=mapped[mapped.include_row].groupby('human_gene_id',sort=True)[samples].median()
   if eid!='GSE148059_C02':x=x.loc[(x!=0).any(axis=1)]
   audits.append(dict(endpoint=eid,source_file=path.relative_to(P).as_posix(),sheet=sheet,raw_rows=len(raw),mapped_rows=int(mapped.human_gene_id.notna().sum()),finite_mapped_rows=int(mapped.include_row.sum()),deduplicated_genes=len(x),target_stat_used=False,target_P_used=False,aggregation='per_sample_median_all_finite_mapped_rows',historic_selection='finite DE statistic and positive baseMean; JCI duplicates prioritized finite statistic, smallest adjusted P, largest absolute statistic' if eid.startswith('original_') else 'ALC original mean-expression-priority mapping'))
  else:audits.append(dict(endpoint=eid,source_file=path.relative_to(P).as_posix(),sheet='',raw_rows=len(x),deduplicated_genes=len(x),target_stat_used=False,target_P_used=False,aggregation='already canonical unique processed matrix unchanged',historic_selection='No target DE statistic used'))
  x.index=x.index.astype(str);x=x.loc[:,samples];assert x.index.is_unique and np.isfinite(x.to_numpy()).all()
  x.to_csv(folder/f'{eid}__statfree.tsv.gz',sep='\t',index_label='human_gene_id')
  for bg in ['inherited','statfree']:
   f=folder/f'{eid}__{bg}.tsv.gz';mat=pd.read_csv(f,sep='\t',index_col=0)
   inventory.append(dict(endpoint=eid,background=bg,n_genes=len(mat),n_samples=len(samples),matrix=f.relative_to(B).as_posix(),sha256=sha(f)))
 put(pd.DataFrame(inventory),'inputs/MATRIX_INVENTORY.tsv');put(pd.DataFrame(audits),'reports/BACKGROUND_SELECTION_AUDIT.tsv')
 put(pd.DataFrame([dict(file=f.relative_to(P).as_posix(),sha256=sha(f)) for f in sorted(set(src))]),'inputs/TARGET_SOURCE_FILES.tsv')
 # Freeze rules, actual members/partitions, input matrices and implementation before disease scoring.
 files=[f for root in ['inputs','config','scripts'] for f in (B/root).rglob('*') if f.is_file() and '__pycache__' not in f.parts]
 rec=[dict(file=f.relative_to(B).as_posix(),sha256=sha(f)) for f in sorted(files)]
 put(pd.DataFrame(rec),'config/LOCKED_FILES.tsv')
 dump(dict(lock_utc=now(),nature='Known-results retrospective correction, not preregistration',rules_sha256=sha(B/'config/REPAIR_RULES.md'),members_sha256=sha(B/'inputs/REPAIRED_MEMBERS.tsv'),partitions_sha256=sha(B/'inputs/REPAIRED_PARTITIONS.tsv'),manifest_sha256=sha(B/'config/LOCKED_FILES.tsv'),corrected_associations_viewed=False),'config/LOCK.json')
 print(pd.DataFrame(inventory)[['endpoint','background','n_genes','n_samples']].to_string(index=False),flush=True)
if __name__=='__main__':run(main)
