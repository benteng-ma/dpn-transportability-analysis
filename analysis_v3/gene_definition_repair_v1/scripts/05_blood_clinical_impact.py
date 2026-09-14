from common import *
import pandas as pd,numpy as np

def main():
 mods=pd.read_csv(B/'inputs/REPAIRED_MEMBERS.tsv',sep='\t',dtype={'gene_id':str});raw=pd.read_csv(B/'inputs/SOURCE_ALL_ROWS.tsv',sep='\t',dtype={'gene_id':str});c=pd.read_csv(B/'inputs/CANONICAL_SOURCE_COMPARISONS.tsv',sep='\t',dtype={'gene_id':str})
 f=mods[mods.module_id.str.startswith('original_')].merge(c[['contrast','gene_id','median_effect']],left_on=['source_contrast','gene_id'],right_on=['contrast','gene_id'],validate='many_to_one')
 q=raw.groupby(['contrast','gene_id']).q.max();f['p_val_adj']=[q.loc[(r.contrast,r.gene_id)] for r in f.itertuples()];f['primary_signature_member']=True;f=f.rename(columns={'symbol':'gene','median_effect':'avg_log2FC','contrast':'contrast_id'})
 put(f,'inputs/REPAIRED_SIGNATURE_ADAPTER.tsv')
 dump(dict(note='p_val_adj in compatibility input is maximum deposited source-row adjusted P, not a newly calculated gene-level P. Membership is already locked; compatibility readers cannot reselect genes.',members_sha256=sha(B/'inputs/REPAIRED_MEMBERS.tsv'),adapter_sha256=sha(B/'inputs/REPAIRED_SIGNATURE_ADAPTER.tsv')),'reports/SIGNATURE_ADAPTER_AUDIT.json')
 # Read-only absolute inputs remain original; ALL declared output directories are redirected locally.
 for task,filename in [('PBMC','08_validate_hdrg_stages_in_human_pbmc_cohorts.py'),('GSE302658','10_validate_hdrg_severity_in_GSE302658.py')]:
  m=load(filename);out=B/'results/historical'/task;out.mkdir(parents=True,exist_ok=True)
  m.TABLES=out;m.FIGURES=B/'figures'/task;m.METADATA=B/'reports'/task;m.METADATA.mkdir(parents=True,exist_ok=True);m.FIGURES.mkdir(parents=True,exist_ok=True);m.SIGNATURES=B/'inputs/REPAIRED_SIGNATURE_ADAPTER.tsv'
  # No source-screening, external actions, or uncontrolled path rebindings occur.
  for name in ['TABLES','FIGURES','METADATA']:assert getattr(m,name).resolve().is_relative_to(B.resolve())
  dump(dict(task=task,source_script=filename,sha256=sha(P/'analysis/scripts'/filename),outputs=[getattr(m,n).relative_to(B).as_posix() for n in ['TABLES','FIGURES','METADATA']],nature='Historical algorithms rerun with locked corrected membership; not independent validation'),'logs/'+task+'_ADAPTER_START.json')
  m.main();dump(dict(task=task,exit_code=0,ended_utc=now()),'logs/'+task+'_ADAPTER_END.json')
  print(task,'completed',flush=True)
if __name__=='__main__':run(main)
