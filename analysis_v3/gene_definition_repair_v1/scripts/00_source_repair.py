from common import *
import pandas as pd,numpy as np,shutil,platform
from collections import defaultdict

def main():
 assert P.name=='phase0_6_human_dpn_stage_projection'
 assert not (B/'config/LOCK.json').exists(),'Locked versions are immutable'
 for d in ['inputs','results','reports','tests','logs','config','figures']:(B/d).mkdir(exist_ok=True)
 assert sha(B4/'DPN_v3_Batch04_REVIEW_BUNDLE_2026-09-10_v1.zip')=='47e716e8048035662c8978658f4ac38c8ed31d7c9374382110bafc0717b39b36'
 protected=pd.read_csv(B4/'inputs_manifest/PROTECTED_BASELINE_SHA256.tsv',sep='\t')
 for r in protected.itertuples():assert sha(P/r.file)==r.sha256
 extra=[dict(file=f.relative_to(P).as_posix(),sha256=sha(f)) for f in B4.rglob('*') if f.is_file() and '__pycache__' not in f.parts and 'environment' not in f.relative_to(B4).parts and 'rendered' not in f.parts]
 put(pd.concat([protected[['file','sha256']],pd.DataFrame(extra)]).drop_duplicates('file'),'inputs/PROTECTED_FILES.tsv')
 m03=load('03_audit_and_extract_hdrg_supplements.py');m12=load('12_validate_hdrg_components_in_human_sural_nerve.py')
 info=pd.read_csv(m12.NCBI/'Homo_sapiens.gene_info.gz',sep='\t',dtype=str,keep_default_na=False)
 exact=defaultdict(set);fold=defaultdict(set);alias=defaultdict(set)
 for r in info.itertuples():
  exact[r.Symbol].add(r.GeneID);fold[r.Symbol.upper()].add(r.GeneID)
  for a in r.Synonyms.split('|'):
   if a!='-':alias[a.upper()].add(r.GeneID)
 def resolve(s):
  for name,lookup,key in [('official_exact',exact,s),('official_casefold',fold,s.upper()),('unique_synonym',alias,s.upper())]:
   if key in lookup:return (next(iter(lookup[key])),name) if len(lookup[key])==1 else ('','ambiguous_'+name)
  return '','unmapped'
 records=[]
 for sheet,cfg in m03.TRANSCRIPT_SHEETS.items():
  x=pd.read_excel(m03.DATA_FILE,sheet_name=sheet,skiprows=1)
  for i,row in x.iterrows():
   raw=row.gene;s,datefix=m03.repair_gene(raw);s='' if s is None else s;gid,method=resolve(s)
   fc=float(row.avg_log2FC);q=float(row.p_val_adj)
   records.append(dict(source_sheet=sheet,excel_row=i+3,contrast=cfg['contrast_id'],orientation=cfg['effect_orientation'],raw_gene=str(raw),symbol=s,date_repaired=datefix,gene_id=gid,mapping=method,effect=fc,q=q,direction='up' if fc>0 else 'down' if fc<0 else 'zero',eligible=np.isfinite(fc) and np.isfinite(q) and q<.05 and abs(fc)>.585))
 raw=pd.DataFrame(records);put(raw,'inputs/SOURCE_ALL_ROWS.tsv')
 canon=[]
 for (contrast,g),f in raw[raw.gene_id.ne('')].groupby(['contrast','gene_id'],sort=True):
  d=set(f.direction);status='eligible' if len(d)==1 and 'zero' not in d and f.eligible.all() else 'within_source_direction_conflict' if len(d)>1 or 'zero' in d else 'threshold_or_finite_disagreement'
  canon.append(dict(contrast=contrast,gene_id=g,symbol=info.set_index('GeneID').at[g,'Symbol'],direction=next(iter(d)) if len(d)==1 else 'conflict',status=status,n_rows=len(f),source_symbols='|'.join(sorted(set(f.symbol))),source_rows='|'.join(map(str,f.excel_row)),median_effect=f.effect.median()))
 c=pd.DataFrame(canon);put(c,'inputs/CANONICAL_SOURCE_COMPARISONS.tsv')
 lookup={k:f.set_index('gene_id') for k,f in c.groupby('contrast')}
 parentnames={'late_allcell_DPN_vs_diabetes':'original_late_allcell','late_neuron_DPN_vs_diabetes':'original_late_neuron','severity_neuron_modhigh_vs_low_nageotte':'original_severity','early_allcell_diabetes_vs_control':'original_early_allcell','xenium_DPN_vs_control':'original_xenium'}
 members=[];parts=[]
 def member(mod,g,r):members.append(dict(module_id=mod,gene_id=g,symbol=r.symbol,direction=r.direction,source_contrast=r.contrast))
 for k,f in lookup.items():
  for g,r in f[f.status=='eligible'].iterrows():member(parentnames[k],g,r)
 families=[('late_allcell','late_allcell_DPN_vs_diabetes','late_neuron_DPN_vs_diabetes','late_shared_concordant_neuronal_core','late_allcell_residual'),('late_neuron','late_neuron_DPN_vs_diabetes','late_allcell_DPN_vs_diabetes','late_shared_concordant_neuronal_core','late_neuron_residual'),('severity','severity_neuron_modhigh_vs_low_nageotte','late_neuron_DPN_vs_diabetes','severity_neuron_shared_concordant_core','severity_neuron_residual')]
 for fam,pk,ok,core,res in families:
  p=lookup[pk];o=lookup[ok]
  for g,r in p[p.status=='eligible'].iterrows():
   od=o.loc[g].direction if g in o.index else 'absent';os=o.loc[g].status if g in o.index else 'absent'
   bin='residual' if g not in o.index else 'excluded_partner_unresolved' if os!='eligible' else 'core' if od==r.direction else 'excluded_cross_source_opposed'
   parts.append(dict(family=fam,parent=parentnames[pk],gene_id=g,symbol=r.symbol,direction=r.direction,bin=bin,partner_direction=od,partner_status=os))
   if bin=='core' and fam!='late_allcell':member(core,g,r)
   if bin=='residual':member(res,g,r)
 members=pd.DataFrame(members).sort_values(['module_id','gene_id']);parts=pd.DataFrame(parts)
 assert not members.duplicated(['module_id','gene_id']).any()
 for fam,f in parts.groupby('family'):
  assert f.gene_id.is_unique
  p=members[members.module_id==f.parent.iloc[0]].set_index('gene_id').direction
  assert f.set_index('gene_id').direction.sort_index().equals(p.sort_index())
 assert set(parts.query("family=='late_allcell' and bin=='core'").gene_id)==set(parts.query("family=='late_neuron' and bin=='core'").gene_id)
 put(members,'inputs/REPAIRED_MEMBERS.tsv');put(parts,'inputs/REPAIRED_PARTITIONS.tsv')
 put(parts.groupby(['family','bin']).size().rename('genes').reset_index(),'results/PARTITION_COUNTS.tsv')
 put(c.groupby(['contrast','status']).size().rename('genes').reset_index(),'results/SOURCE_STATUS_COUNTS.tsv')
 old=pd.read_csv(B4/'inputs/SOURCE_CANONICAL_MAPPING.tsv',sep='\t',dtype={'human_gene_id':str})
 old=old.rename(columns={'human_gene_id':'gene_id'});comp=old[['module_id','gene_id','gene','direction']].merge(members,on=['module_id','gene_id'],how='outer',suffixes=('_old','_new'),indicator=True)
 put(comp,'results/OLD_NEW_MEMBER_ROW_COMPARISON.tsv')
 files=[m03.DATA_FILE,m12.NCBI/'Homo_sapiens.gene_info.gz',m12.SIGNATURES]
 put(pd.DataFrame([dict(file=f.relative_to(P).as_posix(),sha256=sha(f)) for f in files]),'inputs/SOURCE_FILES.tsv')
 dump(dict(python=platform.python_version(),free_disk_bytes=shutil.disk_usage(B).free,source_rows=len(raw),canonical_status=c.status.value_counts().to_dict(),members=members.groupby('module_id').size().to_dict(),gate='PASS',note='Source repair only. Corrected target associations not run yet.'),'reports/SOURCE_REPAIR_AUDIT.json')
 print((B/'reports/SOURCE_REPAIR_AUDIT.json').read_text(),flush=True)
if __name__=='__main__':run(main)
