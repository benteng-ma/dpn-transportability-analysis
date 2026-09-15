from common import *
import pandas as pd,numpy as np
from scipy.stats import spearmanr
s=importlib.util.spec_from_file_location('mods_adapter',B/'scripts/03_historical_neural_cross.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
def main():
 m12=load('12_validate_hdrg_components_in_human_sural_nerve.py');m24=load('24_expand_exploratory_analysis.py');look=m12.build_symbol_lookup(m12.load_gene_info(m12.NCBI/'Homo_sapiens.gene_info.gz'));mods=a.module_frames()
 proteins=pd.read_csv(P/'results/tables/hDRG_proteomics_contrasts_audited_2026-08-27.tsv.gz',sep='\t');pm=m12.add_resolution(proteins,'gene',look);put(pm,'inputs/PROTEIN_CANONICAL_ROW_MAPPING.tsv.gz');summary=[];matched=[]
 plan=[*( (x,'protein_DPN_vs_diabetes','late') for x in m24.LATE_MODULES),*( (x,'protein_severity_modhigh_vs_low_nageotte','severity') for x in m24.SEVERITY_MODULES)]
 for i,(mod,contrast,family) in enumerate(plan):
  target=pm[(pm.contrast_id==contrast)&pm.human_gene_id.notna()].copy();target=target[np.isfinite(target.t)&np.isfinite(target.AveExpr)]
  target=target.groupby('human_gene_id',sort=True)[['t','AveExpr','logFC']].median().reset_index();target['expression_decile']=pd.qcut(target.AveExpr.rank(method='first'),10,labels=False)
  src=mods[mod].rename(columns={'gene_id':'human_gene_id'});mapped=src.merge(target,on='human_gene_id');mapped['module_id']=mod;matched.append(mapped)
  nu=mapped.direction.eq('up').sum();nd=mapped.direction.eq('down').sum();row=dict(module_id=mod,family=family,contrast=contrast,n_up=nu,n_down=nd,status='NOT_EVALUABLE',p=np.nan,statistic=np.nan)
  if nu and nd:
   val,mu,sd,p=m24.expression_matched_protein_null(mapped,target,np.random.default_rng(m24.SEED+i));row.update(status='COMPLETED',p=p,statistic=val,null_mean=mu,null_sd=sd,source_target_spearman=spearmanr(mapped.source_log2fc,mapped.logFC).statistic)
  summary.append(row)
 out=pd.DataFrame(summary);out['q']=out.groupby('family',group_keys=False).p.apply(m12.bh_adjust);old=pd.read_csv(P/'results/expansion_v2/tables/Supplementary_Table_S13_source_transcript_protein_concordance.tsv',sep='\t');put(out.merge(old,on=['module_id','family'],suffixes=('_new','_old'),validate='one_to_one'),'results/PROTEIN_OLD_NEW.tsv');put(pd.concat(matched),'results/PROTEIN_MATCHED_ALL.tsv')
 # Keep all original eigengenes/modules; recompute only membership-dependent descriptive relationships.
 b2=P/'analysis_v3/batch02';mem=pd.read_csv(b2/'results/network/MODULE_MEMBERS.tsv',sep='\t');me=pd.read_csv(b2/'results/network/MODULE_EIGENGENES.tsv',sep='\t').set_index('sample_id');expr=pd.read_csv(b2/'inputs/GSE148059_analysis_expression.tsv.gz',sep='\t',index_col=0);expr.index=expr.index.str.upper();z=(expr.T-expr.T.mean())/expr.T.std(ddof=1)
 scores=pd.read_csv(B/'results/DONOR_SCORES.tsv.gz',sep='\t');sc=scores[(scores.endpoint=='GSE148059_C02')&(scores.version=='repaired_statfree')].pivot(index='sample_id',columns='module_id',values='score').loc[me.index];rel=[]
 for module,f in mem[mem.module!='grey'].groupby('module'):
  genes=set(f.gene.str.upper())
  for mod in sc:
   pg=set(mods[mod].gene.str.upper());shared=genes&pg;remain=genes-pg;ok=len(remain)>=30 and len(remain)/len(genes)>=.5 and sc[mod].notna().all();rel.append(dict(module=module,program=mod,shared_n=len(shared),remaining_n=len(remain),eigengene_score_spearman=spearmanr(me['ME'+module],sc[mod]).statistic if sc[mod].notna().all() else np.nan,nonshared_mean_score_spearman=spearmanr(z[sorted(remain)].mean(axis=1),sc[mod]).statistic if ok else np.nan,status='COMPLETED' if ok else 'NOT_EVALUABLE'))
 old=pd.read_csv(b2/'results/network/MODULE_PROGRAM_RELATIONSHIPS.tsv',sep='\t');put(pd.DataFrame(rel).merge(old,on=['module','program'],suffixes=('_new','_old'),validate='one_to_one'),'results/MODULE_PROGRAM_OLD_NEW.tsv')
 # TPM: same deposited expression source and original sample labels, no new patients.
 eid='original_JCI_sural_DPN_control';meta=pd.read_csv(B/'inputs/matrices'/f'{eid}_metadata.tsv',sep='\t');groups=dict(zip(meta.sample_id,meta.group));u,x,mm,au=m12.load_target('DPN_vs_control',m12.FILES['dpn_control'],'deseq2_results',groups,1,look);mapped,ma,cv=m12.map_modules(mods,u,'DPN_vs_control',m12.DPN_MODULES,look)
 x,au=m12.load_tpm_sensitivity(m12.FILES['quantile_tpm'],groups,look);sc,tests=m12.score_modules('DPN_vs_control',x,groups,'DPN','Control',mapped,m12.DPN_MODULES,'quantile_normalized_TPM_sensitivity');put(sc,'results/JCI_TPM_SCORES.tsv');old=pd.read_csv(P/'results/tables/JCI184075_hDRG_component_sample_tests_2026-08-27.tsv',sep='\t');old=old[old.expression_source=='quantile_normalized_TPM_sensitivity'];put(tests.merge(old,on=['target_id','module_id','expression_source'],suffixes=('_new','_old'),validate='one_to_one'),'results/JCI_TPM_OLD_NEW.tsv')
 # Regulator sensitivity exclusion is based on canonical IDs, not target-association screening.
 prior=pd.read_csv(b2/'results/regulators/REGULON_SNAPSHOT.tsv',sep='\t');vocab=pd.DataFrame({'gene':sorted(set(prior.target.str.upper()))});vocab=m12.add_resolution(vocab,'gene',look);vocab['exclude']=vocab.human_gene_id.isin(set(pd.concat(mods.values()).gene_id));put(vocab,'inputs/REGULATOR_EXCLUSION_VOCABULARY.tsv')
 print('Protein, network relationships and TPM impact completed',flush=True)
if __name__=='__main__':run(main)
