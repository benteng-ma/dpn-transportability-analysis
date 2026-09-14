from common import *
import pandas as pd,numpy as np
from scipy import stats

def main():
 m=load('06_validate_hdrg_signatures_in_independent_human_bulk.py');out=B/'results/historical/hDRG_parents';out.mkdir(exist_ok=True)
 m.TABLES=out;m.METADATA=B/'inputs/hDRG_parent_metadata';m.FIGURES=B/'figures/hDRG_parents';m.SIGNATURES=B/'inputs/REPAIRED_SIGNATURE_ADAPTER.tsv';m.plot_results=lambda *a:None
 m.main()
 old=pd.read_csv(P/'results/tables/independent_human_DPN_bulk_signature_tests_2026-08-27.tsv',sep='\t');new=pd.read_csv(out/'independent_human_DPN_bulk_signature_tests_2026-08-27.tsv',sep='\t')
 put(new.merge(old,on='contrast_id',suffixes=('_new','_old'),validate='one_to_one'),'results/HDRG_PARENT_OLD_NEW.tsv')
 # Historical merged index is recalculated solely for traceability, not validated z inference.
 h=pd.read_csv(B/'results/historical/human_components_2.tsv.gz',sep='\t',dtype={'human_gene_id':str})
 s=pd.read_csv(B/'results/historical/DPN_vs_control_gene_mapping.tsv.gz',sep='\t',dtype={'human_gene_id':str})
 late=['late_shared_concordant_neuronal_core','late_neuron_residual','late_allcell_residual']
 h=h[h.module_id.isin(late)];s=s[s.module_id.isin(late)&s.mapped_to_target_universe]
 keys=['module_id','human_gene_id'];z=h[keys+['gene','direction','stat']].merge(s[keys+['target_wald_oriented']],on=keys,validate='one_to_one')
 sign=z.direction.map({'up':1,'down':-1});z['hdrg_oriented_z']=z.stat*sign;z['sural_oriented_z']=z.target_wald_oriented*sign
 z['historical_combined_index']=(z.hdrg_oriented_z+z.sural_oriented_z)/np.sqrt(2)
 z['conditional_p']=stats.norm.sf(z.historical_combined_index);z['conditional_q']=z.groupby('module_id',group_keys=False).conditional_p.apply(m.bh_adjust)
 z['historical_rule_flag']=(z.hdrg_oriented_z>0)&(z.sural_oriented_z>0)&(z.conditional_q<.1)
 z['inference_status']='CONDITIONAL_INDEX_NOT_GENE_VALIDATION_STAT_PROVENANCE_UNRESOLVED'
 put(z,'results/CONDITIONAL_MERGED_INDEX_ALL.tsv')
 # Membership-only annotation impact with the same archived libraries and historical overlap-filtered family.
 sp=importlib.util.spec_from_file_location('corrected_neural',B/'scripts/03_historical_neural_cross.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a);mods=a.module_frames()
 ann=load('14_annotate_hdrg_transport_components.py');ann.TABLES=B/'results/historical/source_annotation';ann.TABLES.mkdir(exist_ok=True)
 ann.M12.build_source_modules=lambda:(mods,None,{'canonical_lock':sha(B/'config/LOCK.json'),'retrospective':True})
 ann.make_figure=lambda *a:None
 ann.main()
 old=pd.read_csv(P/'results/tables/hDRG_component_functional_annotation_all_terms_2026-08-27.tsv.gz',sep='\t');new=pd.read_csv(ann.TABLES/'hDRG_component_functional_annotation_all_terms_2026-08-27.tsv.gz',sep='\t')
 keys=['library','module_id','direction','term'];j=new.merge(old,on=keys,suffixes=('_new','_old'),how='outer',validate='one_to_one',indicator=True)
 j['interpretation']='historical_overlap_filtered_family_descriptive_only_fixed_GMT_symbols'
 put(j,'results/SOURCE_ANNOTATION_OLD_NEW.tsv.gz')
 dump({'historical_annotation_limitation':'Same literal archived GMT gene vocabulary; BH applied only after overlap >=3 as in legacy code. Not a newly calibrated enrichment or clinical test. No new library or favorable-term selection.', 'merged_index_limitation':'Deposited statistic/shrinkage provenance remains unresolved. Numeric index reproduced for impact only, not formal gene-level validation or a new candidate list.'},'reports/REMAINING_DEPENDENCY_LIMITATIONS.json')
 print('remaining fixed-method dependencies completed',flush=True)
if __name__=='__main__':run(main)
