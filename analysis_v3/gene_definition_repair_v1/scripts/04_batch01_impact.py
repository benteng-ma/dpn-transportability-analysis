from common import *
import pandas as pd,numpy as np
s=importlib.util.spec_from_file_location('adapter_modules',B/'scripts/03_historical_neural_cross.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)

def main():
 p=P/'analysis_v3/scripts/04_run_batch01_frozen_score_associations.py';spec=importlib.util.spec_from_file_location('batch01_original_functions',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 cfg=json.loads(m.CONFIG_PATH.read_text());scores=pd.read_csv(B/'results/DONOR_SCORES.tsv.gz',sep='\t');allr=[]
 for version in ['repaired_inherited','repaired_statfree']:
  for acc in ['GSE24290','GSE148059']:
   eid=acc+('_C01' if acc=='GSE24290' else '_C02')
   ss=scores[(scores.endpoint==eid)&(scores.version==version)].copy();ss['cohort']=acc;ss['mapping_set']='primary_unambiguous';ss['module_display']=ss.module_id.map(m.DISPLAY)
   meta=pd.read_csv(B/'inputs/matrices'/f'{eid}_metadata.tsv',sep='\t');meta['progression_group']=meta.group;meta['dmfd_class']=meta.group;meta['pathology_order']=meta.group.map({'Regenerator':0,'Intermediate':1,'Degenerator':2})
   out=m.run_group_analyses(ss[['cohort','mapping_set','module_id','module_display','sample_id','patient_id','score']],meta[['sample_id','patient_id','progression_group','dmfd_class','pathology_order']],cfg);out['version']=version;allr.append(out)
   put(out,f'results/historical/{acc}_{version}_Batch01.tsv');print('Batch01',acc,version,flush=True)
 old=pd.read_csv(P/'analysis_v3/results/PATHOLOGY_ASSOCIATION_RESULTS.tsv',sep='\t');old=old[old.mapping_set=='primary_unambiguous'];new=pd.concat(allr)
 comp=new.merge(old,on=['cohort','module_id','mapping_set','test_family'],suffixes=('_new','_old'),validate='many_to_one');put(comp,'results/BATCH01_PATHOLOGY_OLD_NEW.tsv')
 # Actual five paired keys are inherited from the archived input and matched S/T identities, not guessed across studies.
 mods=m.resolve_modules(a.module_frames(),m.M12.build_symbol_lookup(m.M12.load_gene_info(m.NCBI)))
 x,meta=m.jci_paired_matrix(m.M12.build_symbol_lookup(m.M12.load_gene_info(m.NCBI)))
 ss,cv,mp=m.score_matrix('JCI184075_paired',x,meta,mods,m.PAIRED_MODULES,None);out,diff=m.run_paired(ss,cfg)
 put(ss,'results/PAIRED_SCORES.tsv');put(cv,'results/PAIRED_COVERAGE.tsv');put(mp,'results/PAIRED_MAPPING.tsv');put(diff,'results/PAIRED_DIFFERENCES.tsv');put(out,'results/PAIRED_RESULTS.tsv')
 old=pd.read_csv(P/'analysis_v3/results/PAIRED_NERVE_RESULTS.tsv',sep='\t');put(out.merge(old,on=['cohort','module_id','mapping_set','test_family'],suffixes=('_new','_old'),validate='one_to_one'),'results/PAIRED_OLD_NEW.tsv')
 # Official-only was a source spelling sensitivity. Canonical output names are official by construction; do not pretend this is the old sensitivity.
 dump(dict(note='Primary pathology and paired families rerun with original 100000 permutations/5000 bootstrap. The old official-spelling sensitivity is not equivalent after canonicalization and is not reclassified as validation.',source_script_sha256=sha(p),config_sha256=sha(m.CONFIG_PATH)),'reports/BATCH01_ADAPTER_AUDIT.json')
if __name__=='__main__':run(main)
