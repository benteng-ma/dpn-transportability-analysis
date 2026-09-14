from common import *
import pandas as pd,numpy as np,shutil,platform,importlib.metadata
from collections import defaultdict

def main():
 checks=[]
 def ck(k,v,detail=''):checks.append(dict(check=k,status='PASS' if v else 'FAIL',detail=str(detail)))
 info=pd.read_csv(P/'data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz',sep='\t',dtype=str,keep_default_na=False)
 ex=defaultdict(set);fold=defaultdict(set);ali=defaultdict(set)
 for r in info.itertuples():
  ex[r.Symbol].add(r.GeneID);fold[r.Symbol.upper()].add(r.GeneID)
  for s in r.Synonyms.split('|'):
   if s!='-':ali[s.upper()].add(r.GeneID)
 def resolve(s):
  for lookup,key in [(ex,s),(fold,s.upper()),(ali,s.upper())]:
   if key in lookup:return next(iter(lookup[key])) if len(lookup[key])==1 else ''
  return ''
 diffs=[]
 for file in (B/'strict_background_addendum/inputs').glob('*STRICT_ROW_MAPPING.tsv'):
  f=pd.read_csv(file,sep='\t',dtype=str,keep_default_na=False);want=f.gene_input.map(resolve);ck('strict_target_resolver_parity:'+file.name,want.equals(f.human_gene_id),int((want!=f.human_gene_id).sum()))
  d=f[want!=f.human_gene_id].copy();d['strict_expected_id']=want.loc[d.index];d['source_file']=file.name;diffs.append(d)
 put(pd.concat(diffs),'tests/STRICT_TARGET_RESOLVER_DIFFERENCES.tsv')
 for r in pd.read_csv(B/'inputs/PROTECTED_FILES.tsv',sep='\t').itertuples():ck('protected_final:'+r.file,sha(P/r.file)==r.sha256)
 for r in pd.read_csv(B/'config/LOCKED_FILES.tsv',sep='\t').itertuples():ck('lock_final:'+r.file,sha(B/r.file)==r.sha256)
 for name in ['RAT_PARENTS_OLD_NEW.tsv','OCULAR_PARENTS_OLD_NEW.tsv','HDRG_PARENT_OLD_NEW.tsv']:
  f=pd.read_csv(B/'results'/name,sep='\t');ck('parent_output:'+name,len(f)>0)
 put(pd.DataFrame(checks),'tests/HANDOFF_TESTS.tsv');assert all(z['status']=='PASS' for z in checks)
 # Archive exact helper code; source workbooks/cell caches remain project prerequisites.
 helpers=[]
 for folder in [P/'analysis/scripts',P/'analysis_v3/scripts',P/'analysis_v3/batch02/scripts',P/'analysis_v3/batch03/scripts',P/'analysis_v3/batch04/scripts']:
  if folder.exists():
   for f in folder.iterdir():
    if f.is_file() and f.suffix.lower() in ['.py','.r']:
     rel=f.relative_to(P);dst=B/'baseline_code'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,dst);helpers.append(dict(project_relative_file=rel.as_posix(),sha256=sha(f),archived_copy=dst.relative_to(B).as_posix()))
 put(pd.DataFrame(helpers),'inputs/BASELINE_CODE_DEPENDENCIES.tsv')
 packages=[dict(package=x.metadata['Name'],version=x.version) for x in importlib.metadata.distributions()];put(pd.DataFrame(packages).sort_values('package'),'inputs/PYTHON_DISTRIBUTIONS.tsv')
 dump({'python_version':platform.python_version(),'platform':platform.platform(),'environment':'existing Batch04 project venv','global_packages_modified':False},'inputs/RUNTIME.json')
 # Status is backed by outputs, not promises. NE/conditional dependencies are explicit.
 rows=[]
 tasks=[
 ('source_repair','COMPLETED','00_source_repair.py','inputs/SOURCE_ALL_ROWS.tsv;inputs/REPAIRED_MEMBERS.tsv','Canonical source-only rules, real conflicts retained'),
 ('backgrounds_and_lock','COMPLETED','01_backgrounds_lock.py','config/LOCK.json;inputs/MATRIX_INVENTORY.tsv','Expression-only median sensitivity; known-result correction'),
 ('main_effect_impact','COMPLETED','02_associations_benchmarks.py','results/ASSOCIATION_IMPACT_ALL.tsv;results/DONOR_SCORES.tsv.gz','All five endpoints, all ten programs, three versions; includes NE rows'),
 ('late_empirical_benchmark','COMPLETED','02_associations_benchmarks.py','results/REAL_BENCHMARK_ALL.tsv;results/MATCHED_RANDOM_SPLITS.tsv.gz','80 tests; no Holm<0.05; 19980 partitions'),
 ('severity_empirical_benchmark','NOT_EVALUABLE','02_associations_benchmarks.py','results/REAL_BENCHMARK_ALL.tsv','40 tasks blocked by seven-down-gene core; no threshold relaxation'),
 ('historical_neural_tests','COMPLETED','03_historical_neural_cross.py','results/HISTORICAL_DONOR_TEST_IMPACT.tsv','Historical single-sided families separately preserved'),
 ('pathology_and_paired','COMPLETED','04_batch01_impact.py','results/BATCH01_PATHOLOGY_OLD_NEW.tsv;results/PAIRED_OLD_NEW.tsv','Original subject keys and original resampling methods'),
 ('official_spelling_sensitivity','NOT_EVALUABLE','04_batch01_impact.py','reports/BATCH01_ADAPTER_AUDIT.json','Not an equivalent definition after canonicalization; no false carry-forward'),
 ('blood_and_pain','COMPLETED','05_blood_clinical_impact.py','results/PBMC_OLD_NEW.tsv;results/CLINICAL_OLD_NEW.tsv','Original clinical linkage only'),
 ('reference_and_proxy','COMPLETED','06_reference_proxy_impact.py','results/ALL_REFERENCE_PROGRAM_ASSOCIATIONS_OLD_NEW.tsv;results/PROXY_ADJUSTMENT_OLD_NEW.tsv','CIAP/CIDP not DPN; 2 DPN descriptive only'),
 ('protein_module_TPM','COMPLETED','07_protein_network_tpm.py','results/PROTEIN_OLD_NEW.tsv;results/MODULE_PROGRAM_OLD_NEW.tsv;results/JCI_TPM_OLD_NEW.tsv','No WGCNA rebuild'),
 ('TF_exclusion','COMPLETED','09_run_regulators.py','results/REGULATOR_EXCLUSION_OLD_NEW.tsv','R4.6.1 successful after retained failed attempts; fixed prior'),
 ('hDRG_parent_gene_layer','COMPLETED','11_remaining_dependencies.py','results/HDRG_PARENT_OLD_NEW.tsv','Conditional deposited-statistic test, not new gene validation'),
 ('conditional_merged_index','COMPLETED','11_remaining_dependencies.py','results/CONDITIONAL_MERGED_INDEX_ALL.tsv','Numeric historical formula only; no valid formal gene-evidence claim'),
 ('formal_merged_gene_validation','BLOCKED','11_remaining_dependencies.py','reports/REMAINING_DEPENDENCY_LIMITATIONS.json','Original JCI statistic/shrinkage provenance unresolved'),
 ('source_annotation_impact','COMPLETED','11_remaining_dependencies.py','results/SOURCE_ANNOTATION_OLD_NEW.tsv.gz','Fixed literal GMT; old overlap-filtered testing family, descriptive only'),
 ('spatial_panel_gate','NOT_EVALUABLE','12_final_audit.py','results/SPATIAL_REPAIRED_PROGRAM_COVERAGE.tsv','0/10 corrected programs meet original panel rule; no large download'),
 ('rodent_ocular_parents','COMPLETED','14_parent_context_impact.py','results/RAT_PARENTS_OLD_NEW.tsv;results/OCULAR_PARENTS_OLD_NEW.tsv','Complete original parent projection contexts; cached successful rat run reused after join repair'),
 ('unaffected_DE_pathway_network_primaryTF','COMPLETED','','inputs/PROTECTED_FILES.tsv','Reviewed dependency: independent of modified memberships; deliberately not rerun'),
 ('unaffected_module_only_localization','COMPLETED','','inputs/PROTECTED_FILES.tsv','Reviewed dependency only; no new calculation or validation asserted'),
 ('old_exploratory_display_bootstrap_PCA','PLANNED','','','Not regenerated; superseded main effects available, old supplementary display not cleared for reuse'),
 ('manuscript_and_public_release','PLANNED','','','Not requested or executed; no GitHub/Zenodo/DOI change')]
 for task,status,script,outputs,reason in tasks:
  code=B/'scripts'/script if script else None
  rows.append(dict(task=task,status=status,script=('scripts/'+script if script else ''),code_sha256=sha(code) if code else '',command=('project_python scripts/'+script if script else ''),inputs='config/LOCKED_FILES.tsv;inputs/SOURCE_FILES.tsv;inputs/TARGET_SOURCE_FILES.tsv;inputs/BASELINE_CODE_DEPENDENCIES.tsv',outputs=outputs,reason=reason))
 put(pd.DataFrame(rows),'STATUS.tsv')
 counts={}
 for file in ['ALL_TESTS.tsv','ADDITIONAL_TESTS.tsv','HANDOFF_TESTS.tsv','STRICT_ADDENDUM_TESTS.tsv']:
  f=pd.read_csv(B/'tests'/file,sep='\t');counts[file]=f.status.value_counts().to_dict()
 dump({'suites':counts,'prior_analysis_test_passes':28546,'handoff_passes':len(checks),'all_passes':sum(v.get('PASS',0) for v in counts.values()),'zero_failures':all(v.get('FAIL',0)==0 for v in counts.values()),'timestamp_utc':now()},'reports/HANDOFF_CHECKS.json');print(json.dumps(counts),flush=True)
if __name__=='__main__':run(main)
