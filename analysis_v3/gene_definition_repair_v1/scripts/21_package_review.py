from common import *
import pandas as pd,zipfile,io
from pypdf import PdfReader

def main():
 S=B/'strict_background_addendum';hc=json.loads((B/'reports/HANDOFF_CHECKS.json').read_text());assert hc['zero_failures']
 for row in pd.read_csv(B/'inputs/PROTECTED_FILES.tsv',sep='\t').itertuples():assert sha(P/row.file)==row.sha256
 pdf=B/'reports/REPAIR_REVIEW_REPORT.pdf';reader=PdfReader(pdf);assert len(reader.pages)==21
 text='\n'.join(page.extract_text() for page in reader.pages);assert '37975' in text.replace(',','') and 'Implementation discrepancy' in text and '15:59:02' in text
 assert '\ufffd' not in text
 dump({'pdf_sha256':sha(pdf),'pages':21,'rendered_pages_visually_reviewed':list(range(1,22)),'review_method':'All page PNGs inspected in six contact sheets; figures and tables visible, black text, no running header, three-line tables, no clipped content observed.','text_extraction_checks':'pass; strict addendum and current test count present','render_engine':'bundled Poppler at 95 dpi','review_figures':['Figure_R1_partition','Figure_R2_neural_effects','Figure_R3_ALC_effects','Figure_R4_benchmarks']},'reports/ARTIFACT_QA.json')
 status=pd.read_csv(B/'STATUS.tsv',sep='\t').fillna('')
 def finalpath(value):
  out=[]
  for path in value.split(';'):
   alt=S/'consolidated_results'/Path(path).name
   out.append(alt.relative_to(B).as_posix() if path.startswith('results/') and alt.exists() else path)
  return ';'.join(out)
 status['authoritative_outputs']=status.outputs.map(finalpath)
 extra=pd.DataFrame([{'task':'strict_background_addendum','status':'COMPLETED','script':'scripts/18_strict_background_addendum.py;scripts/19_strict_dependents.py;scripts/20_strict_checks.py','code_sha256':sha(B/'scripts/18_strict_background_addendum.py'),'command':'see logs/18,19,20','inputs':'strict_background_addendum/config/LOCKED_FILES.tsv','outputs':'strict_background_addendum/consolidated_results;reports/STRICT_ADDENDUM_SUMMARY.json','reason':'Four failed initial resolver checks retained. Strict ambiguity-safe background locked and affected dependencies rerun. Original members unchanged.','authoritative_outputs':'strict_background_addendum/consolidated_results'}])
 if not status.task.eq('strict_background_addendum').any():status=pd.concat([status,extra],ignore_index=True)
 put(status,'STATUS.tsv')
 index=[]
 for f in (B/'results').rglob('*'):
  if f.is_file():
   alt=S/'consolidated_results'/f.name;override=alt.exists() and f.parent==B/'results';index.append(dict(initial_file=f.relative_to(B).as_posix(),authoritative_file=(alt if override else f).relative_to(B).as_posix(),strict_replacement=override))
 put(pd.DataFrame(index),'FINAL_VIEW_INDEX.tsv')
 dump({'started_utc':now(),'code_sha256':sha(__file__),'no_external_publication':True},'logs/PACKAGE_EXECUTION_START.json')
 archive=B/'DPN_Gene_Definition_Repair_REVIEW_BUNDLE_2026-09-11_v1.zip';assert not archive.exists(),'Do not overwrite delivered archive'
 def include(f):
  rel=f.relative_to(B);return f.is_file() and 'rendered' not in rel.parts and '__pycache__' not in rel.parts and f.suffix!='.zip' and f.name not in ['MANIFEST_SHA256.tsv','PACKAGE_VERIFICATION.json','ZIP_SHA256.txt'] and not f.name.startswith('21_package_review_') and not (rel.parts[0]=='figures' and f.suffix=='.pdf')
 files=sorted(f for f in B.rglob('*') if include(f));manifest=pd.DataFrame([dict(file=f.relative_to(B).as_posix(),bytes=f.stat().st_size,sha256=sha(f)) for f in files]);put(manifest,'MANIFEST_SHA256.tsv')
 with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files+[B/'MANIFEST_SHA256.tsv']:z.write(f,f.relative_to(B).as_posix())
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None;mf=pd.read_csv(io.BytesIO(z.read('MANIFEST_SHA256.tsv')),sep='\t')
  for row in mf.itertuples():assert hashlib.sha256(z.read(row.file)).hexdigest()==row.sha256,row.file
  assert len(z.namelist())==len(mf)+1
 digest=sha(archive);result={'archive':archive.name,'bytes':archive.stat().st_size,'sha256':digest,'payloads_verified':len(manifest),'CRC':'PASS','all_payload_SHA256':'PASS','protected_files_unchanged':len(pd.read_csv(B/'inputs/PROTECTED_FILES.tsv',sep='\t')),'pdf_pages':21,'current_test_passes':hc['all_passes'],'historical_resolver_failures_retained':4,'completed_utc':now()};dump(result,'reports/PACKAGE_VERIFICATION.json')
 (B/'ZIP_SHA256.txt').write_text(digest+'  '+archive.name+'\n',encoding='utf-8');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':run(main)
