from pathlib import Path
import sys,os,json,hashlib,datetime,traceback,importlib.util
sys.dont_write_bytecode=True
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
B=Path(__file__).resolve().parents[1];P=B.parents[1];B4=P/'analysis_v3/batch04'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for v in iter(lambda:f.read(1048576),b''):h.update(v)
 return h.hexdigest()
def seed(*keys):return int.from_bytes(hashlib.sha256(json.dumps([20260910,*keys],separators=(',',':')).encode()).digest()[:16],'big')
def dump(x,p):
 p=B/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
def put(x,p):
 p=B/p;p.parent.mkdir(parents=True,exist_ok=True);x.to_csv(p,sep='\t',index=False)
def load(n):
 p=P/'analysis/scripts'/n;s=importlib.util.spec_from_file_location(n[:-3],p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def run(f):
 start=now();code=0
 try:f()
 except BaseException:code=1;traceback.print_exc();raise
 finally:dump(dict(start_utc=start,end_utc=now(),exit_code=code,command=['python','scripts/'+Path(sys.argv[0]).name,*sys.argv[1:]],code_sha256=sha(sys.argv[0])),f'logs/{Path(sys.argv[0]).stem}_{start.replace(":","-")}.json')
