from common import *
import subprocess,shutil
def main():
 # Reuse the installed ABI-compatible runtime and the unchanged DPN package library.
 r=Path(os.environ['LOCALAPPDATA'])/'Programs/R/R-4.6.1/bin/Rscript.exe';assert r.exists()
 if (B/'logs/REGULATOR_STDOUT.txt').exists():shutil.copyfile(B/'logs/REGULATOR_STDOUT.txt',B/'logs'/('REGULATOR_PREVIOUS_ATTEMPT_'+now().replace(':','-')+'.txt'))
 cmd=[str(r),str(B/'scripts/08_regulator_exclusion.R'),str(P/'analysis_v3/batch02'),str(B)]
 with (B/'logs/REGULATOR_STDOUT.txt').open('w',encoding='utf-8') as f:out=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 dump(dict(exit_code=out.returncode,script_sha256=sha(B/'scripts/08_regulator_exclusion.R'),r_version_path='R-4.6.1',finished_utc=now()),'logs/REGULATOR_EXECUTION.json')
 assert out.returncode==0,'See REGULATOR_STDOUT.txt'
 print('Regulator exclusion complete',flush=True)
if __name__=='__main__':run(main)
