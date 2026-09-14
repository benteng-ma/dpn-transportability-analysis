from common import *
import pandas as pd

def main():
 for name,file,table,keys in [
  ('rat_parents','04_project_hdrg_stages_to_gse176017.py','GSE176017_human_stage_projection_tests_2026-08-27.tsv',['contrast_id','earlier_group','later_group']),
  ('ocular_parents','07_project_hdrg_stages_to_diabetic_tg_cornea.py','ocular_hDRG_stage_projection_tests_2026-08-27.tsv',['dataset','contrast_id'])]:
  m=load(file);out=B/'results/historical'/name;out.mkdir(exist_ok=True);m.TABLES=out;m.FIGURES=B/'figures'/name;m.FIGURES.mkdir(exist_ok=True);m.SIGNATURES=B/'inputs/REPAIRED_SIGNATURE_ADAPTER.tsv'
  if name=='ocular_parents':m.plot_matrix=lambda *a:None
  if not (out/table).exists():m.main()
  nn=pd.read_csv(out/table,sep='\t');oo=pd.read_csv(P/'results/tables'/table,sep='\t')
  put(nn.merge(oo,on=keys,how='outer',suffixes=('_new','_old'),validate='one_to_one'),'results/'+name.upper()+'_OLD_NEW.tsv')
  print(name,'completed',flush=True)
if __name__=='__main__':run(main)
