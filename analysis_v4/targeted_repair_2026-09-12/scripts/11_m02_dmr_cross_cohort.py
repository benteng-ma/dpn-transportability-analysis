from pathlib import Path
import json,re
import numpy as np,pandas as pd

R = Path(__file__).resolve().parents[1]
O=R/"07_resource_recovery/M02_results"
def load(c):
 d=pd.read_csv(O/f"M02_{c}_DMR_ALL.tsv",sep="\t"); z=d.coord.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$");d[["chr","start","end"]]=z;d[["start","end"]]=d[["start","end"]].astype(int);d["dmr_id"]=[f"{c}_DMR_{i+1:04d}" for i in range(len(d))];return d
g,e=load("PROPGER"),load("PROPENG")
pairs=[]
for ch,a in g.groupby("chr"):
 b=e[e.chr.eq(ch)]
 for x in a.itertuples(index=False):
  q=b[(b.start<=x.end)&(b.end>=x.start)]
  for y in q.itertuples(index=False):
   pairs.append({"PROPGER_dmr":x.dmr_id,"PROPENG_dmr":y.dmr_id,"chr":ch,"overlap_start":max(x.start,y.start),"overlap_end":min(x.end,y.end),"PROPGER_meandiff":x.meandiff,"PROPENG_meandiff":y.meandiff,"direction_concordant":np.sign(x.meandiff)==np.sign(y.meandiff)})
pp=pd.DataFrame(pairs);pp.to_csv(O/"M02_DMR_CROSS_COHORT_OVERLAPS.tsv",sep="\t",index=False)
s={"PROPGER_DMRs":len(g),"PROPENG_DMRs":len(e),"overlap_pairs":len(pp),"PROPGER_DMRs_with_overlap":int(pp.PROPGER_dmr.nunique()) if len(pp) else 0,"PROPENG_DMRs_with_overlap":int(pp.PROPENG_dmr.nunique()) if len(pp) else 0,"direction_concordant_pairs":int(pp.direction_concordant.sum()) if len(pp) else 0}
(O/"M02_DMR_CROSS_COHORT_SUMMARY.json").write_text(json.dumps(s,indent=2),encoding="utf-8");print(json.dumps(s,indent=2))
