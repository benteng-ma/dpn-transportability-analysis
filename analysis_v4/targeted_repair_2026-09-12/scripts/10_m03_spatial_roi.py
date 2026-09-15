from __future__ import annotations

import os
from pathlib import Path
from collections import defaultdict
import csv,gzip,hashlib,itertools,json,re
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats
import matplotlib.pyplot as plt

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
ROOT=P/"analysis_v4/targeted_repair_2026-09-12"; OLD=P/"analysis_v4/all_extensions_2026-09-11"
INP=ROOT/"07_resource_recovery/M03_selected_inputs"; OUT=ROOT/"07_resource_recovery/M03_results"; OUT.mkdir(parents=True,exist_ok=True)
MEM=pd.read_csv(P/"analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv",sep="\t",dtype={"gene_id":str})
INFO=P/"data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz"
LABELS={"original_early_allcell":"P1","original_late_allcell":"P2","original_late_neuron":"P3","original_severity":"P4","original_xenium":"P5","late_shared_concordant_neuronal_core":"P6","late_neuron_residual":"P7","late_allcell_residual":"P8","severity_neuron_shared_concordant_core":"P9","severity_neuron_residual":"P10"}
SEED=20260912

def resolver():
    d=pd.read_csv(INFO,sep="\t",dtype=str,keep_default_na=False)[["GeneID","Symbol","Synonyms"]]
    a,b,c=defaultdict(set),defaultdict(set),defaultdict(set)
    for r in d.itertuples(index=False):
        a[r.Symbol].add(r.GeneID); b[r.Symbol.upper()].add(r.GeneID)
        if r.Synonyms and r.Synonyms!="-":
            for x in r.Synonyms.split("|"): c[x.strip().upper()].add(r.GeneID)
    exact={k:next(iter(v)) for k,v in a.items() if len(v)==1}; fold={k:next(iter(v)) for k,v in b.items() if len(v)==1}
    syn={k:next(iter(v)) for k,v in c.items() if len(v)==1 and k not in fold}
    return exact,fold,syn
R=resolver()
def resolve(x):
    x=str(x).strip()
    if x in R[0]: return R[0][x],"official_exact"
    if x.upper() in R[1]: return R[1][x.upper()],"official_casefold"
    if x.upper() in R[2]: return R[2][x.upper()],"unique_synonym"
    return None,"unresolved_or_ambiguous"
def dec(x): return x.decode() if isinstance(x,(bytes,np.bytes_)) else str(x)
def norm_label(x): return re.sub(r"\s+"," ",str(x).strip().lower())
def roi(x):
    x=norm_label(x)
    if "also touching" in x:return "mixed_excluded"
    if x.startswith("nageotte nodule"):return "Nageotte"
    if re.search(r"neurons? near nageotte",x):return "Neuronal"
    return "other_excluded"
def signflip(x):
    x=np.asarray(x,float); obs=x.mean(); null=[]
    for s in itertools.product([-1,1],repeat=len(x)): null.append(np.mean(x*np.asarray(s)))
    return obs,float(np.mean(np.abs(null)>=abs(obs)-1e-15)),len(null)
def boot(x,key):
    x=np.asarray(x,float); seed=int(hashlib.sha256((key+str(SEED)).encode()).hexdigest()[:8],16); r=np.random.default_rng(seed)
    z=np.asarray([r.choice(x,len(x),replace=True).mean() for _ in range(2000)])
    return np.quantile(z,[.025,.975])
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.empty(len(p)); z=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q[o]=np.minimum(z,1);return q

meta=pd.read_csv(OLD/"00_admin/PUBLIC_SAMPLE_METADATA.tsv",sep="\t",dtype=str)
meta=meta[meta.resource.eq("GSE295206")].set_index("geo_accession")
cross=pd.read_csv(OLD/"04_M03_spatial/results/M03_DONOR_DRG_SECTION_MAP.tsv",sep="\t",dtype=str).set_index("Sample ID")
spot_rows=[]; cov_rows=[]; sample_rows=[]; label_rows=[]
for hp in sorted(INP.glob("*filtered_feature_bc_matrix.h5")):
    gsm=hp.name.split("_")[0]; title=meta.at[gsm,"title"]; cm=cross.loc[title]
    sf=next(x for x in INP.glob(gsm+"*.csv.gz") if "tissue_positions" not in x.name)
    cf=next(INP.glob(gsm+"*tissue_positions*.csv.gz"))
    selected=[]
    with gzip.open(sf,"rt",newline="") as f:
        for r in csv.reader(f):
            if len(r)>=2: selected.append((r[0],r[1],roi(r[1])))
    label_rows += [{"geo_accession":gsm,"raw_label":x[1],"normalized_label":norm_label(x[1]),"roi":x[2]} for x in selected]
    with h5py.File(hp,"r") as h:
        m=h["matrix"]; shape=tuple(m["shape"][:]); bars=[dec(x) for x in m["barcodes"][:]]; names=[dec(x) for x in m["features/name"][:]]
        mat=sparse.csc_matrix((m["data"][:],m["indices"][:],m["indptr"][:]),shape=shape)
    bidx={b:i for i,b in enumerate(bars)}; selected=[x for x in selected if x[0] in bidx and x[2] in ("Nageotte","Neuronal")]
    gids=[]; methods=[]
    for x in names:
        g,me=resolve(x); gids.append(g); methods.append(me)
    ok=[i for i,g in enumerate(gids) if g is not None]
    ug=sorted(set(gids[i] for i in ok),key=lambda z:int(z)); gi={g:i for i,g in enumerate(ug)}
    agg=sparse.csr_matrix((np.ones(len(ok)),([gi[gids[i]] for i in ok],ok)),shape=(len(ug),len(names)))@mat
    nz=np.asarray(agg.sum(axis=1)).ravel()>0; agg=agg[nz]; ug=np.asarray(ug)[nz]
    cols=[bidx[x[0]] for x in selected]; dense=agg[:,cols].toarray().astype(float)
    ranks=stats.rankdata(dense,axis=0,method="average")/dense.shape[0]
    gene_index={g:i for i,g in enumerate(ug)}
    for module,fr in MEM.groupby("module_id",sort=True):
        up=[gene_index[g] for g in fr.loc[fr.direction.eq("up"),"gene_id"] if g in gene_index]
        dn=[gene_index[g] for g in fr.loc[fr.direction.eq("down"),"gene_id"] if g in gene_index]
        eligible=len(up)>=10 and len(dn)>=10
        cov_rows.append({"geo_accession":gsm,"sample_id":title,"module_id":module,"program":LABELS[module],"measured_up":len(up),"measured_down":len(dn),"eligible":eligible})
        if not eligible: continue
        U=ranks[up].mean(0);D=ranks[dn].mean(0)
        for k,x in enumerate(selected): spot_rows.append({"geo_accession":gsm,"sample_id":title,"donor_id":cm["Donor ID"],"DRG_ID":cm["DRG ID"],"barcode":x[0],"roi":x[2],"module_id":module,"program":LABELS[module],"U":U[k],"D":D[k],"S":U[k]-D[k]})
    with gzip.open(cf,"rt",newline="") as f: coord=list(csv.reader(f))
    selected_b={x[0] for x in selected}
    sample_rows.append({"geo_accession":gsm,"sample_id":title,"donor_id":cm["Donor ID"],"DRG_ID":cm["DRG ID"],"h5_spots":len(bars),"coordinate_rows":len(coord),"selected_valid_spots":len(selected),"selected_all_in_h5":len(selected_b)==len(selected),"mapped_feature_rows":len(ok),"eligible_background_geneids":len(ug)})

spots=pd.DataFrame(spot_rows); cov=pd.DataFrame(cov_rows); samples=pd.DataFrame(sample_rows); labels=pd.DataFrame(label_rows)
spots.to_csv(OUT/"M03_SPOT_SCORES.tsv.gz",sep="\t",index=False,compression="gzip");cov.to_csv(OUT/"M03_PROGRAM_COVERAGE.tsv",sep="\t",index=False);samples.to_csv(OUT/"M03_SECTION_QC.tsv",sep="\t",index=False)
labels.groupby(["raw_label","normalized_label","roi"]).size().reset_index(name="spots").to_csv(OUT/"M03_LABEL_AUDIT.tsv",sep="\t",index=False)
sec=spots.groupby(["sample_id","donor_id","DRG_ID","module_id","program","roi"])[["S","U","D"]].agg(["mean","count"]).reset_index()
sec.columns=["_".join(x).rstrip("_") if isinstance(x,tuple) else x for x in sec.columns]
wide=sec.pivot(index=["sample_id","donor_id","DRG_ID","module_id","program"],columns="roi",values=["S_mean","S_count","U_mean","D_mean"]).reset_index()
wide.columns=["_".join(x).rstrip("_") if isinstance(x,tuple) else x for x in wide.columns]
wide["S_delta"]=wide["S_mean_Nageotte"]-wide["S_mean_Neuronal"];wide["U_delta"]=wide["U_mean_Nageotte"]-wide["U_mean_Neuronal"];wide["D_delta"]=wide["D_mean_Nageotte"]-wide["D_mean_Neuronal"]
wide["section_evaluable"]=(wide.S_count_Nageotte>=5)&(wide.S_count_Neuronal>=5)
wide.to_csv(OUT/"M03_SECTION_ROI_EFFECTS.tsv",sep="\t",index=False)
ev=wide[wide.section_evaluable].copy()
drg=ev.groupby(["donor_id","DRG_ID","module_id","program"])[["S_delta","U_delta","D_delta"]].mean().reset_index();drg.to_csv(OUT/"M03_DRG_EFFECTS.tsv",sep="\t",index=False)
don=drg.groupby(["donor_id","module_id","program"])[["S_delta","U_delta","D_delta"]].mean().reset_index();don.to_csv(OUT/"M03_DONOR_EFFECTS.tsv",sep="\t",index=False)
tests=[]
for (m,p),fr in don.groupby(["module_id","program"]):
    if len(fr)<3: tests.append({"module_id":m,"program":p,"status":"NOT_EVALUABLE","n_donors":len(fr)});continue
    eff,pv,nperm=signflip(fr.S_delta);lo,hi=boot(fr.S_delta,m)
    tests.append({"module_id":m,"program":p,"status":"COMPLETED_EXPLORATORY","n_donors":len(fr),"effect":eff,"CI_low":lo,"CI_high":hi,"P":pv,"permutations":nperm})
t=pd.DataFrame(tests); ok=t.status.eq("COMPLETED_EXPLORATORY");t.loc[ok,"q"]=bh(t.loc[ok,"P"]);t.to_csv(OUT/"M03_PROGRAM_TESTS.tsv",sep="\t",index=False)
fig,ax=plt.subplots(figsize=(7,4.8)); tt=t[ok].sort_values("effect");y=np.arange(len(tt));ax.errorbar(tt.effect,y,xerr=[tt.effect-tt.CI_low,tt.CI_high-tt.effect],fmt='o',color='#2f4858',ecolor='#7a8a93',capsize=2);ax.axvline(0,color='black',lw=.8);ax.set_yticks(y,tt.program);ax.set_xlabel('Donor weighted Nageotte minus neuronal score');fig.tight_layout();fig.savefig(ROOT/"07_resource_recovery/M03_program_effects.png",dpi=300);plt.close(fig)
summary={"sections":int(samples.sample_id.nunique()),"DRGs":int(samples.DRG_ID.nunique()),"donors":int(samples.donor_id.nunique()),"selected_spots":int(spots[['sample_id','barcode']].drop_duplicates().shape[0]),"evaluable_programs":int(ok.sum()),"q_lt_005":int((t.q<.05).sum()),"gradient":"NOT_EVALUABLE_NO_VERIFIED_SPOT_DIAMETER_CALIBRATION","diagnosis":"diabetes_history_not_explicit_DPN","independence":"UNKNOWN"}
(ROOT/"logs/10_m03_spatial_roi.json").write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))
