from __future__ import annotations

import os

import csv, gzip, hashlib, itertools, json, re
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

P = Path(os.environ.get("DPN_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
T = P / "analysis_v4/targeted_repair_2026-09-12"
C = P / "analysis_v4/closeout_20260912T062015Z"
INP = T / "07_resource_recovery/M03_selected_inputs"
SRC = T / "07_resource_recovery/M03_results"
OUT = C / "04_M03_spatial"
MEM = pd.read_csv(P / "analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv", sep="\t", dtype={"gene_id": str})
INFO = P / "data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz"
LABELS = {"original_early_allcell":"P1","original_late_allcell":"P2","original_late_neuron":"P3","original_severity":"P4","original_xenium":"P5","late_shared_concordant_neuronal_core":"P6","late_neuron_residual":"P7","late_allcell_residual":"P8","severity_neuron_shared_concordant_core":"P9","severity_neuron_residual":"P10"}

def dec(x): return x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x)
def norm(x): return re.sub(r"\s+", " ", str(x).strip().lower())
def roi(x):
    x = norm(x)
    if "also touching" in x: return "mixed_excluded"
    if x.startswith("nageotte nodule"): return "Nageotte"
    if re.search(r"neurons? near nageotte", x): return "Neuronal"
    return "other_excluded"
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.empty(len(p)); z=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]; q[o]=np.minimum(z,1); return q
def signflip(x):
    x=np.asarray(x,float); obs=x.mean(); null=[np.mean(x*np.asarray(s)) for s in itertools.product([-1,1],repeat=len(x))]
    return obs, float(np.mean(np.abs(null)>=abs(obs)-1e-15)), len(null)

def resolver():
    d=pd.read_csv(INFO,sep="\t",dtype=str,keep_default_na=False)[["GeneID","Symbol","Synonyms"]]
    a,b,c=defaultdict(set),defaultdict(set),defaultdict(set)
    for r in d.itertuples(index=False):
        a[r.Symbol].add(r.GeneID); b[r.Symbol.upper()].add(r.GeneID)
        if r.Synonyms and r.Synonyms!="-":
            for x in r.Synonyms.split("|"): c[x.strip().upper()].add(r.GeneID)
    return ({k:next(iter(v)) for k,v in a.items() if len(v)==1}, {k:next(iter(v)) for k,v in b.items() if len(v)==1}, {k:next(iter(v)) for k,v in c.items() if len(v)==1 and k not in b})
R=resolver()
def resolve(x):
    if x in R[0]: return R[0][x]
    if x.upper() in R[1]: return R[1][x.upper()]
    if x.upper() in R[2]: return R[2][x.upper()]
    return None

OUT.mkdir(parents=True,exist_ok=True)
don = pd.read_csv(SRC/"M03_DONOR_EFFECTS.tsv",sep="\t")
std = don.rename(columns={"S_delta":"effect","U_delta":"up_effect","D_delta":"down_effect"})
std["eligible"] = True; std["reason"] = ""
std = std[["donor_id","program","effect","up_effect","down_effect","eligible","reason"]]
# P9 was structurally ineligible and must remain explicit.
donors=sorted(don.donor_id.unique())
std=pd.concat([std,pd.DataFrame({"donor_id":donors,"program":"P9","effect":np.nan,"up_effect":np.nan,"down_effect":np.nan,"eligible":False,"reason":"RNA full score ineligible: seven source-down members"})],ignore_index=True)
std.to_csv(OUT/"SPATIAL_DONOR_EFFECTS.tsv",sep="\t",index=False,na_rep="")

# Two-arm family: U and D are tested separately across every eligible program.
tests=[]
for program,fr in std[std.eligible].groupby("program",sort=True):
    for arm,col in [("source_up","up_effect"),("source_down","down_effect")]:
        eff,p,n=signflip(fr[col])
        tests.append({"program":program,"arm":arm,"n_donors":len(fr),"effect":eff,"P_two_sided":p,"permutations":n})
arm=pd.DataFrame(tests); arm["q_BH_18_tests"]=bh(arm.P_two_sided); arm.to_csv(OUT/"M03_TWO_ARM_TESTS.tsv",sep="\t",index=False)

# Recompute gene-level ROI rank contributions from the real H5 and author ROI labels.
qc=pd.read_csv(SRC/"M03_SECTION_QC.tsv",sep="\t",dtype=str).set_index("geo_accession")
contrib=[]; coords=[]
for hp in sorted(INP.glob("*filtered_feature_bc_matrix.h5")):
    gsm=hp.name.split("_")[0]
    sf=next(x for x in INP.glob(gsm+"*.csv.gz") if "tissue_positions" not in x.name)
    cf=next(INP.glob(gsm+"*tissue_positions*.csv.gz"))
    labels=[]
    with gzip.open(sf,"rt",newline="") as f:
        for r in csv.reader(f):
            if len(r)>=2: labels.append((r[0],roi(r[1])))
    label_map=dict(labels)
    with gzip.open(cf,"rt",newline="") as f:
        for r in csv.reader(f):
            if len(r)>=6 and r[0] in label_map:
                coords.append({"geo_accession":gsm,"barcode":r[0],"array_row":int(r[2]),"array_col":int(r[3]),"pixel_row":float(r[4]),"pixel_col":float(r[5]),"roi":label_map[r[0]]})
    with h5py.File(hp,"r") as h:
        m=h["matrix"]; shape=tuple(m["shape"][:]); bars=[dec(x) for x in m["barcodes"][:]]; names=[dec(x) for x in m["features/name"][:]]
        mat=sparse.csc_matrix((m["data"][:],m["indices"][:],m["indptr"][:]),shape=shape)
    bidx={b:i for i,b in enumerate(bars)}
    selected=[(b,r) for b,r in labels if b in bidx and r in ("Nageotte","Neuronal")]
    gids=[resolve(x) for x in names]; ok=[i for i,g in enumerate(gids) if g is not None]
    ug=sorted(set(gids[i] for i in ok),key=lambda z:int(z)); gi={g:i for i,g in enumerate(ug)}
    agg=sparse.csr_matrix((np.ones(len(ok)),([gi[gids[i]] for i in ok],ok)),shape=(len(ug),len(names)))@mat
    nz=np.asarray(agg.sum(axis=1)).ravel()>0; agg=agg[nz]; ug=np.asarray(ug)[nz]
    dense=agg[:,[bidx[x[0]] for x in selected]].toarray().astype(float)
    ranks=stats.rankdata(dense,axis=0,method="average")/dense.shape[0]
    gi={g:i for i,g in enumerate(ug)}; rois=np.asarray([x[1] for x in selected])
    q=qc.loc[gsm]
    for m,fr in MEM.groupby("module_id",sort=True):
        for r in fr.itertuples(index=False):
            if r.gene_id not in gi: continue
            v=ranks[gi[r.gene_id]]
            delta=float(v[rois=="Nageotte"].mean()-v[rois=="Neuronal"].mean())
            contrib.append({"geo_accession":gsm,"sample_id":q.sample_id,"donor_id":q.donor_id,"DRG_ID":q.DRG_ID,"program":LABELS[m],"module_id":m,"gene_id":r.gene_id,"direction":r.direction,"section_rank_delta":delta,"signed_score_contribution":delta if r.direction=="up" else -delta})
con=pd.DataFrame(contrib)
con.to_csv(OUT/"M03_SECTION_GENE_CONTRIBUTIONS.tsv.gz",sep="\t",index=False,compression="gzip")
drg=con.groupby(["donor_id","DRG_ID","program","module_id","gene_id","direction"],as_index=False)[["section_rank_delta","signed_score_contribution"]].mean()
dg=drg.groupby(["donor_id","program","module_id","gene_id","direction"],as_index=False)[["section_rank_delta","signed_score_contribution"]].mean()
dg.to_csv(OUT/"M03_DONOR_GENE_CONTRIBUTIONS.tsv.gz",sep="\t",index=False,compression="gzip")
summ=dg.groupby(["program","module_id","gene_id","direction"],as_index=False).agg(mean_rank_delta=("section_rank_delta","mean"),mean_signed_score_contribution=("signed_score_contribution","mean"),n_donors=("donor_id","nunique"))
summ["absolute_contribution"]=summ.mean_signed_score_contribution.abs()
summ.to_csv(OUT/"M03_GENE_CONTRIBUTION_SUMMARY.tsv.gz",sep="\t",index=False,compression="gzip")

cd=pd.DataFrame(coords); cd.to_csv(OUT/"M03_REAL_COORDINATES_WITH_AUTHOR_ROI.tsv.gz",sep="\t",index=False,compression="gzip")
valid=cd[cd.roi.isin(["Nageotte","Neuronal"])].copy()
counts=valid.groupby("geo_accession").size().sort_index(); representative=counts.index[0]

def coordinate_panel(ax,d,title=None):
    colors={"Nageotte":"#C44E52","Neuronal":"#4C72B0","mixed_excluded":"#999999","other_excluded":"#D9D9D9"}
    for lab,g in d.groupby("roi"):
        ax.scatter(g.pixel_col,g.pixel_row,s=12,c=colors.get(lab,"#D9D9D9"),label=lab,alpha=.85,linewidths=0)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    if title: ax.set_title(title,fontsize=8)

fig,ax=plt.subplots(figsize=(5,4)); coordinate_panel(ax,cd[cd.geo_accession==representative]); ax.legend(frameon=False,fontsize=7,loc="best"); fig.tight_layout(); fig.savefig(C/"07_figures/M03_representative_real_coordinates.png",dpi=450); plt.close(fig)
gsms=sorted(cd.geo_accession.unique()); fig,axes=plt.subplots(4,4,figsize=(10,10))
for ax,gsm in zip(axes.flat,gsms): coordinate_panel(ax,cd[cd.geo_accession==gsm],gsm)
fig.tight_layout(); fig.savefig(C/"07_figures/M03_all_sections_real_coordinates.png",dpi=450); plt.close(fig)

order=[f"P{i}" for i in range(1,11) if i!=9]
fig,axes=plt.subplots(1,2,figsize=(10,4.5))
for i,p in enumerate(order):
    x=don[don.program==p].S_delta.to_numpy(); axes[0].scatter(np.repeat(i,len(x)),x,s=25,color="#333333",alpha=.8)
    axes[0].plot([i-.22,i+.22],[x.mean(),x.mean()],color="#C44E52",lw=2)
axes[0].axhline(0,color="black",lw=.7); axes[0].set_xticks(range(len(order)),order); axes[0].set_ylabel("Nageotte minus neuronal score")
pivot=arm.pivot(index="program",columns="arm",values="effect").reindex(order)
x=np.arange(len(order)); axes[1].bar(x-.18,pivot.source_up,.36,color="#4C72B0",label="ΔU"); axes[1].bar(x+.18,-pivot.source_down,.36,color="#DD8452",label="−ΔD")
axes[1].axhline(0,color="black",lw=.7); axes[1].set_xticks(x,order); axes[1].set_ylabel("Contribution to ΔS"); axes[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(C/"07_figures/M03_donor_distribution_and_arms.png",dpi=450); plt.close(fig)

(OUT/"M03_CLOSEOUT_SUMMARY.json").write_text(json.dumps({"donors":len(donors),"sections":int(qc.sample_id.nunique()),"programs_evaluable":9,"two_arm_tests":len(arm),"gene_contribution_rows":len(con),"representative_section_rule":"lexicographically first section among sections with valid author ROI coordinates; not selected by effect","representative_section":representative,"distance_gradient":"NOT_EVALUABLE: no histology image, scalefactors_json, verified spot diameter or lesion boundary among 80 tar members","diagnosis_boundary":"donors have diabetes history; explicit DPN diagnosis unavailable"},indent=2),encoding="utf-8")
print("M03 closeout complete")
