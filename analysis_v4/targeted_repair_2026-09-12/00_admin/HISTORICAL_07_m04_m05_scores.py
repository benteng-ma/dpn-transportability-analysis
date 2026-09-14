from __future__ import annotations

import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


P = Path(r"E:/1 Codex project/00_Cross_Project_Integration/osfu_redesign_preflight_2026-08-25/phase0_6_human_dpn_stage_projection")
O = P / "analysis_v4/all_extensions_2026-09-11"
S = P / "analysis_v3/gene_definition_repair_v1/strict_background_addendum"
MEMBERS = pd.read_csv(S / "inputs/REPAIRED_MEMBERS.tsv", sep="\t", dtype={"gene_id": str})
RNG_SEED = 20260912
LABELS={'original_early_allcell':'P1','original_late_allcell':'P2','original_late_neuron':'P3','original_severity':'P4','original_xenium':'P5','late_shared_concordant_neuronal_core':'P6','late_neuron_residual':'P7','late_allcell_residual':'P8','severity_neuron_shared_concordant_core':'P9','severity_neuron_residual':'P10'}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""): h.update(b)
    return h.hexdigest()


def score_matrix(expr: pd.DataFrame, resource: str):
    rank = expr.rank(axis=0, method="average", pct=True)
    rows, cov = [], []
    genes = set(rank.index.astype(str))
    for mid, f in MEMBERS.groupby("module_id", sort=True):
        up = f.loc[f.direction.eq("up") & f.gene_id.isin(genes), "gene_id"].tolist()
        down = f.loc[f.direction.eq("down") & f.gene_id.isin(genes), "gene_id"].tolist()
        eligible = len(up) >= 10 and len(down) >= 10
        cov.append({"resource": resource, "module_id": mid, "measured_up": len(up), "measured_down": len(down),
                    "total_up": int(f.direction.eq("up").sum()), "total_down": int(f.direction.eq("down").sum()), "eligible": eligible,
                    "reason": "" if eligible else "fewer than 10 measured genes in at least one arm"})
        if not eligible: continue
        U, D = rank.loc[up].mean(), rank.loc[down].mean()
        for sample in expr.columns:
            rows.append({"resource": resource, "module_id": mid, "sample_id": sample, "U": U[sample], "D": D[sample], "S": U[sample]-D[sample]})
    return pd.DataFrame(rows), pd.DataFrame(cov), rank


def all_allocations(n: int, k: int):
    if math.comb(n, k) <= 100000:
        for chosen in itertools.combinations(range(n), k):
            y = np.zeros(n, int); y[list(chosen)] = 1; yield y
    else:
        rng = np.random.default_rng(RNG_SEED+n+k)
        for _ in range(9999):
            y = np.zeros(n, int); y[rng.choice(n, k, replace=False)] = 1; yield y


def perm_p(values, y):
    values, y = np.asarray(values, float), np.asarray(y, int)
    obs = values[y == 1].mean() - values[y == 0].mean()
    null = np.array([values[a == 1].mean() - values[a == 0].mean() for a in all_allocations(len(y), int(y.sum()))])
    exact = math.comb(len(y), int(y.sum())) <= 100000
    hits = int((np.abs(null) >= abs(obs)-1e-12).sum())
    return obs, hits/len(null) if exact else (hits+1)/(len(null)+1), exact, len(null)


def sex_stratified_p(values, y, sex):
    values, y, sex = np.asarray(values,float), np.asarray(y,int), np.asarray(sex,str)
    groups = [np.where(sex == s)[0] for s in sorted(set(sex))]
    combos = [list(itertools.combinations(g, int(y[g].sum()))) for g in groups]
    total = math.prod(len(c) for c in combos)
    obs = values[y==1].mean()-values[y==0].mean()
    if total <= 100000:
        null=[]
        for cc in itertools.product(*combos):
            yy=np.zeros(len(y),int); yy[[i for z in cc for i in z]]=1
            null.append(values[yy==1].mean()-values[yy==0].mean())
        null=np.asarray(null); return float(np.mean(np.abs(null)>=abs(obs)-1e-12)), True, total
    rng=np.random.default_rng(RNG_SEED+77); hits=0
    for _ in range(9999):
        yy=np.zeros(len(y),int)
        for g in groups: yy[rng.choice(g,int(y[g].sum()),replace=False)]=1
        d=values[yy==1].mean()-values[yy==0].mean();hits+=abs(d)>=abs(obs)-1e-12
    return (hits+1)/10000,False,total


def boot_ci(values, y, key):
    rng=np.random.default_rng(int(hashlib.sha256((str(key)+str(RNG_SEED)).encode()).hexdigest()[:8],16))
    a=np.asarray(values,float); y=np.asarray(y,int); i0=np.where(y==0)[0];i1=np.where(y==1)[0]
    d=np.empty(2000)
    for i in range(2000):d[i]=a[rng.choice(i1,len(i1),True)].mean()-a[rng.choice(i0,len(i0),True)].mean()
    return np.quantile(d,[.025,.975])


def add_bh(d, pcol, qcol, family=None):
    d[qcol]=np.nan
    groups=[(None,d)] if family is None else d.groupby(family,dropna=False)
    for _,g in groups:
        ok=g[pcol].notna()
        if ok.any(): d.loc[g.index[ok],qcol]=multipletests(g.loc[ok,pcol],method="fdr_bh")[1]


def contributions(rank, meta, resource, comparisons):
    mods=pd.read_csv(P/'analysis_v3/batch02/results/network/MODULE_MEMBERS.tsv',sep='\t',dtype=str)
    gidcol='gene_id' if 'gene_id' in mods else ('human_gene_id' if 'human_gene_id' in mods else mods.columns[0]);modcol='module' if 'module' in mods else ('module_color' if 'module_color' in mods else mods.columns[-1]);blockmap=dict(zip(mods[gidcol].astype(str),mods[modcol].astype(str)))
    gene_rows=[];block_rows=[]
    for cname,sub,group_col,positive in comparisons:
        sm=meta.copy() if sub=='ALL' else meta.query(sub).copy();samples=[s for s in sm.sample_id if s in rank.columns];sm=sm.set_index('sample_id').loc[samples];y=sm[group_col].eq(positive).astype(int)
        for mid,f in MEMBERS.groupby('module_id',sort=True):
            fm=f[f.gene_id.isin(rank.index)].copy();up=fm.loc[fm.direction.eq('up'),'gene_id'].tolist();down=fm.loc[fm.direction.eq('down'),'gene_id'].tolist()
            if len(up)<10 or len(down)<10:continue
            dr=rank.loc[fm.gene_id,samples].T.groupby(y).mean().T;delta=dr[1]-dr[0]
            du=float(delta.loc[up].mean());dd=float(delta.loc[down].mean());full=du-dd
            fm['block']=fm.gene_id.map(blockmap).fillna('unassigned')
            for r in fm.itertuples():
                v=float(delta.loc[r.gene_id]);n=len(up) if r.direction=='up' else len(down);con=v/n if r.direction=='up' else -v/n;change=(du-v)/(len(up)-1) if r.direction=='up' else (v-dd)/(len(down)-1)
                gene_rows.append({'resource':resource,'comparison':cname,'module_id':mid,'gene_id':r.gene_id,'symbol':r.symbol,'direction':r.direction,'block':r.block,'delta_percentile_rank':v,'signed_additive_contribution':con,'delta_S_full':full,'change_in_delta_S_leave_one_gene_out':change})
            for bl,bf in fm.groupby('block'):
                ru=[g for g in up if g not in set(bf.gene_id)];rd=[g for g in down if g not in set(bf.gene_id)]
                if len(ru)>=10 and len(rd)>=10:
                    after=float(delta.loc[ru].mean()-delta.loc[rd].mean());status='COMPLETED'
                else:after=np.nan;status='NOT_EVALUABLE'
                block_rows.append({'resource':resource,'comparison':cname,'module_id':mid,'block':bl,'removed_genes':len(bf),'remaining_up':len(ru),'remaining_down':len(rd),'delta_S_full':full,'delta_S_after_block_removal':after,'change_in_delta_S':after-full if np.isfinite(after) else np.nan,'status':status})
    return pd.DataFrame(gene_rows),pd.DataFrame(block_rows)


def effects(scores, meta, resource, comparisons):
    out=[]
    merged=scores.merge(meta,on="sample_id",how="left",validate="many_to_one")
    for cname, sub, group_col, positive in comparisons:
        submeta=meta.copy() if sub == "ALL" else meta.query(sub).copy(); mm=merged[merged.sample_id.isin(submeta.sample_id)]
        for mid,g in mm.groupby("module_id"):
            y=g[group_col].eq(positive).astype(int).to_numpy()
            if len(set(y))<2: continue
            for comp in ["S","U","D"]:
                delta,p,exact,nperm=perm_p(g[comp],y);lo,hi=boot_ci(g[comp],y,(resource,cname,mid,comp))
                row={"resource":resource,"comparison":cname,"module_id":mid,"component":comp,"n_positive":int(y.sum()),"n_negative":int((1-y).sum()),
                     "effect_positive_minus_negative":delta,"ci_low":lo,"ci_high":hi,"p":p,"exact":exact,"permutations":nperm,"status":"COMPLETED"}
                if resource=="GSE250152" and "sex" in g:
                    ps,exs,space=sex_stratified_p(g[comp],y,g.sex.fillna("missing"));row.update({"p_sex_stratified":ps,"sex_stratified_exact":exs,"sex_stratified_space":space})
                out.append(row)
    return pd.DataFrame(out),merged


def main():
    for d in [O/'05_M04_morton/results', O/'05_M04_morton/figures', O/'06_M05_muscle/results', O/'06_M05_muscle/figures']:
        d.mkdir(parents=True, exist_ok=True)
    locked={"locked_utc":datetime.now(timezone.utc).isoformat(),"known_after_v3_1":True,"members_sha256":sha(S/'inputs/REPAIRED_MEMBERS.tsv'),
            "score":"within-sample percentile ranks; mean(up) minus mean(down); at least 10 genes per arm",
            "inference":"same 2000 group-bootstrap and fixed-label two-sided permutation implementation as M01, per instruction section 4H",
            "M04_primary":"Morton neuroma minus control nerves; full-score programs one BH family; U/D secondary one BH family; sex-stratified permutation sensitivity; anatomy remains fully confounded",
            "M05_primary":"DPN including severe DPN minus MNC separately in AH and MG; all muscle by program full-score tests one BH family; U/D secondary family",
            "M05_paired":"AH minus MG donor differences, compared DPN minus MNC; separate all-program family; at least 3 per group remains small-sample exploration",
            "DE_status":"not run: full companion config is missing and no usable R command-line/locked count-engine environment exists; counts were not converted or rounded",
            "no_outcome_tuning":True}
    (O/'00_admin/LOCKED_PROTOCOL_M04_M05.json').write_text(json.dumps(locked,indent=2),encoding='utf-8')

    # M04: Ensembl IDs to unambiguous NCBI GeneIDs, duplicate mappings summed before rank scoring.
    m4p=O/'01_inputs/public_downloads/GSE250152/GSE250152_ens_genes_toc_all.csv.gz'
    x4=pd.read_csv(m4p,index_col=0);x4.index=x4.index.astype(str).str.replace(r'\.\d+$','',regex=True)
    mp=pd.read_csv(P/'results/tables/GSE302658_ensembl_to_ncbi_gene_mapping_2026-08-27.tsv.gz',sep='\t',dtype=str)
    mp=mp[mp.ensembl_mapping_is_one_to_one.eq('True')].drop_duplicates('ensembl_gene_id').set_index('ensembl_gene_id')
    x4=x4.loc[x4.index.intersection(mp.index)].copy();x4['gene_id']=mp.loc[x4.index,'human_gene_id'].to_numpy();x4=x4.groupby('gene_id').sum(numeric_only=True)
    md=pd.read_csv(O/'00_admin/PUBLIC_SAMPLE_METADATA.tsv',sep='\t');m4=md[md.resource.eq('GSE250152')].copy()
    m4meta=pd.DataFrame({'sample_id':m4.title,'group':m4.disease_state,'sex':m4.sex,'tissue':m4.tissue,'geo_accession':m4.geo_accession})
    sc4,cov4,rank4=score_matrix(x4,'GSE250152')
    ef4,merged4=effects(sc4,m4meta,'GSE250152',[("Morton minus control","ALL","group","Morton’s neuroma")])
    gc4,bi4=contributions(rank4,m4meta,'GSE250152',[("Morton minus control","ALL","group","Morton’s neuroma")])
    add_bh(ef4[ef4.component.eq('S')].copy(),'p','q')
    sidx=ef4.component.eq('S');ef4.loc[sidx,'q_family']=multipletests(ef4.loc[sidx,'p'],method='fdr_bh')[1]
    aidx=ef4.component.ne('S');ef4.loc[aidx,'q_family']=multipletests(ef4.loc[aidx,'p'],method='fdr_bh')[1]
    sc4.to_csv(O/'05_M04_morton/results/M04_SCORES.tsv.gz',sep='\t',index=False);cov4.to_csv(O/'05_M04_morton/results/M04_COVERAGE.tsv',sep='\t',index=False);ef4.to_csv(O/'05_M04_morton/results/M04_EFFECTS.tsv',sep='\t',index=False)
    gc4.to_csv(O/'05_M04_morton/results/M04_GENE_CONTRIBUTIONS.tsv.gz',sep='\t',index=False);bi4.to_csv(O/'05_M04_morton/results/M04_BLOCK_INFLUENCE.tsv',sep='\t',index=False)
    m4meta.to_csv(O/'05_M04_morton/results/M04_SAMPLE_MAP.tsv',sep='\t',index=False)

    # M05: source rows are gene-symbol counts; resolve only unique current symbols from locked members.
    m5p=O/'01_inputs/public_downloads/GSE143979/GSE143979_merged_gene_name_expression.txt.gz'
    raw5=pd.read_csv(m5p,sep='\t');raw5=raw5.drop(columns=['length']);symmap=MEMBERS[['symbol','gene_id']].drop_duplicates();amb=symmap.groupby('symbol').gene_id.nunique();symmap=symmap[symmap.symbol.map(amb).eq(1)].drop_duplicates('symbol').set_index('symbol')
    raw5=raw5[raw5.name.isin(symmap.index)].copy();raw5['gene_id']=symmap.loc[raw5.name,'gene_id'].to_numpy();x5=raw5.drop(columns='name').groupby('gene_id').sum(numeric_only=True)
    # Convert deposited lower-case library names to direct GEO titles.
    m5=md[md.resource.eq('GSE143979')].copy(); colmap={}
    for c in x5.columns:
        z=c.lower();dm=re.match(r'(mb\d+)_(ah|mg\d*)_',z)
        if dm: colmap[c]=dm.group(1).upper()+' '+('MG' if dm.group(2).startswith('mg') else 'AH')
    x5=x5.rename(columns=colmap); assert set(x5.columns)==set(m5.title)
    m5meta=pd.DataFrame({'sample_id':m5.title,'donor':m5.donor_key,'muscle':m5.tissue_key,'diagnosis':m5.diagnosis,'geo_accession':m5.geo_accession})
    m5meta['group']=np.where(m5meta.diagnosis.str.startswith('DPN'),'DPN','MNC')
    sc5,cov5,rank5=score_matrix(x5,'GSE143979')
    ef5,merged5=effects(sc5,m5meta,'GSE143979',[("AH DPN minus MNC","muscle == 'AH'","group","DPN"),("MG DPN minus MNC","muscle == 'MG'","group","DPN")])
    gc5,bi5=contributions(rank5,m5meta,'GSE143979',[("AH DPN minus MNC","muscle == 'AH'","group","DPN"),("MG DPN minus MNC","muscle == 'MG'","group","DPN")])
    sidx=ef5.component.eq('S');ef5.loc[sidx,'q_family']=multipletests(ef5.loc[sidx,'p'],method='fdr_bh')[1]
    aidx=ef5.component.ne('S');ef5.loc[aidx,'q_family']=multipletests(ef5.loc[aidx,'p'],method='fdr_bh')[1]
    # Complete donor pairs: program-wise AH-MG difference, then DPN-MNC exact comparison.
    paired=[]
    z=merged5.pivot(index=['module_id','donor','group'],columns='muscle',values=['S','U','D']).dropna().reset_index()
    for _,r in z.iterrows():
        paired.append({'module_id':r['module_id'].iloc[0] if hasattr(r['module_id'],'iloc') else r['module_id'],'donor':r['donor'].iloc[0] if hasattr(r['donor'],'iloc') else r['donor'],'group':r['group'].iloc[0] if hasattr(r['group'],'iloc') else r['group'],
                       'S_AH_minus_MG':r[('S','AH')]-r[('S','MG')],'U_AH_minus_MG':r[('U','AH')]-r[('U','MG')],'D_AH_minus_MG':r[('D','AH')]-r[('D','MG')]})
    pairdf=pd.DataFrame(paired);pe=[]
    for mid,g in pairdf.groupby('module_id'):
        y=g.group.eq('DPN').astype(int).to_numpy()
        for comp in ['S','U','D']:
            d,pv,ex,npv=perm_p(g[f'{comp}_AH_minus_MG'],y);lo,hi=boot_ci(g[f'{comp}_AH_minus_MG'],y,('paired',mid,comp))
            pe.append({'comparison':'DPN minus MNC difference of AH minus MG','module_id':mid,'component':comp,'n_DPN':int(y.sum()),'n_MNC':int((1-y).sum()),'effect':d,'ci_low':lo,'ci_high':hi,'p':pv,'exact':ex,'permutations':npv,'status':'COMPLETED_SMALL_SAMPLE'})
    pe=pd.DataFrame(pe); sidx=pe.component.eq('S');pe.loc[sidx,'q_family']=multipletests(pe.loc[sidx,'p'],method='fdr_bh')[1];aidx=~sidx;pe.loc[aidx,'q_family']=multipletests(pe.loc[aidx,'p'],method='fdr_bh')[1]
    sc5.to_csv(O/'06_M05_muscle/results/M05_SCORES.tsv.gz',sep='\t',index=False);cov5.to_csv(O/'06_M05_muscle/results/M05_COVERAGE.tsv',sep='\t',index=False);ef5.to_csv(O/'06_M05_muscle/results/M05_EFFECTS_BY_MUSCLE.tsv',sep='\t',index=False);pairdf.to_csv(O/'06_M05_muscle/results/M05_PAIRED_DIFFERENCES.tsv',sep='\t',index=False);pe.to_csv(O/'06_M05_muscle/results/M05_PAIRED_INTERACTION.tsv',sep='\t',index=False);m5meta.to_csv(O/'06_M05_muscle/results/M05_DONOR_MUSCLE_MAP.tsv',sep='\t',index=False)
    gc5.to_csv(O/'06_M05_muscle/results/M05_GENE_CONTRIBUTIONS.tsv.gz',sep='\t',index=False);bi5.to_csv(O/'06_M05_muscle/results/M05_BLOCK_INFLUENCE.tsv',sep='\t',index=False)

    # Compact empirical displays; full numeric sources remain TSV.
    for name,ef,folder in [('M04',ef4,O/'05_M04_morton/figures'),('M05',ef5,O/'06_M05_muscle/figures')]:
        folder.mkdir(parents=True,exist_ok=True);d=ef[ef.component.eq('S')].copy();d['lab']=d.module_id.map(LABELS)
        fig,axes=plt.subplots(1,d.comparison.nunique(),figsize=(6.8*max(1,d.comparison.nunique()),4.8),squeeze=False)
        for ax,(cmp,g) in zip(axes.ravel(),d.groupby('comparison')):
            yy=np.arange(len(g));ax.errorbar(g.effect_positive_minus_negative,yy,xerr=[g.effect_positive_minus_negative-g.ci_low,g.ci_high-g.effect_positive_minus_negative],fmt='o',color='#2f5d8a',ecolor='#8097ad',capsize=2);ax.axvline(0,color='#555',lw=1);ax.set_yticks(yy,g.lab);ax.invert_yaxis();ax.set_title(cmp);ax.set_xlabel('Fixed score difference (rank scale)');ax.spines[['top','right']].set_visible(False)
        fig.tight_layout()
        for ext,dpi in [('png',300),('tiff',600),('pdf',None)]:fig.savefig(folder/f'{name}_Figure_1_program_effects.{ext}',dpi=dpi,bbox_inches='tight')
        plt.close(fig)
    summary={'utc':datetime.now(timezone.utc).isoformat(),'M04_eligible_programs':int(cov4.eligible.sum()),'M04_full_score_q_lt_0_05':int((ef4.loc[ef4.component.eq('S'),'q_family']<.05).sum()),'M05_eligible_programs':int(cov5.eligible.sum()),'M05_full_score_q_lt_0_05':int((ef5.loc[ef5.component.eq('S'),'q_family']<.05).sum()),'M05_paired_n':int(pairdf.donor.nunique()),'M05_paired_q_lt_0_05':int((pe.loc[pe.component.eq('S'),'q_family']<.05).sum())}
    (O/'logs/07_m04_m05_scores.json').write_text(json.dumps({**summary,'exit_code':0},indent=2),encoding='utf-8');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    import re
    main()
