from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


P = Path(r"E:/1 Codex project/00_Cross_Project_Integration/osfu_redesign_preflight_2026-08-25/phase0_6_human_dpn_stage_projection")
OLD = P / "analysis_v4/all_extensions_2026-09-11"
O = P / "analysis_v4/targeted_repair_2026-09-12"
STRICT = P / "analysis_v3/gene_definition_repair_v1/strict_background_addendum"
MEMBERS = pd.read_csv(STRICT / "inputs/REPAIRED_MEMBERS.tsv", sep="\t", dtype={"gene_id": str})
GENE_INFO = P / "data/raw/NCBI_orthology_2026-08-27/Homo_sapiens.gene_info.gz"
NETWORK = P / "analysis_v3/batch02/results/network/MODULE_MEMBERS.tsv"
SEED = 20260912
LABELS = {
    "original_early_allcell": "P1", "original_late_allcell": "P2", "original_late_neuron": "P3",
    "original_severity": "P4", "original_xenium": "P5", "late_shared_concordant_neuronal_core": "P6",
    "late_neuron_residual": "P7", "late_allcell_residual": "P8",
    "severity_neuron_shared_concordant_core": "P9", "severity_neuron_residual": "P10",
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def put(frame: pd.DataFrame, rel: str):
    path = O / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def bh(values):
    x = pd.to_numeric(values, errors="coerce")
    out = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    if ok.any():
        out[ok] = multipletests(x[ok], method="fdr_bh")[1]
    return out


def build_symbol_resolver():
    info = pd.read_csv(GENE_INFO, sep="\t", dtype=str, keep_default_na=False)[["GeneID", "Symbol", "Synonyms"]]
    exact_candidates, folded_candidates, synonym_candidates = defaultdict(set), defaultdict(set), defaultdict(set)
    id_to_symbol = {}
    for row in info.itertuples(index=False):
        gid, sym = str(row.GeneID), str(row.Symbol)
        exact_candidates[sym].add(gid); folded_candidates[sym.upper()].add(gid); id_to_symbol[gid] = sym
        if row.Synonyms and row.Synonyms != "-":
            for alias in str(row.Synonyms).split("|"):
                alias = alias.strip()
                if alias and alias != "-": synonym_candidates[alias.upper()].add(gid)
    exact = {k: next(iter(v)) for k, v in exact_candidates.items() if len(v) == 1}
    folded = {k: next(iter(v)) for k, v in folded_candidates.items() if len(v) == 1}
    synonym = {k: next(iter(v)) for k, v in synonym_candidates.items() if len(v) == 1 and k not in folded}
    return exact, folded, synonym, id_to_symbol, info


def resolve_symbol(value, resolver):
    exact, folded, synonym, id_to_symbol, _ = resolver
    symbol = str(value).strip()
    if symbol in exact: gid, method = exact[symbol], "unique_official_exact"
    elif symbol.upper() in folded: gid, method = folded[symbol.upper()], "unique_official_casefold"
    elif symbol.upper() in synonym: gid, method = synonym[symbol.upper()], "unique_synonym"
    else: return None, None, "unresolved_or_ambiguous"
    return gid, id_to_symbol.get(gid), method


def allocations(n, k):
    total = math.comb(n, k)
    if total <= 100000:
        for chosen in itertools.combinations(range(n), k):
            y = np.zeros(n, int); y[list(chosen)] = 1; yield y
    else:
        rng = np.random.default_rng(SEED + n + k)
        for _ in range(9999):
            y = np.zeros(n, int); y[rng.choice(n, k, replace=False)] = 1; yield y


def perm_p(values, labels):
    values, labels = np.asarray(values, float), np.asarray(labels, int)
    obs = values[labels == 1].mean() - values[labels == 0].mean()
    null = np.asarray([values[a == 1].mean() - values[a == 0].mean() for a in allocations(len(labels), int(labels.sum()))])
    exact = math.comb(len(labels), int(labels.sum())) <= 100000
    hits = int((np.abs(null) >= abs(obs) - 1e-12).sum())
    p = hits / len(null) if exact else (hits + 1) / (len(null) + 1)
    return float(obs), float(p), bool(exact), int(len(null))


def bootstrap_ci(values, labels, key):
    rng = np.random.default_rng(int(hashlib.sha256((str(key) + str(SEED)).encode()).hexdigest()[:8], 16))
    values, labels = np.asarray(values, float), np.asarray(labels, int)
    i0, i1 = np.where(labels == 0)[0], np.where(labels == 1)[0]
    delta = np.empty(2000)
    for i in range(2000):
        delta[i] = values[rng.choice(i1, len(i1), True)].mean() - values[rng.choice(i0, len(i0), True)].mean()
    return np.quantile(delta, [0.025, 0.975])


def sex_stratified_p(values, labels, sex):
    values, labels, sex = np.asarray(values, float), np.asarray(labels, int), np.asarray(sex, str)
    groups = [np.where(sex == s)[0] for s in sorted(set(sex))]
    choices = [list(itertools.combinations(g, int(labels[g].sum()))) for g in groups]
    total = math.prod(len(x) for x in choices)
    obs = values[labels == 1].mean() - values[labels == 0].mean()
    if total <= 100000:
        null = []
        for selected in itertools.product(*choices):
            y = np.zeros(len(labels), int); y[[i for group in selected for i in group]] = 1
            null.append(values[y == 1].mean() - values[y == 0].mean())
        return float(np.mean(np.abs(np.asarray(null)) >= abs(obs) - 1e-12)), True, int(total)
    rng, hits = np.random.default_rng(SEED + 77), 0
    for _ in range(9999):
        y = np.zeros(len(labels), int)
        for g in groups: y[rng.choice(g, int(labels[g].sum()), replace=False)] = 1
        hits += abs(values[y == 1].mean() - values[y == 0].mean()) >= abs(obs) - 1e-12
    return (hits + 1) / 10000, False, int(total)


def score_matrix(expr: pd.DataFrame, resource: str):
    rank = expr.rank(axis=0, method="average", pct=True)
    genes, score_rows, coverage = set(rank.index.astype(str)), [], []
    for module, frame in MEMBERS.groupby("module_id", sort=True):
        up = frame.loc[frame.direction.eq("up") & frame.gene_id.isin(genes), "gene_id"].tolist()
        down = frame.loc[frame.direction.eq("down") & frame.gene_id.isin(genes), "gene_id"].tolist()
        eligible = len(up) >= 10 and len(down) >= 10
        coverage.append({"resource": resource, "module_id": module, "measured_up": len(up), "measured_down": len(down),
                         "total_up": int(frame.direction.eq("up").sum()), "total_down": int(frame.direction.eq("down").sum()),
                         "eligible": eligible, "reason": "" if eligible else "fewer than 10 measured genes in at least one arm"})
        if not eligible: continue
        U, D = rank.loc[up].mean(), rank.loc[down].mean()
        for sample in expr.columns:
            score_rows.append({"resource": resource, "module_id": module, "sample_id": sample,
                               "U": U[sample], "D": D[sample], "S": U[sample] - D[sample]})
    return pd.DataFrame(score_rows), pd.DataFrame(coverage), rank


def group_effects(scores, metadata, resource, comparisons, sex_sensitivity=False):
    merged, rows = scores.merge(metadata, on="sample_id", validate="many_to_one"), []
    for comparison, query, group_col, positive in comparisons:
        eligible_samples = metadata.sample_id if query == "ALL" else metadata.query(query).sample_id
        sub = merged[merged.sample_id.isin(eligible_samples)]
        for module, frame in sub.groupby("module_id"):
            labels = frame[group_col].eq(positive).astype(int).to_numpy()
            for component in ("S", "U", "D"):
                effect, p, exact, permutations = perm_p(frame[component], labels)
                lo, hi = bootstrap_ci(frame[component], labels, (resource, comparison, module, component))
                row = {"resource": resource, "comparison": comparison, "module_id": module, "component": component,
                       "n_positive": int(labels.sum()), "n_negative": int((1 - labels).sum()),
                       "effect_positive_minus_negative": effect, "ci_low": lo, "ci_high": hi, "p": p,
                       "exact": exact, "permutations": permutations, "status": "COMPLETED"}
                if sex_sensitivity:
                    ps, exs, space = sex_stratified_p(frame[component], labels, frame.sex.fillna("missing"))
                    row.update({"p_sex_stratified": ps, "sex_stratified_exact": exs, "sex_stratified_space": space})
                rows.append(row)
    out = pd.DataFrame(rows)
    out["q_family"] = np.nan
    full = out.component.eq("S")
    out.loc[full, "q_family"] = bh(out.loc[full, "p"])
    out.loc[~full, "q_family"] = bh(out.loc[~full, "p"])
    return out, merged


def paired_m05(scores, meta):
    merged = scores.merge(meta, on="sample_id", validate="many_to_one")
    wide = merged.pivot(index=["module_id", "donor", "group"], columns="muscle", values=["S", "U", "D"]).dropna().reset_index()
    donor_rows = []
    for row in wide.itertuples(index=False, name=None):
        record = dict(zip(wide.columns, row))
        module, donor, group = record[("module_id", "")], record[("donor", "")], record[("group", "")]
        donor_rows.append({"module_id": module, "donor": donor, "group": group,
                           **{f"{c}_AH_minus_MG": record[(c, "AH")] - record[(c, "MG")] for c in ("S", "U", "D")}})
    donors, tests = pd.DataFrame(donor_rows), []
    for module, frame in donors.groupby("module_id"):
        labels = frame.group.eq("DPN").astype(int).to_numpy()
        for component in ("S", "U", "D"):
            values = frame[f"{component}_AH_minus_MG"]
            effect, p, exact, permutations = perm_p(values, labels)
            lo, hi = bootstrap_ci(values, labels, ("paired", module, component))
            tests.append({"comparison": "DPN minus MNC difference of AH minus MG", "module_id": module,
                          "component": component, "n_DPN": int(labels.sum()), "n_MNC": int((1-labels).sum()),
                          "effect": effect, "ci_low": lo, "ci_high": hi, "p": p, "exact": exact,
                          "permutations": permutations, "status": "COMPLETED_SMALL_SAMPLE"})
    tests = pd.DataFrame(tests); tests["q_family"] = np.nan
    full = tests.component.eq("S"); tests.loc[full, "q_family"] = bh(tests.loc[full, "p"]); tests.loc[~full, "q_family"] = bh(tests.loc[~full, "p"])
    return donors, tests


def module_block_map(resolver):
    modules = pd.read_csv(NETWORK, sep="\t", dtype=str)
    assert list(modules.columns) == ["gene", "module"], "MODULE_MEMBERS.tsv schema changed"
    resolved = modules.gene.map(lambda x: resolve_symbol(x, resolver))
    modules[["gene_id", "current_symbol", "mapping_method"]] = pd.DataFrame(resolved.tolist(), index=modules.index)
    modules["is_non_grey"] = modules.module.ne("grey")
    modules["candidate"] = modules.gene_id.notna() & modules.is_non_grey
    conflicts = set(modules.loc[modules.candidate].groupby("gene_id").module.nunique().loc[lambda x: x > 1].index)
    modules["include_in_block_map"] = modules.candidate & ~modules.gene_id.isin(conflicts)
    modules["exclusion_reason"] = np.select(
        [modules.gene_id.isna(), modules.module.eq("grey"), modules.gene_id.isin(conflicts)],
        ["unresolved_or_ambiguous_symbol", "grey_not_a_network_block", "GeneID_maps_to_multiple_modules"], default="included")
    block_map = modules.loc[modules.include_in_block_map].drop_duplicates(["gene_id", "module"]).set_index("gene_id").module.to_dict()
    return modules, block_map


def gene_and_block_influence(rank, meta, comparisons, resource, block_map):
    gene_rows, block_rows = [], []
    for comparison, query, group_col, positive in comparisons:
        m = meta if query == "ALL" else meta.query(query)
        samples = [s for s in m.sample_id if s in rank.columns]
        labels = m.set_index("sample_id").loc[samples, group_col].eq(positive).astype(int)
        for module, members in MEMBERS.groupby("module_id", sort=True):
            measured = members[members.gene_id.isin(rank.index)].copy()
            up, down = measured.loc[measured.direction.eq("up"), "gene_id"].tolist(), measured.loc[measured.direction.eq("down"), "gene_id"].tolist()
            if len(up) < 10 or len(down) < 10: continue
            delta = rank.loc[measured.gene_id, samples].T.groupby(labels.to_numpy()).mean().T
            delta = delta[1] - delta[0]
            du, dd = float(delta.loc[up].mean()), float(delta.loc[down].mean()); full = du - dd
            measured["block"] = measured.gene_id.map(block_map)
            for row in measured.itertuples():
                value = float(delta.loc[row.gene_id]); n = len(up) if row.direction == "up" else len(down)
                signed = value / n if row.direction == "up" else -value / n
                gene_rows.append({"resource": resource, "comparison": comparison, "module_id": module, "gene_id": row.gene_id,
                                  "symbol": row.symbol, "direction": row.direction, "block": row.block if pd.notna(row.block) else "unassigned",
                                  "delta_percentile_rank": value, "signed_additive_contribution": signed, "delta_S_full": full})
            for block, block_members in measured.dropna(subset=["block"]).groupby("block"):
                remove = set(block_members.gene_id); remaining_up = [g for g in up if g not in remove]; remaining_down = [g for g in down if g not in remove]
                if len(remaining_up) >= 10 and len(remaining_down) >= 10:
                    after = float(delta.loc[remaining_up].mean() - delta.loc[remaining_down].mean()); status = "COMPLETED"
                else: after, status = np.nan, "NOT_EVALUABLE"
                block_rows.append({"resource": resource, "comparison": comparison, "module_id": module, "block": block,
                                   "removed_genes": len(remove), "remaining_up": len(remaining_up), "remaining_down": len(remaining_down),
                                   "delta_S_full": full, "delta_S_after_block_removal": after,
                                   "change_in_delta_S": after-full if np.isfinite(after) else np.nan, "status": status})
    return pd.DataFrame(gene_rows), pd.DataFrame(block_rows)


def m01_block_influence(block_map):
    gene = pd.read_csv(OLD / "02_M01_arms/results/GENE_CONTRIBUTIONS.tsv.gz", sep="\t", dtype={"gene_id": str})
    gene["block_repaired"] = gene.gene_id.map(block_map)
    rows = []
    for (endpoint, module), frame in gene.groupby(["endpoint", "module_id"]):
        up = frame.loc[frame.direction.eq("up"), "gene_id"].tolist(); down = frame.loc[frame.direction.eq("down"), "gene_id"].tolist()
        delta = frame.drop_duplicates("gene_id").set_index("gene_id").delta_rank.astype(float)
        full = float(delta.loc[up].mean() - delta.loc[down].mean())
        for block, bf in frame.dropna(subset=["block_repaired"]).groupby("block_repaired"):
            remove = set(bf.gene_id); ru = [g for g in up if g not in remove]; rd = [g for g in down if g not in remove]
            if len(ru) >= 10 and len(rd) >= 10:
                after, status = float(delta.loc[ru].mean() - delta.loc[rd].mean()), "COMPLETED"
            else: after, status = np.nan, "NOT_EVALUABLE"
            rows.append({"resource": endpoint, "comparison": endpoint, "module_id": module, "block": block,
                         "removed_genes": len(remove), "remaining_up": len(ru), "remaining_down": len(rd),
                         "delta_S_full": full, "delta_S_after_block_removal": after,
                         "change_in_delta_S": after-full if np.isfinite(after) else np.nan, "status": status})
    return pd.DataFrame(rows)


def plot_effects(frame, path, title):
    data = frame[frame.component.eq("S")].copy(); data["program"] = data.module_id.map(LABELS)
    comparisons = list(data.comparison.unique())
    fig, axes = plt.subplots(1, len(comparisons), figsize=(6.4 * len(comparisons), 4.6), squeeze=False)
    for ax, comparison in zip(axes.ravel(), comparisons):
        part = data[data.comparison.eq(comparison)].sort_values("program", key=lambda x: x.str.extract(r"(\d+)")[0].astype(int))
        y = np.arange(len(part)); effect = part.effect_positive_minus_negative.to_numpy(float)
        ax.errorbar(effect, y, xerr=[effect-part.ci_low.to_numpy(float), part.ci_high.to_numpy(float)-effect], fmt="o", color="#2f5d8a", ecolor="#8097ad", capsize=2)
        ax.axvline(0, color="#555555", lw=1); ax.set_yticks(y, part.program); ax.invert_yaxis(); ax.set_title(comparison); ax.set_xlabel("Corrected fixed score difference")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(title, y=1.01, fontsize=11); fig.tight_layout()
    for suffix, dpi in [("png", 300), ("tiff", 600), ("pdf", None)]: fig.savefig(path.with_suffix("."+suffix), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    resolver = build_symbol_resolver()
    modules, block_map = module_block_map(resolver)
    put(modules, "03_module_block_repair/results/MODULE_SYMBOL_TO_GENEID_AUDIT.tsv")
    summary = modules.groupby(["module", "exclusion_reason"], dropna=False).size().reset_index(name="rows")
    put(summary, "03_module_block_repair/results/MODULE_MAPPING_SUMMARY.tsv")

    metadata = pd.read_csv(OLD / "00_admin/PUBLIC_SAMPLE_METADATA.tsv", sep="\t")

    # M05: resolve all raw rows independently of program membership and rank over the complete eligible universe.
    raw5 = pd.read_csv(OLD / "01_inputs/public_downloads/GSE143979/GSE143979_merged_gene_name_expression.txt.gz", sep="\t")
    sample_cols5 = [c for c in raw5.columns if c not in {"name", "length"}]
    assert len(sample_cols5) == 15 and "length" not in sample_cols5
    resolved5 = raw5.name.map(lambda x: resolve_symbol(x, resolver))
    map5 = raw5[["name", "length"]].copy(); map5[["gene_id", "current_symbol", "mapping_method"]] = pd.DataFrame(resolved5.tolist(), index=raw5.index)
    numeric5 = raw5[sample_cols5].apply(pd.to_numeric, errors="coerce")
    map5["finite_all_libraries"] = np.isfinite(numeric5).all(axis=1); map5["all_zero"] = numeric5.eq(0).all(axis=1)
    map5["include_row"] = map5.gene_id.notna() & map5.finite_all_libraries & ~map5.all_zero
    map5["exclusion_reason"] = np.select([map5.gene_id.isna(), ~map5.finite_all_libraries, map5.all_zero], ["unresolved_or_ambiguous", "nonfinite", "all_zero"], default="included")
    map5.to_csv(O / "02_M05_background_repair/results/M05_ROW_MAPPING.tsv.gz", sep="\t", index=False)
    kept5 = numeric5.loc[map5.include_row].copy(); kept5["gene_id"] = map5.loc[map5.include_row, "gene_id"].to_numpy()
    x5 = kept5.groupby("gene_id", sort=True)[sample_cols5].median()
    x5_de = kept5.groupby("gene_id", sort=True)[sample_cols5].sum()
    colmap = {}
    import re
    for col in x5.columns:
        match = re.match(r"(mb\d+)_(ah|mg\d*)_", col.lower())
        if match: colmap[col] = match.group(1).upper() + " " + ("MG" if match.group(2).startswith("mg") else "AH")
    x5 = x5.rename(columns=colmap)
    x5_de = x5_de.rename(columns=colmap)
    m5 = metadata[metadata.resource.eq("GSE143979")].copy()
    m5meta = pd.DataFrame({"sample_id": m5.title, "donor": m5.donor_key, "muscle": m5.tissue_key,
                           "diagnosis": m5.diagnosis, "geo_accession": m5.geo_accession})
    m5meta["group"] = np.where(m5meta.diagnosis.str.startswith("DPN"), "DPN", "MNC")
    assert set(x5.columns) == set(m5meta.sample_id)
    x5.to_csv(O / "01_inputs/M05_ELIGIBLE_COUNTS.tsv.gz", sep="\t", index_label="gene_id")
    x5_de.to_csv(O / "01_inputs/M05_DE_COUNTS_SUM.tsv.gz", sep="\t", index_label="gene_id")
    put(m5meta, "01_inputs/M05_SAMPLE_METADATA.tsv")
    scores5, coverage5, rank5 = score_matrix(x5, "GSE143979")
    comparisons5 = [("AH DPN minus MNC", "muscle == 'AH'", "group", "DPN"), ("MG DPN minus MNC", "muscle == 'MG'", "group", "DPN")]
    effects5, _ = group_effects(scores5, m5meta, "GSE143979", comparisons5)
    donors5, paired5 = paired_m05(scores5, m5meta)
    genes5, blocks5 = gene_and_block_influence(rank5, m5meta, comparisons5, "GSE143979", block_map)
    put(scores5, "02_M05_background_repair/results/M05_CORRECTED_SCORES.tsv"); put(coverage5, "02_M05_background_repair/results/M05_CORRECTED_COVERAGE.tsv")
    put(effects5, "02_M05_background_repair/results/M05_CORRECTED_EFFECTS_BY_MUSCLE.tsv"); put(donors5, "02_M05_background_repair/results/M05_CORRECTED_PAIRED_DIFFERENCES.tsv")
    put(paired5, "02_M05_background_repair/results/M05_CORRECTED_PAIRED_INTERACTION.tsv")
    genes5.to_csv(O / "02_M05_background_repair/results/M05_CORRECTED_GENE_CONTRIBUTIONS.tsv.gz", sep="\t", index=False)
    put(blocks5, "03_module_block_repair/results/M05_CORRECTED_BLOCK_INFLUENCE.tsv")
    background5 = pd.DataFrame([{
        "raw_feature_rows": len(raw5), "rna_library_columns": len(sample_cols5), "resolved_rows": int(map5.gene_id.notna().sum()),
        "all_zero_rows": int(map5.all_zero.sum()), "included_rows_before_duplicate_aggregation": int(map5.include_row.sum()),
        "eligible_unique_geneids": len(x5), "program_geneids_in_background": len(set(x5.index) & set(MEMBERS.gene_id)),
        "nonprogram_geneids_in_background": len(set(x5.index) - set(MEMBERS.gene_id)), "historical_program_limited_background": 3728,
        "aggregation": "median_by_GeneID", "ranking": "average_percentile_over_all_eligible_geneids"}])
    put(background5, "02_M05_background_repair/results/M05_BACKGROUND_AUDIT.tsv")

    old_scores5 = pd.read_csv(OLD / "06_M05_muscle/results/M05_SCORES.tsv.gz", sep="\t")
    old_effects5 = pd.read_csv(OLD / "06_M05_muscle/results/M05_EFFECTS_BY_MUSCLE.tsv", sep="\t")
    score_cmp5 = old_scores5.merge(scores5, on=["resource", "module_id", "sample_id"], suffixes=("_old", "_corrected"))
    for c in ("U", "D", "S"): score_cmp5[f"delta_{c}_corrected_minus_old"] = score_cmp5[f"{c}_corrected"] - score_cmp5[f"{c}_old"]
    effect_cmp5 = old_effects5.merge(effects5, on=["resource", "comparison", "module_id", "component"], suffixes=("_old", "_corrected"))
    effect_cmp5["delta_effect_corrected_minus_old"] = effect_cmp5.effect_positive_minus_negative_corrected - effect_cmp5.effect_positive_minus_negative_old
    put(score_cmp5, "02_M05_background_repair/results/M05_OLD_NEW_SCORE_COMPARISON.tsv"); put(effect_cmp5, "02_M05_background_repair/results/M05_OLD_NEW_EFFECT_COMPARISON.tsv")

    # M04: audit general Ensembl mapping and repair zero/duplicate handling only because the prior path violates the locked strict-background rule.
    raw4 = pd.read_csv(OLD / "01_inputs/public_downloads/GSE250152/GSE250152_ens_genes_toc_all.csv.gz", index_col=0)
    original_ids = raw4.index.astype(str); base_ids = original_ids.str.replace(r"\.\d+$", "", regex=True)
    mapping = pd.read_csv(P / "results/tables/GSE302658_ensembl_to_ncbi_gene_mapping_2026-08-27.tsv.gz", sep="\t", dtype=str)
    mapping["one_to_one"] = mapping.ensembl_mapping_is_one_to_one.eq("True")
    map_unique = mapping[mapping.one_to_one].groupby("ensembl_gene_id").filter(lambda x: x.human_gene_id.nunique() == 1).drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
    audit4 = pd.DataFrame({"raw_ensembl_id": original_ids, "ensembl_gene_id": base_ids})
    audit4["gene_id"] = audit4.ensembl_gene_id.map(map_unique.human_gene_id); audit4["current_symbol"] = audit4.ensembl_gene_id.map(map_unique.current_symbol)
    audit4["finite_all_libraries"] = np.isfinite(raw4.to_numpy()).all(axis=1); audit4["all_zero"] = raw4.eq(0).all(axis=1).to_numpy()
    audit4["include_row"] = audit4.gene_id.notna() & audit4.finite_all_libraries & ~audit4.all_zero
    audit4["exclusion_reason"] = np.select([audit4.gene_id.isna(), ~audit4.finite_all_libraries, audit4.all_zero], ["not_unique_one_to_one_mapping", "nonfinite", "all_zero"], default="included")
    audit4.to_csv(O / "01_inputs/M04_ENSEMBL_ROW_MAPPING_AUDIT.tsv.gz", sep="\t", index=False)
    keep4 = audit4.include_row.to_numpy(bool)
    kept4 = raw4.iloc[keep4].copy(); kept4["gene_id"] = audit4.loc[keep4, "gene_id"].to_numpy()
    x4 = kept4.groupby("gene_id", sort=True).median(numeric_only=True)
    x4_de = kept4.groupby("gene_id", sort=True).sum(numeric_only=True)
    m4 = metadata[metadata.resource.eq("GSE250152")].copy()
    m4meta = pd.DataFrame({"sample_id": m4.title, "group": m4.disease_state, "sex": m4.sex, "tissue": m4.tissue, "geo_accession": m4.geo_accession})
    x4.to_csv(O / "01_inputs/M04_ELIGIBLE_COUNTS.tsv.gz", sep="\t", index_label="gene_id")
    x4_de.to_csv(O / "01_inputs/M04_DE_COUNTS_SUM.tsv.gz", sep="\t", index_label="gene_id")
    put(m4meta, "01_inputs/M04_SAMPLE_METADATA.tsv")
    scores4, coverage4, rank4 = score_matrix(x4, "GSE250152")
    comparisons4 = [("Morton minus control", "ALL", "group", "Morton’s neuroma")]
    effects4, _ = group_effects(scores4, m4meta, "GSE250152", comparisons4, sex_sensitivity=True)
    genes4, blocks4 = gene_and_block_influence(rank4, m4meta, comparisons4, "GSE250152", block_map)
    put(scores4, "01_inputs/M04_CORRECTED_SCORES.tsv"); put(coverage4, "01_inputs/M04_CORRECTED_COVERAGE.tsv"); put(effects4, "01_inputs/M04_CORRECTED_EFFECTS.tsv")
    genes4.to_csv(O / "01_inputs/M04_CORRECTED_GENE_CONTRIBUTIONS.tsv.gz", sep="\t", index=False); put(blocks4, "03_module_block_repair/results/M04_CORRECTED_BLOCK_INFLUENCE.tsv")
    background4 = pd.DataFrame([{"raw_feature_rows": len(raw4), "mapped_unique_one_to_one_rows": int(audit4.gene_id.notna().sum()),
                                 "all_zero_rows": int(audit4.all_zero.sum()), "included_rows_before_duplicate_aggregation": int(audit4.include_row.sum()),
                                 "eligible_unique_geneids": len(x4), "mapping_rows_total": len(mapping),
                                 "mapping_present_in_GSE302658_false": int(mapping.present_in_GSE302658.eq("False").sum()),
                                 "mapping_conclusion": "general annotation table, not a blood-expression subset",
                                 "repair_trigger": "historical path retained all-zero rows and summed duplicate GeneIDs"}])
    put(background4, "01_inputs/M04_MAPPING_AND_BACKGROUND_AUDIT.tsv")
    old_scores4 = pd.read_csv(OLD / "05_M04_morton/results/M04_SCORES.tsv.gz", sep="\t")
    old_effects4 = pd.read_csv(OLD / "05_M04_morton/results/M04_EFFECTS.tsv", sep="\t")
    score_cmp4 = old_scores4.merge(scores4, on=["resource", "module_id", "sample_id"], suffixes=("_old", "_corrected"))
    for c in ("U", "D", "S"): score_cmp4[f"delta_{c}_corrected_minus_old"] = score_cmp4[f"{c}_corrected"] - score_cmp4[f"{c}_old"]
    effect_cmp4 = old_effects4.merge(effects4, on=["resource", "comparison", "module_id", "component"], suffixes=("_old", "_corrected"))
    effect_cmp4["delta_effect_corrected_minus_old"] = effect_cmp4.effect_positive_minus_negative_corrected - effect_cmp4.effect_positive_minus_negative_old
    put(score_cmp4, "01_inputs/M04_OLD_NEW_SCORE_COMPARISON.tsv"); put(effect_cmp4, "01_inputs/M04_OLD_NEW_EFFECT_COMPARISON.tsv")

    # M01 block influence only; score decomposition is otherwise untouched.
    blocks1 = m01_block_influence(block_map); put(blocks1, "03_module_block_repair/results/M01_CORRECTED_BLOCK_INFLUENCE.tsv")

    # M03 display-only section counts from distinct sample IDs at the DRG level.
    old_map3 = pd.read_csv(OLD / "04_M03_spatial/results/M03_DONOR_DRG_SECTION_MAP.tsv", sep="\t")
    corrected3 = old_map3.groupby(["Donor ID", "Donor #", "DRG ID", "Sex", "Age"], as_index=False).agg(Sections=("Sample ID", "nunique"), Sample_IDs=("Sample ID", lambda x: "|".join(sorted(set(x)))))
    put(corrected3, "04_M03_count_repair/results/M03_CORRECTED_DONOR_DRG_SECTION_COUNT.tsv")
    count_check = pd.DataFrame([{"distinct_sections": old_map3["Sample ID"].nunique(), "distinct_DRGs": old_map3["DRG ID"].nunique(),
                                "distinct_donors": old_map3["Donor ID"].nunique(), "sum_corrected_DRG_sections": int(corrected3.Sections.sum()),
                                "historical_display_sum": int(old_map3.drop_duplicates(["Donor ID", "DRG ID"])["# Samples"].sum())}])
    put(count_check, "04_M03_count_repair/results/M03_COUNT_RECONCILIATION.tsv")

    # Cross-module sensitivity: replace only the historical M04/M05 primary records, retain scope of 133 tests.
    cross = pd.read_csv(OLD / "10_integration/CROSS_MODULE_PRIMARY_BH_SENSITIVITY.tsv", sep="\t")
    repaired_rows = []
    for mod, effects in [("M04", effects4), ("M05", effects5)]:
        for row in effects[effects.component.eq("S")].itertuples():
            repaired_rows.append({"module": mod, "endpoint": row.comparison, "program": row.module_id, "test": "full score", "p": row.p})
    repaired_rows = pd.DataFrame(repaired_rows)
    retained = cross[~cross.module.isin(["M04", "M05"])][["module", "endpoint", "program", "test", "p"]]
    updated_cross = pd.concat([retained, repaired_rows], ignore_index=True)
    assert len(updated_cross) == 133
    updated_cross["q_BH_cross_module_sensitivity"] = bh(updated_cross.p)
    put(cross, "08_integration/HISTORICAL_CROSS_MODULE_PRIMARY_BH_SENSITIVITY.tsv")
    put(updated_cross, "08_integration/CORRECTED_CROSS_MODULE_PRIMARY_BH_SENSITIVITY.tsv")

    plot_effects(effects5, O / "02_M05_background_repair/figures/M05_corrected_program_effects", "GSE143979 corrected full-background scores")
    plot_effects(effects4, O / "02_M05_background_repair/figures/M04_corrected_program_effects", "GSE250152 corrected strict-background scores")

    # Semantic checks, including independent score and BH recomputation.
    program_gene_ids = set(MEMBERS.gene_id)
    tests = []
    def check(name, passed, observed, expected): tests.append({"check": name, "pass": bool(passed), "observed": observed, "expected": expected})
    check("M05_has_15_RNA_libraries", len(sample_cols5) == 15, len(sample_cols5), 15)
    check("M05_length_not_counted_as_library", "length" not in sample_cols5, "length" in sample_cols5, False)
    check("M05_background_contains_nonprogram_genes", len(set(x5.index)-program_gene_ids) > 0, len(set(x5.index)-program_gene_ids), ">0")
    check("M05_background_exceeds_program_gene_count", len(x5) > len(set(x5.index)&program_gene_ids), f"{len(x5)} vs {len(set(x5.index)&program_gene_ids)}", "background > program genes")
    check("M05_background_all_finite", np.isfinite(x5.to_numpy()).all(), True, True)
    check("M05_background_no_all_zero_GeneIDs", not x5.eq(0).all(axis=1).any(), int(x5.eq(0).all(axis=1).sum()), 0)
    check("M05_unique_GeneID_rows", x5.index.is_unique, x5.index.duplicated().sum(), 0)
    check("M05_P9_remains_ineligible", not bool(coverage5.loc[coverage5.module_id.eq("severity_neuron_shared_concordant_core"), "eligible"].iloc[0]), False, False)
    sample = scores5.iloc[0]; members = MEMBERS[MEMBERS.module_id.eq(sample.module_id)]; up = members.loc[members.direction.eq("up") & members.gene_id.isin(rank5.index), "gene_id"]; down = members.loc[members.direction.eq("down") & members.gene_id.isin(rank5.index), "gene_id"]
    independent = float(rank5.loc[up, sample.sample_id].mean() - rank5.loc[down, sample.sample_id].mean())
    check("M05_sample_score_independent_recalculation", abs(independent-sample.S) < 1e-12, abs(independent-sample.S), "<1e-12")
    q5 = bh(effects5.loc[effects5.component.eq("S"), "p"]); check("M05_joint_full_score_BH_recomputed", np.allclose(q5, effects5.loc[effects5.component.eq("S"), "q_family"]), float(np.max(abs(q5-effects5.loc[effects5.component.eq("S"), "q_family"]))), 0)
    check("module_schema_exact", list(pd.read_csv(NETWORK, sep="\t", nrows=1).columns) == ["gene", "module"], list(pd.read_csv(NETWORK, sep="\t", nrows=1).columns), ["gene", "module"])
    check("module_non_grey_GeneID_mapping_nonempty", len(block_map) > 0, len(block_map), ">0")
    completed_blocks = sum((d.status == "COMPLETED").sum() for d in [blocks1, blocks4, blocks5])
    check("block_deletion_has_evaluable_results", completed_blocks > 0, int(completed_blocks), ">0")
    check("M03_sections_reconcile", int(corrected3.Sections.sum()) == 16, int(corrected3.Sections.sum()), 16)
    check("M03_DRGs_reconcile", len(corrected3) == 7, len(corrected3), 7)
    check("M03_donors_reconcile", corrected3["Donor ID"].nunique() == 6, corrected3["Donor ID"].nunique(), 6)
    check("cross_module_family_size_preserved", len(updated_cross) == 133, len(updated_cross), 133)
    put(pd.DataFrame(tests), "tests/SEMANTIC_TESTS.tsv")

    status = pd.DataFrame([
        {"task": "restore_configuration_template", "status": "COMPLETED", "detail": "7,566-byte attachment copied with receipt-time hash; not backdated"},
        {"task": "restore_SOURCE_VERIFICATION_NOTES", "status": "BLOCKED_RESOURCE", "detail": "content described in user message but original file bytes not attached or located"},
        {"task": "M05_full_background_score_repair", "status": "COMPLETED", "detail": "strict canonical mapping, all-zero exclusion, median GeneID aggregation, full eligible rank universe"},
        {"task": "M04_mapping_background_audit_and_repair", "status": "COMPLETED", "detail": "mapping is general; repaired all-zero retention and duplicate summation"},
        {"task": "network_block_mapping_and_deletion", "status": "COMPLETED", "detail": "explicit gene/module schema and canonical symbol mapping; no WGCNA rebuild"},
        {"task": "M03_section_count_display", "status": "COMPLETED", "detail": "distinct Sample ID per DRG"},
        {"task": "whole_gene_DE", "status": "RUNNING", "detail": "awaiting verified project-local R/Bioconductor environment"},
        {"task": "M02_input_recovery", "status": "PLANNED", "detail": "resume/official-source audit after targeted P0 calculations"},
        {"task": "M03_input_recovery", "status": "PLANNED", "detail": "single official archive strategy after targeted P0 calculations"},
    ])
    put(status, "STATUS.tsv")
    run = {"completed_utc": datetime.now(timezone.utc).isoformat(), "semantic_tests": len(tests), "semantic_tests_pass": int(pd.DataFrame(tests)["pass"].sum()),
           "M05_background_geneids": len(x5), "M05_nonprogram_geneids": len(set(x5.index)-program_gene_ids),
           "network_block_geneids": len(block_map), "evaluable_block_deletions": int(completed_blocks),
           "M04_full_q_lt_005": int((effects4.loc[effects4.component.eq("S"), "q_family"] < 0.05).sum()),
           "M05_full_q_lt_005": int((effects5.loc[effects5.component.eq("S"), "q_family"] < 0.05).sum()),
           "cross_module_q_lt_005": int((updated_cross.q_BH_cross_module_sensitivity < 0.05).sum()), "exit_code": 0}
    (O / "logs/00_run_targeted_repairs.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(json.dumps(run, indent=2))


if __name__ == "__main__":
    main()
