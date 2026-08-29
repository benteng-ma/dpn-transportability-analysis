#!/usr/bin/env python3
"""Project frozen human hDRG disease-stage signatures into GSE176017 animals."""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
RAW_NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
SIGNATURES = TABLES / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"
COUNTS = TABLES / "GSE176017_animal_pseudobulk_raw_counts_2026-08-27.tsv.gz"
ANIMAL_META = PHASE / "metadata" / "GSE176017_animal_metadata_2026-08-27.tsv"
RANDOM_SEED = 20260827
N_NULL = 2000

PRIMARY_TESTS = {
    "early_allcell_diabetes_vs_control": ("Normal", "Diabetes_no_allodynia"),
    "late_allcell_DPN_vs_diabetes": ("Diabetes_no_allodynia", "Painful_DPN"),
    "late_neuron_DPN_vs_diabetes": ("Diabetes_no_allodynia", "Painful_DPN"),
    "severity_neuron_modhigh_vs_low_nageotte": ("Diabetes_no_allodynia", "Painful_DPN"),
}
SECONDARY_TESTS = {
    "xenium_DPN_vs_control": ("Normal", "Painful_DPN"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_gene_info(path: Path) -> pd.DataFrame:
    info = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, na_filter=False)
    return info[["#tax_id", "GeneID", "Symbol", "Synonyms"]].copy()


def build_symbol_lookup(info: pd.DataFrame) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    official_exact: dict[str, str] = {}
    official_folded_candidates: defaultdict[str, set[str]] = defaultdict(set)
    synonym_candidates: defaultdict[str, set[str]] = defaultdict(set)
    for row in info.itertuples(index=False):
        gene_id = str(row.GeneID)
        symbol = str(row.Symbol)
        official_exact[symbol] = gene_id
        official_folded_candidates[symbol.upper()].add(gene_id)
        synonyms = str(row.Synonyms)
        if synonyms and synonyms != "-":
            for synonym in synonyms.split("|"):
                synonym = synonym.strip()
                if synonym and synonym != "-":
                    synonym_candidates[synonym.upper()].add(gene_id)
    official_folded = {
        symbol: next(iter(ids)) for symbol, ids in official_folded_candidates.items() if len(ids) == 1
    }
    synonym_unique = {
        symbol: next(iter(ids))
        for symbol, ids in synonym_candidates.items()
        if len(ids) == 1 and symbol not in official_folded
    }
    return official_exact, official_folded, synonym_unique


def resolve_symbol(
    symbol: str,
    official_exact: dict[str, str],
    official_folded: dict[str, str],
    synonym_unique: dict[str, str],
) -> tuple[str | None, str]:
    if symbol in official_exact:
        return official_exact[symbol], "official_exact"
    folded = symbol.upper()
    if folded in official_folded:
        return official_folded[folded], "official_casefold"
    if folded in synonym_unique:
        return synonym_unique[folded], "unique_synonym"
    return None, "unresolved_or_ambiguous"


def load_one_to_one_orthologues(path: Path) -> tuple[dict[str, str], dict]:
    pairs: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n").split("\t")
        expected = ["#tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID"]
        if header != expected:
            raise RuntimeError(f"Unexpected orthologue header: {header}")
        for line in handle:
            tax_id, gene_id, relationship, other_tax, other_gene = line.rstrip("\n").split("\t")
            if relationship != "Ortholog":
                continue
            if tax_id == "9606" and other_tax == "10116":
                pairs.add((gene_id, other_gene))
            elif tax_id == "10116" and other_tax == "9606":
                pairs.add((other_gene, gene_id))
    human_to_rat_all: defaultdict[str, set[str]] = defaultdict(set)
    rat_to_human_all: defaultdict[str, set[str]] = defaultdict(set)
    for human, rat in pairs:
        human_to_rat_all[human].add(rat)
        rat_to_human_all[rat].add(human)
    one_to_one = {
        human: next(iter(rats))
        for human, rats in human_to_rat_all.items()
        if len(rats) == 1 and len(rat_to_human_all[next(iter(rats))]) == 1
    }
    qc = {
        "human_rat_unique_pairs": len(pairs),
        "human_gene_ids_with_any_rat_orthologue": len(human_to_rat_all),
        "rat_gene_ids_with_any_human_orthologue": len(rat_to_human_all),
        "reciprocal_one_to_one_pairs": len(one_to_one),
    }
    return one_to_one, qc


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    values = pvalues.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return pd.Series(result, index=pvalues.index)


def exact_label_permutation(values: pd.Series, groups: pd.Series, early: str, later: str) -> tuple[float, int]:
    subset = values[groups.isin([early, later])]
    subset_groups = groups.loc[subset.index]
    n_later = int((subset_groups == later).sum())
    observed = float(subset[subset_groups == later].mean() - subset[subset_groups == early].mean())
    indices = list(range(len(subset)))
    array = subset.to_numpy(dtype=float)
    null = []
    for later_indices in itertools.combinations(indices, n_later):
        later_mask = np.zeros(len(array), dtype=bool)
        later_mask[list(later_indices)] = True
        null.append(float(array[later_mask].mean() - array[~later_mask].mean()))
    p_two = sum(abs(value) >= abs(observed) - 1e-12 for value in null) / len(null)
    return float(p_two), len(null)


def hedges_g(values: pd.Series, groups: pd.Series, early: str, later: str) -> float:
    x = values[groups == later].to_numpy(dtype=float)
    y = values[groups == early].to_numpy(dtype=float)
    if len(x) < 2 or len(y) < 2:
        return np.nan
    df = len(x) + len(y) - 2
    pooled_var = ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df
    if pooled_var <= 0:
        return np.nan
    d = (x.mean() - y.mean()) / math.sqrt(pooled_var)
    correction = 1 - 3 / (4 * df - 1) if df > 1 else np.nan
    return float(d * correction)


def leave_one_out(values: pd.Series, groups: pd.Series, early: str, later: str) -> tuple[bool, float, float]:
    keep = groups.isin([early, later])
    values = values[keep]
    groups = groups[keep]
    differences: list[float] = []
    for omitted in values.index:
        remaining_values = values.drop(index=omitted)
        remaining_groups = groups.drop(index=omitted)
        if (remaining_groups == early).sum() == 0 or (remaining_groups == later).sum() == 0:
            continue
        differences.append(
            float(
                remaining_values[remaining_groups == later].mean()
                - remaining_values[remaining_groups == early].mean()
            )
        )
    return bool(differences and all(value > 0 for value in differences)), min(differences), max(differences)


def matched_null_delta(
    rank_matrix: pd.DataFrame,
    group_series: pd.Series,
    early: str,
    later: str,
    up_ids: list[str],
    down_ids: list[str],
    deciles: pd.Series,
    signature_ids: set[str],
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    observed_scores = rank_matrix.loc[up_ids].mean(axis=0) - rank_matrix.loc[down_ids].mean(axis=0)
    observed_delta = float(
        observed_scores[group_series == later].mean() - observed_scores[group_series == early].mean()
    )
    background_by_decile: dict[int, np.ndarray] = {}
    for decile in sorted(deciles.dropna().unique()):
        candidates = deciles[(deciles == decile) & (~deciles.index.isin(signature_ids))].index.to_numpy()
        if len(candidates) == 0:
            candidates = deciles[deciles == decile].index.to_numpy()
        background_by_decile[int(decile)] = candidates

    def decile_counts(gene_ids: list[str]) -> dict[int, int]:
        return {int(key): int(value) for key, value in deciles.loc[gene_ids].value_counts().to_dict().items()}

    up_counts = decile_counts(up_ids)
    down_counts = decile_counts(down_ids)
    null_deltas = np.empty(N_NULL, dtype=float)
    for iteration in range(N_NULL):
        sampled_up: list[str] = []
        sampled_down: list[str] = []
        for decile, count in up_counts.items():
            pool = background_by_decile[decile]
            sampled_up.extend(rng.choice(pool, size=count, replace=count > len(pool)).tolist())
        for decile, count in down_counts.items():
            pool = background_by_decile[decile]
            sampled_down.extend(rng.choice(pool, size=count, replace=count > len(pool)).tolist())
        null_scores = rank_matrix.loc[sampled_up].mean(axis=0) - rank_matrix.loc[sampled_down].mean(axis=0)
        null_deltas[iteration] = float(
            null_scores[group_series == later].mean() - null_scores[group_series == early].mean()
        )
    p_directional = (1 + int((null_deltas >= observed_delta).sum())) / (N_NULL + 1)
    p_two_sided = (1 + int((np.abs(null_deltas) >= abs(observed_delta)).sum())) / (N_NULL + 1)
    return float(p_directional), float(p_two_sided), N_NULL


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    human_info = load_gene_info(RAW_NCBI / "Homo_sapiens.gene_info.gz")
    rat_info = load_gene_info(RAW_NCBI / "Rattus_norvegicus.gene_info.gz")
    h_exact, h_folded, h_synonym = build_symbol_lookup(human_info)
    r_exact, r_folded, r_synonym = build_symbol_lookup(rat_info)
    rat_id_to_symbol = rat_info.drop_duplicates("GeneID").set_index("GeneID")["Symbol"].to_dict()
    one_to_one, orth_qc = load_one_to_one_orthologues(RAW_NCBI / "gene_orthologs.gz")

    metadata = pd.read_csv(ANIMAL_META, sep="\t").set_index("sample_id")
    counts_symbol = pd.read_csv(COUNTS, sep="\t", index_col=0)
    counts_symbol = counts_symbol[metadata.index]
    library_sizes = counts_symbol.sum(axis=0)

    rat_gene_counts: defaultdict[str, np.ndarray] = defaultdict(
        lambda: np.zeros(counts_symbol.shape[1], dtype=np.int64)
    )
    rat_matrix_mapping_rows: list[dict] = []
    for symbol, row in counts_symbol.iterrows():
        gene_id, method = resolve_symbol(str(symbol), r_exact, r_folded, r_synonym)
        rat_matrix_mapping_rows.append(
            {"matrix_symbol": symbol, "rat_gene_id": gene_id, "rat_mapping_method": method}
        )
        if gene_id is not None:
            rat_gene_counts[gene_id] += row.to_numpy(dtype=np.int64)
    counts = pd.DataFrame.from_dict(rat_gene_counts, orient="index", columns=counts_symbol.columns)
    counts.index.name = "rat_gene_id"
    cpm = counts.divide(library_sizes, axis=1) * 1_000_000
    expressed = (cpm >= 1).sum(axis=1) >= 2

    rat_one_to_one_ids = set(one_to_one.values())
    background_ids = sorted(set(counts.index[expressed]) & rat_one_to_one_ids)
    log_cpm = np.log2(cpm.loc[background_ids] + 1)
    rank_matrix = log_cpm.rank(axis=0, method="average", pct=True) * 2 - 1
    mean_expression = log_cpm.mean(axis=1)
    deciles = pd.qcut(mean_expression.rank(method="first"), q=10, labels=False).astype(int)

    signatures = pd.read_csv(SIGNATURES, sep="\t")
    mapping_rows: list[dict] = []
    for row in signatures.itertuples(index=False):
        human_gene_id, human_method = resolve_symbol(str(row.gene), h_exact, h_folded, h_synonym)
        rat_gene_id = one_to_one.get(human_gene_id) if human_gene_id is not None else None
        mapping_rows.append(
            {
                "contrast_id": row.contrast_id,
                "human_source_symbol": row.gene,
                "human_gene_id": human_gene_id,
                "human_symbol_mapping_method": human_method,
                "human_direction": row.direction,
                "human_avg_log2FC": row.avg_log2FC,
                "rat_gene_id": rat_gene_id,
                "rat_official_symbol": rat_id_to_symbol.get(rat_gene_id) if rat_gene_id else None,
                "reciprocal_one_to_one": rat_gene_id is not None,
                "present_in_rat_matrix": rat_gene_id in counts.index if rat_gene_id else False,
                "passes_rat_expression_filter": rat_gene_id in background_ids if rat_gene_id else False,
            }
        )
    mapping = pd.DataFrame(mapping_rows)
    mapping_output = TABLES / "hDRG_to_rat_one_to_one_signature_mapping_2026-08-27.tsv"
    mapping.to_csv(mapping_output, sep="\t", index=False)
    pd.DataFrame(rat_matrix_mapping_rows).to_csv(
        TABLES / "GSE176017_rat_symbol_to_geneid_mapping_2026-08-27.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    sample_scores: list[dict] = []
    test_rows: list[dict] = []
    rng = np.random.default_rng(RANDOM_SEED)
    all_tests = {**PRIMARY_TESTS, **SECONDARY_TESTS}
    for contrast_id, (early, later) in all_tests.items():
        members = mapping[
            (mapping["contrast_id"] == contrast_id) & mapping["passes_rat_expression_filter"]
        ].copy()
        up = members[members["human_direction"] == "up"]
        down = members[members["human_direction"] == "down"]
        up_ids = up["rat_gene_id"].drop_duplicates().tolist()
        down_ids = down["rat_gene_id"].drop_duplicates().tolist()
        minimum_direction_size = 10 if contrast_id in PRIMARY_TESTS else 5
        if len(up_ids) < minimum_direction_size or len(down_ids) < minimum_direction_size:
            raise RuntimeError(
                f"Insufficient mapped genes for {contrast_id}: up={len(up_ids)}, down={len(down_ids)}, "
                f"minimum={minimum_direction_size}"
            )
        unweighted = rank_matrix.loc[up_ids].mean(axis=0) - rank_matrix.loc[down_ids].mean(axis=0)

        up_weights = up.drop_duplicates("rat_gene_id").set_index("rat_gene_id")["human_avg_log2FC"].abs().clip(upper=3)
        down_weights = down.drop_duplicates("rat_gene_id").set_index("rat_gene_id")["human_avg_log2FC"].abs().clip(upper=3)
        weighted_up = rank_matrix.loc[up_ids].multiply(up_weights, axis=0).sum(axis=0) / up_weights.sum()
        weighted_down = rank_matrix.loc[down_ids].multiply(down_weights, axis=0).sum(axis=0) / down_weights.sum()
        weighted = weighted_up - weighted_down

        for sample in metadata.index:
            sample_scores.append(
                {
                    "sample_id": sample,
                    "group": metadata.loc[sample, "group"],
                    "contrast_id": contrast_id,
                    "score_role": "primary" if contrast_id in PRIMARY_TESTS else "secondary",
                    "unweighted_rank_score": float(unweighted[sample]),
                    "weighted_rank_score": float(weighted[sample]),
                    "mapped_up_genes": len(up_ids),
                    "mapped_down_genes": len(down_ids),
                }
            )

        groups = metadata["group"]
        difference = float(unweighted[groups == later].mean() - unweighted[groups == early].mean())
        exact_p, permutation_count = exact_label_permutation(unweighted, groups, early, later)
        g = hedges_g(unweighted, groups, early, later)
        loo_stable, loo_min, loo_max = leave_one_out(unweighted, groups, early, later)
        null_directional, null_two_sided, null_count = matched_null_delta(
            rank_matrix,
            groups,
            early,
            later,
            up_ids,
            down_ids,
            deciles,
            set(up_ids) | set(down_ids),
            rng,
        )
        source_n = int((signatures["contrast_id"] == contrast_id).sum())
        test_rows.append(
            {
                "contrast_id": contrast_id,
                "test_role": "primary" if contrast_id in PRIMARY_TESTS else "secondary",
                "earlier_group": early,
                "later_group": later,
                "n_earlier": int((groups == early).sum()),
                "n_later": int((groups == later).sum()),
                "source_signature_genes": source_n,
                "mapped_expressed_up": len(up_ids),
                "mapped_expressed_down": len(down_ids),
                "mapping_coverage_fraction": (len(up_ids) + len(down_ids)) / source_n,
                "earlier_mean_score": float(unweighted[groups == early].mean()),
                "later_mean_score": float(unweighted[groups == later].mean()),
                "later_minus_earlier": difference,
                "hedges_g": g,
                "exact_two_sided_p": exact_p,
                "exact_permutation_count": permutation_count,
                "loo_all_positive": loo_stable,
                "loo_min_difference": loo_min,
                "loo_max_difference": loo_max,
                "matched_null_directional_p": null_directional,
                "matched_null_two_sided_p": null_two_sided,
                "matched_null_iterations": null_count,
            }
        )

    scores = pd.DataFrame(sample_scores)
    tests = pd.DataFrame(test_rows)
    primary_mask = tests["test_role"] == "primary"
    secondary_mask = tests["test_role"] == "secondary"
    tests.loc[primary_mask, "exact_BH_Q"] = bh_adjust(tests.loc[primary_mask, "exact_two_sided_p"])
    tests.loc[primary_mask, "matched_null_BH_Q"] = bh_adjust(
        tests.loc[primary_mask, "matched_null_directional_p"]
    )
    if secondary_mask.any():
        tests.loc[secondary_mask, "exact_BH_Q"] = tests.loc[secondary_mask, "exact_two_sided_p"]
        tests.loc[secondary_mask, "matched_null_BH_Q"] = tests.loc[
            secondary_mask, "matched_null_directional_p"
        ]
    tests["positive_direction"] = tests["later_minus_earlier"] > 0
    tests["effect_at_least_0_5"] = tests["hedges_g"] >= 0.5
    tests["matched_null_q_below_0_10"] = tests["matched_null_BH_Q"] < 0.10
    tests["provisionally_supportive"] = (
        tests["positive_direction"]
        & tests["effect_at_least_0_5"]
        & tests["loo_all_positive"]
        & tests["matched_null_q_below_0_10"]
    )

    early_pass = bool(
        tests.loc[tests["contrast_id"] == "early_allcell_diabetes_vs_control", "provisionally_supportive"].iloc[0]
    )
    neuron_late_pass = bool(
        tests.loc[
            tests["contrast_id"].isin(
                ["late_neuron_DPN_vs_diabetes", "severity_neuron_modhigh_vs_low_nageotte"]
            ),
            "provisionally_supportive",
        ].any()
    )
    gate_b = "PROVISIONALLY_SUPPORTIVE" if early_pass and neuron_late_pass else "FAIL"

    scores_output = TABLES / "GSE176017_human_stage_projection_scores_2026-08-27.tsv"
    tests_output = TABLES / "GSE176017_human_stage_projection_tests_2026-08-27.tsv"
    scores.to_csv(scores_output, sep="\t", index=False)
    tests.to_csv(tests_output, sep="\t", index=False)

    plot_contrasts = list(PRIMARY_TESTS)
    display_names = {
        "early_allcell_diabetes_vs_control": "Human early diabetes\n(all hDRG cells)",
        "late_allcell_DPN_vs_diabetes": "Human late DPN\n(all hDRG cells)",
        "late_neuron_DPN_vs_diabetes": "Human late DPN\n(neurons)",
        "severity_neuron_modhigh_vs_low_nageotte": "Human Nageotte severity\n(neurons)",
    }
    group_order = ["Normal", "Diabetes_no_allodynia", "Painful_DPN"]
    group_labels = ["Normal", "DM/no allodynia", "Painful DPN"]
    colors = ["#4C78A8", "#F2CF5B", "#D9534F"]
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.8), constrained_layout=True)
    jitter_rng = np.random.default_rng(RANDOM_SEED)
    for axis, contrast_id in zip(axes, plot_contrasts):
        subset = scores[scores["contrast_id"] == contrast_id]
        for x, (group, label, color) in enumerate(zip(group_order, group_labels, colors)):
            values = subset.loc[subset["group"] == group, "unweighted_rank_score"].to_numpy()
            jitter = jitter_rng.uniform(-0.08, 0.08, size=len(values))
            axis.scatter(np.full(len(values), x) + jitter, values, s=46, color=color, edgecolor="black", linewidth=0.55, zorder=3)
            axis.hlines(values.mean(), x - 0.20, x + 0.20, color="black", linewidth=2.0, zorder=4)
        stat = tests[tests["contrast_id"] == contrast_id].iloc[0]
        axis.axhline(0, color="#777777", linewidth=0.8, linestyle="--")
        axis.set_title(display_names[contrast_id], fontsize=10.5)
        axis.set_xticks(range(3), group_labels, rotation=28, ha="right", fontsize=8)
        axis.text(
            0.03,
            0.97,
            f"Δ={stat['later_minus_earlier']:.3f}\ng={stat['hedges_g']:.2f}\nexact P={stat['exact_two_sided_p']:.3f}\nLOO={'+' if stat['loo_all_positive'] else '−'}",
            transform=axis.transAxes,
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#CCCCCC", "alpha": 0.9},
        )
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Human-stage projection score\n(later disease state → positive)")
    fig.suptitle("Animal-level projection of human hDRG disease stages into GSE176017", fontsize=13, fontweight="bold")
    figure_png = FIGURES / "GSE176017_human_stage_projection_2026-08-27.png"
    figure_pdf = FIGURES / "GSE176017_human_stage_projection_2026-08-27.pdf"
    fig.savefig(figure_png, dpi=300, bbox_inches="tight")
    fig.savefig(figure_pdf, bbox_inches="tight")
    plt.close(fig)

    coverage = (
        mapping.groupby("contrast_id")
        .agg(
            source_genes=("human_source_symbol", "size"),
            human_symbol_resolved=("human_gene_id", lambda x: int(x.notna().sum())),
            reciprocal_one_to_one=("reciprocal_one_to_one", "sum"),
            rat_matrix_present=("present_in_rat_matrix", "sum"),
            rat_expressed=("passes_rat_expression_filter", "sum"),
        )
        .reset_index()
    )
    coverage["rat_expressed_fraction"] = coverage["rat_expressed"] / coverage["source_genes"]
    coverage_output = TABLES / "hDRG_to_GSE176017_signature_coverage_2026-08-27.tsv"
    coverage.to_csv(coverage_output, sep="\t", index=False)

    qc = {
        "status": "PASS",
        "gate_B": gate_b,
        "random_seed": RANDOM_SEED,
        "matched_null_iterations": N_NULL,
        "signature_sha256": sha256(SIGNATURES),
        "counts_sha256": sha256(COUNTS),
        "animal_metadata_sha256": sha256(ANIMAL_META),
        "ncbi_files": {
            name: {"path": str(RAW_NCBI / name), "sha256": sha256(RAW_NCBI / name)}
            for name in ["gene_orthologs.gz", "Homo_sapiens.gene_info.gz", "Rattus_norvegicus.gene_info.gz"]
        },
        "orthologue_qc": orth_qc,
        "rat_matrix_symbols": int(len(counts_symbol)),
        "rat_matrix_symbols_resolved_to_geneid": int(pd.DataFrame(rat_matrix_mapping_rows)["rat_gene_id"].notna().sum()),
        "rat_geneids_after_aggregation": int(len(counts)),
        "expressed_one_to_one_background_genes": len(background_ids),
        "early_transition_supportive": early_pass,
        "at_least_one_neuron_late_transition_supportive": neuron_late_pass,
        "exact_test_resolution_warning": (
            "There are 6 label assignments for n=2 vs n=2 (minimum absolute-statistic two-sided P=2/6=1/3) "
            "and 15 assignments for n=4 vs n=2 (P-value granularity=1/15). "
            "They are reported as small-sample uncertainty, not used as a conventional significance gate."
        ),
        "outputs": [
            str(mapping_output),
            str(scores_output),
            str(tests_output),
            str(coverage_output),
            str(figure_png),
            str(figure_pdf),
        ],
    }
    qc_output = TABLES / "GSE176017_human_stage_projection_qc_2026-08-27.json"
    qc_output.write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))
    print("\nPrimary tests:\n" + tests.to_string(index=False))


if __name__ == "__main__":
    main()
