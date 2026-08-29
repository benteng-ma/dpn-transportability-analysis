#!/usr/bin/env python3
"""Validate frozen human hDRG stage signatures in an independent human DPN bulk cohort."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PHASE = Path(__file__).resolve().parents[2]
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
METADATA = PHASE / "metadata"
RAW = PHASE / "data" / "raw" / "human_DPN_bulk_PMC8933403"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
SIGNATURES = TABLES / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"
DATE = "2026-08-27"
RANDOM_SEED = 20260827
N_NULL = 10_000

PRIMARY = {
    "late_neuron_DPN_vs_diabetes",
    "severity_neuron_modhigh_vs_low_nageotte",
}
CONTEXTUAL = {
    "early_allcell_diabetes_vs_control",
    "late_allcell_DPN_vs_diabetes",
    "xenium_DPN_vs_control",
}
DISPLAY = {
    "early_allcell_diabetes_vs_control": "Early all-cell",
    "late_allcell_DPN_vs_diabetes": "Late all-cell",
    "late_neuron_DPN_vs_diabetes": "Late neuron",
    "severity_neuron_modhigh_vs_low_nageotte": "Neuron severity",
    "xenium_DPN_vs_control": "Xenium spatial",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    array = values.to_numpy(dtype=float)
    order = np.argsort(array)
    ranked = array[order]
    adjusted = ranked * len(array) / np.arange(1, len(array) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return pd.Series(result, index=values.index)


def load_gene_info(path: Path) -> pd.DataFrame:
    info = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, na_filter=False)
    return info[["GeneID", "Symbol", "Synonyms"]].copy()


def build_symbol_lookup(
    info: pd.DataFrame,
) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    official_exact: dict[str, str] = {}
    official_folded_candidates: defaultdict[str, set[str]] = defaultdict(set)
    synonym_candidates: defaultdict[str, set[str]] = defaultdict(set)
    id_to_symbol: dict[str, str] = {}
    for row in info.itertuples(index=False):
        gene_id = str(row.GeneID)
        symbol = str(row.Symbol)
        official_exact[symbol] = gene_id
        official_folded_candidates[symbol.upper()].add(gene_id)
        id_to_symbol[gene_id] = symbol
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
    return official_exact, official_folded, synonym_unique, id_to_symbol


def repair_excel_gene(value: object) -> tuple[str, bool]:
    if isinstance(value, (datetime, date, pd.Timestamp)):
        month = int(value.month)
        day = int(value.day)
        if month == 3:
            return f"MARCH{day}", True
        if month == 9:
            return f"SEPT{day}", True
        return str(value), True
    return str(value).strip(), False


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


def add_resolution(
    frame: pd.DataFrame,
    gene_column: str,
    lookups: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> pd.DataFrame:
    official_exact, official_folded, synonym_unique, id_to_symbol = lookups
    repaired = frame[gene_column].map(repair_excel_gene)
    frame = frame.copy()
    frame["gene_input"] = [item[0] for item in repaired]
    frame["gene_recovered_from_excel_date"] = [item[1] for item in repaired]
    resolved = frame["gene_input"].map(
        lambda symbol: resolve_symbol(symbol, official_exact, official_folded, synonym_unique)
    )
    frame["human_gene_id"] = [item[0] for item in resolved]
    frame["mapping_method"] = [item[1] for item in resolved]
    frame["current_human_symbol"] = frame["human_gene_id"].map(id_to_symbol)
    return frame


def load_target(
    lookups: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    de_path = RAW / "41598_2022_8100_MOESM5_ESM.xlsx"
    expr_path = RAW / "41598_2022_8100_MOESM6_ESM.xlsx"
    de = pd.read_excel(de_path, sheet_name="S4 All_genes", engine="openpyxl")
    de = de.rename(columns={de.columns[0]: "source_gene"})
    de["source_row"] = np.arange(2, len(de) + 2)
    de = add_resolution(de, "source_gene", lookups)
    for column in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "weight"]:
        de[column] = pd.to_numeric(de[column], errors="coerce")
    de["finite_stat"] = np.isfinite(de["stat"])
    de["padj_sort"] = de["padj"].fillna(np.inf)
    de["abs_stat"] = de["stat"].abs().fillna(-np.inf)
    resolved_de = de[de["human_gene_id"].notna()].copy()
    resolved_de = resolved_de.sort_values(
        ["human_gene_id", "finite_stat", "padj_sort", "abs_stat", "source_row"],
        ascending=[True, False, True, False, True],
    )
    resolved_de = resolved_de.drop_duplicates("human_gene_id", keep="first")
    universe = resolved_de[
        resolved_de["finite_stat"] & np.isfinite(resolved_de["baseMean"]) & (resolved_de["baseMean"] > 0)
    ].copy()
    expression = pd.read_excel(
        expr_path, sheet_name="S5 Normalized_data D_vs_C", engine="openpyxl"
    )
    expression = expression.rename(columns={expression.columns[0]: "source_gene"})
    expression["source_row"] = np.arange(2, len(expression) + 2)
    expression = add_resolution(expression, "source_gene", lookups)
    donors = [column for column in expression.columns if column.startswith("DPN") or column.startswith("Con")]
    for donor in donors:
        expression[donor] = pd.to_numeric(expression[donor], errors="coerce")
    expression["nonmissing_donor_count"] = expression[donors].notna().sum(axis=1)
    resolved_expr = expression[expression["human_gene_id"].notna()].copy()
    resolved_expr = resolved_expr.sort_values(
        ["human_gene_id", "nonmissing_donor_count", "source_row"],
        ascending=[True, False, True],
    ).drop_duplicates("human_gene_id", keep="first")
    expression_matrix = resolved_expr.set_index("human_gene_id")[donors]
    common = universe["human_gene_id"].isin(expression_matrix.index)
    universe = universe[common].copy()
    expression_matrix = expression_matrix.loc[universe["human_gene_id"]]
    expression_matrix.index.name = "human_gene_id"
    audit = {
        "all_gene_rows_deposited": int(len(de)),
        "de_rows_resolved_to_human_gene_id": int(de["human_gene_id"].notna().sum()),
        "de_unique_resolved_gene_ids": int(len(resolved_de)),
        "finite_positive_basemean_universe_with_expression": int(len(universe)),
        "excel_date_cells_recovered_de": int(de["gene_recovered_from_excel_date"].sum()),
        "excel_date_cells_recovered_expression": int(
            expression["gene_recovered_from_excel_date"].sum()
        ),
        "donor_columns": donors,
        "de_workbook_sha256": sha256(de_path),
        "expression_workbook_sha256": sha256(expr_path),
    }
    return universe, expression_matrix, audit


def load_metadata(donor_columns: list[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    path = RAW / "41598_2022_8100_MOESM2_ESM.xlsx"
    metadata = pd.read_excel(path, sheet_name="S1 Demographics", engine="openpyxl")
    metadata = metadata.rename(columns={"Donor": "donor"})
    metadata["group"] = metadata["Condition"].map({"Diabetic": "DPN", "Healthy": "Control"})
    metadata["Age"] = pd.to_numeric(metadata["Age"], errors="coerce")
    metadata["Sex"] = metadata["Sex"].astype(str).str.upper().str.strip()
    metadata = metadata[metadata["donor"].isin(donor_columns)].copy()
    metadata = metadata.set_index("donor").loc[donor_columns].reset_index()
    if set(metadata["donor"]) != set(donor_columns) or metadata["group"].isna().any():
        raise RuntimeError("Donor metadata does not map exactly to the normalized-expression columns")
    counts = metadata.groupby(["Sex", "group"]).size().to_dict()
    audit = {
        "metadata_workbook_sha256": sha256(path),
        "donor_count": int(len(metadata)),
        "group_counts": metadata["group"].value_counts().to_dict(),
        "sex_by_group_counts": {f"{sex}_{group}": int(count) for (sex, group), count in counts.items()},
        "age_by_group": metadata.groupby("group")["Age"].agg(["count", "mean", "min", "max"]).to_dict("index"),
    }
    return metadata, audit


def load_mapped_signatures(
    universe: pd.DataFrame,
    lookups: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    signatures = pd.read_csv(SIGNATURES, sep="\t")
    signatures = signatures[signatures["contrast_id"].isin(PRIMARY | CONTEXTUAL)].copy()
    signatures = add_resolution(signatures, "gene", lookups)
    target_columns = [
        "human_gene_id",
        "current_human_symbol",
        "baseMean",
        "log2FoldChange",
        "stat",
        "padj",
        "source_gene",
    ]
    target = universe[target_columns].rename(
        columns={
            "current_human_symbol": "target_current_symbol",
            "log2FoldChange": "target_log2FoldChange",
            "stat": "target_stat",
            "padj": "target_padj",
            "source_gene": "target_source_symbol",
        }
    )
    mapped = signatures.merge(target, on="human_gene_id", how="left")
    mapped["mapped_to_target_universe"] = mapped["target_stat"].notna()
    mapped = mapped.sort_values(["contrast_id", "direction", "gene"])
    unique_mapped = (
        mapped[mapped["mapped_to_target_universe"]]
        .sort_values(["contrast_id", "direction", "p_val_adj", "gene"])
        .drop_duplicates(["contrast_id", "direction", "human_gene_id"], keep="first")
    )
    return mapped, unique_mapped


def expression_matched_null(
    observed_up: pd.DataFrame,
    observed_down: pd.DataFrame,
    universe: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    stat_by_decile = {
        int(decile): frame["stat"].to_numpy(dtype=float)
        for decile, frame in universe.groupby("expression_decile")
    }
    up_counts = observed_up["expression_decile"].value_counts().to_dict()
    down_counts = observed_down["expression_decile"].value_counts().to_dict()
    null = np.empty(N_NULL, dtype=float)
    for iteration in range(N_NULL):
        up_sum = 0.0
        down_sum = 0.0
        for decile, count in up_counts.items():
            candidates = stat_by_decile[int(decile)]
            up_sum += float(rng.choice(candidates, size=int(count), replace=False).sum())
        for decile, count in down_counts.items():
            candidates = stat_by_decile[int(decile)]
            down_sum += float(rng.choice(candidates, size=int(count), replace=False).sum())
        null[iteration] = up_sum / len(observed_up) - down_sum / len(observed_down)
    observed = float(observed_up["stat"].mean() - observed_down["stat"].mean())
    p = (1 + int(np.sum(null >= observed - 1e-12))) / (N_NULL + 1)
    return null, float(p)


def hedges_g(scores: pd.Series, groups: pd.Series) -> float:
    x = scores[groups == "DPN"].to_numpy(dtype=float)
    y = scores[groups == "Control"].to_numpy(dtype=float)
    df = len(x) + len(y) - 2
    pooled = ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df
    if pooled <= 0:
        return float("nan")
    correction = 1 - 3 / (4 * df - 1)
    return float(((x.mean() - y.mean()) / math.sqrt(pooled)) * correction)


def leave_one_out(scores: pd.Series, groups: pd.Series) -> tuple[bool, float, float]:
    differences = []
    for donor in scores.index:
        kept = scores.index != donor
        current_scores = scores[kept]
        current_groups = groups[kept]
        differences.append(
            float(current_scores[current_groups == "DPN"].mean() - current_scores[current_groups == "Control"].mean())
        )
    return bool(all(value > 0 for value in differences)), min(differences), max(differences)


def sex_stratified_exact_permutation(
    scores: pd.Series, metadata: pd.DataFrame
) -> tuple[float, int, np.ndarray]:
    metadata = metadata.set_index("donor").loc[scores.index]
    observed = float(scores[metadata["group"] == "DPN"].mean() - scores[metadata["group"] == "Control"].mean())
    sex_options: list[list[tuple[str, ...]]] = []
    for sex in sorted(metadata["Sex"].unique()):
        donors = metadata.index[metadata["Sex"] == sex].tolist()
        n_dpn = int(((metadata["Sex"] == sex) & (metadata["group"] == "DPN")).sum())
        sex_options.append(list(itertools.combinations(donors, n_dpn)))
    null = []
    all_donors = set(metadata.index)
    for chosen_by_sex in itertools.product(*sex_options):
        dpn = set().union(*[set(chosen) for chosen in chosen_by_sex])
        control = all_donors - dpn
        null.append(float(scores[list(dpn)].mean() - scores[list(control)].mean()))
    null_array = np.asarray(null, dtype=float)
    p = float(np.mean(null_array >= observed - 1e-12))
    return p, len(null_array), null_array


def ols_condition_effect(scores: pd.Series, metadata: pd.DataFrame, include_age: bool) -> dict[str, float]:
    meta = metadata.set_index("donor").loc[scores.index]
    columns = [np.ones(len(meta)), (meta["group"] == "DPN").astype(float).to_numpy(), (meta["Sex"] == "M").astype(float).to_numpy()]
    if include_age:
        age = meta["Age"].to_numpy(dtype=float)
        columns.append((age - age.mean()) / age.std(ddof=1))
    design = np.column_stack(columns)
    outcome = scores.to_numpy(dtype=float)
    beta, _, rank, _ = np.linalg.lstsq(design, outcome, rcond=None)
    residual = outcome - design @ beta
    df = len(outcome) - rank
    sigma2 = float(residual @ residual / df)
    covariance = sigma2 * np.linalg.inv(design.T @ design)
    standard_error = float(math.sqrt(covariance[1, 1]))
    t_value = float(beta[1] / standard_error)
    p_two_sided = float(2 * stats.t.sf(abs(t_value), df))
    return {
        "coefficient": float(beta[1]),
        "standard_error": standard_error,
        "t": t_value,
        "df": int(df),
        "two_sided_p": p_two_sided,
    }


def plot_results(scores_long: pd.DataFrame, tests: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    primary_order = ["late_neuron_DPN_vs_diabetes", "severity_neuron_modhigh_vs_low_nageotte"]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.1), gridspec_kw={"width_ratios": [1, 1, 1.35]})
    rng = np.random.default_rng(RANDOM_SEED)
    colors = {"Control": "#4C78A8", "DPN": "#D1495B"}
    markers = {"F": "o", "M": "^"}
    for axis, contrast in zip(axes[:2], primary_order):
        subset = scores_long[scores_long["contrast_id"] == contrast].copy()
        test = tests.set_index("contrast_id").loc[contrast]
        for x, group in enumerate(["Control", "DPN"]):
            group_data = subset[subset["group"] == group]
            for row in group_data.itertuples(index=False):
                jitter = float(rng.uniform(-0.08, 0.08))
                axis.scatter(
                    x + jitter,
                    row.score,
                    s=58,
                    marker=markers[row.Sex],
                    color=colors[group],
                    edgecolor="white",
                    linewidth=0.7,
                    zorder=3,
                )
            mean = group_data["score"].mean()
            sem = group_data["score"].sem()
            axis.errorbar(x, mean, yerr=sem, color="black", capsize=4, marker="_", markersize=20, lw=1.3, zorder=4)
        axis.axhline(0, color="#B8B8B8", lw=0.8, zorder=0)
        axis.set_xticks([0, 1], ["Control\n(n=7)", "DPN\n(n=5)"])
        axis.set_ylabel("Within-donor rank signature score")
        axis.set_title(DISPLAY[contrast], fontsize=11, fontweight="bold")
        axis.text(
            0.03,
            0.97,
            f"Δ={test['donor_score_difference']:.3f}\nsex-stratified Q={test['donor_exact_bh_q']:.3g}",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
        )
        axis.spines[["top", "right"]].set_visible(False)
    ordered = [
        "early_allcell_diabetes_vs_control",
        "late_allcell_DPN_vs_diabetes",
        "late_neuron_DPN_vs_diabetes",
        "severity_neuron_modhigh_vs_low_nageotte",
        "xenium_DPN_vs_control",
    ]
    test_index = tests.set_index("contrast_id").loc[ordered]
    z = (test_index["gene_concordance"] - test_index["matched_null_mean"]) / test_index["matched_null_sd"]
    bar_colors = ["#72A0C1" if contrast not in PRIMARY else "#2A9D8F" for contrast in ordered]
    axes[2].barh(np.arange(len(ordered)), z, color=bar_colors)
    axes[2].axvline(0, color="#555555", lw=0.8)
    axes[2].set_yticks(np.arange(len(ordered)), [DISPLAY[item] for item in ordered])
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Gene-level concordance vs matched null (SD)")
    axes[2].set_title("Cross-cohort directional concordance", fontsize=11, fontweight="bold")
    for y, contrast in enumerate(ordered):
        row = test_index.loc[contrast]
        q = row["gene_matched_bh_q"]
        axes[2].text(z.loc[contrast] + (0.12 if z.loc[contrast] >= 0 else -0.12), y, f"Q={q:.3g}", va="center", ha="left" if z.loc[contrast] >= 0 else "right", fontsize=8)
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.suptitle("Independent human hDRG bulk validation: painful DPN vs non-diabetic control", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.01, "Circles: female donors; triangles: male donors. Error bars show mean ± SEM.", ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.94])
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"human_DPN_bulk_signature_validation_{DATE}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    gene_info_path = NCBI / "Homo_sapiens.gene_info.gz"
    lookups = build_symbol_lookup(load_gene_info(gene_info_path))
    universe, expression, target_audit = load_target(lookups)
    metadata, metadata_audit = load_metadata(expression.columns.tolist())
    mapped_all, mapped_unique = load_mapped_signatures(universe, lookups)

    universe = universe.copy()
    universe["expression_decile"] = pd.qcut(
        universe["baseMean"].rank(method="first"), q=10, labels=False
    ).astype(int)
    universe_index = universe.set_index("human_gene_id")
    centered_ranks = expression.rank(axis=0, method="average", pct=True) - 0.5
    groups = metadata.set_index("donor")["group"].loc[centered_ranks.columns]

    rng = np.random.default_rng(RANDOM_SEED)
    tests: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    null_summaries: dict[str, object] = {}
    for contrast in sorted(PRIMARY | CONTEXTUAL):
        signature = mapped_unique[mapped_unique["contrast_id"] == contrast]
        up_ids = signature.loc[signature["direction"] == "up", "human_gene_id"].tolist()
        down_ids = signature.loc[signature["direction"] == "down", "human_gene_id"].tolist()
        minimum = 10 if contrast in PRIMARY or contrast in {"early_allcell_diabetes_vs_control", "late_allcell_DPN_vs_diabetes"} else 5
        if len(up_ids) < minimum or len(down_ids) < minimum:
            raise RuntimeError(f"Insufficient mapping for {contrast}: up={len(up_ids)}, down={len(down_ids)}, minimum={minimum}")
        observed_up = universe_index.loc[up_ids].reset_index()
        observed_down = universe_index.loc[down_ids].reset_index()
        null, matched_p = expression_matched_null(observed_up, observed_down, universe, rng)
        gene_concordance = float(observed_up["stat"].mean() - observed_down["stat"].mean())
        donor_scores = centered_ranks.loc[up_ids].mean(axis=0) - centered_ranks.loc[down_ids].mean(axis=0)
        donor_difference = float(donor_scores[groups == "DPN"].mean() - donor_scores[groups == "Control"].mean())
        donor_p, permutation_count, donor_null = sex_stratified_exact_permutation(donor_scores, metadata)
        loo_all_positive, loo_min, loo_max = leave_one_out(donor_scores, groups)
        model_sex = ols_condition_effect(donor_scores, metadata, include_age=False)
        model_sex_age = ols_condition_effect(donor_scores, metadata, include_age=True)
        tests.append(
            {
                "contrast_id": contrast,
                "test_family": "primary_neuronal" if contrast in PRIMARY else "contextual_secondary",
                "mapped_up_n": len(up_ids),
                "mapped_down_n": len(down_ids),
                "gene_concordance": gene_concordance,
                "up_mean_target_stat": float(observed_up["stat"].mean()),
                "down_mean_target_stat": float(observed_down["stat"].mean()),
                "up_median_target_log2fc": float(observed_up["log2FoldChange"].median()),
                "down_median_target_log2fc": float(observed_down["log2FoldChange"].median()),
                "up_fraction_directionally_concordant": float((observed_up["log2FoldChange"] > 0).mean()),
                "down_fraction_directionally_concordant": float((observed_down["log2FoldChange"] < 0).mean()),
                "matched_null_mean": float(null.mean()),
                "matched_null_sd": float(null.std(ddof=1)),
                "gene_matched_one_sided_p": matched_p,
                "donor_control_mean": float(donor_scores[groups == "Control"].mean()),
                "donor_dpn_mean": float(donor_scores[groups == "DPN"].mean()),
                "donor_score_difference": donor_difference,
                "donor_hedges_g": hedges_g(donor_scores, groups),
                "donor_sex_stratified_exact_one_sided_p": donor_p,
                "donor_exact_permutation_count": permutation_count,
                "loo_all_positive": loo_all_positive,
                "loo_min_difference": loo_min,
                "loo_max_difference": loo_max,
                "ols_dpn_coefficient_sex_adjusted": model_sex["coefficient"],
                "ols_dpn_p_sex_adjusted": model_sex["two_sided_p"],
                "ols_dpn_coefficient_sex_age_adjusted": model_sex_age["coefficient"],
                "ols_dpn_p_sex_age_adjusted": model_sex_age["two_sided_p"],
            }
        )
        null_summaries[contrast] = {
            "gene_matched_null_quantiles": {str(q): float(np.quantile(null, q)) for q in [0, 0.025, 0.5, 0.975, 1]},
            "donor_sex_stratified_null_quantiles": {str(q): float(np.quantile(donor_null, q)) for q in [0, 0.025, 0.5, 0.975, 1]},
        }
        meta_index = metadata.set_index("donor")
        for donor, score in donor_scores.items():
            score_rows.append(
                {
                    "contrast_id": contrast,
                    "donor": donor,
                    "group": meta_index.loc[donor, "group"],
                    "Sex": meta_index.loc[donor, "Sex"],
                    "Age": meta_index.loc[donor, "Age"],
                    "score": float(score),
                }
            )

    tests_df = pd.DataFrame(tests)
    tests_df["gene_matched_bh_q"] = np.nan
    tests_df["donor_exact_bh_q"] = np.nan
    for family in ["primary_neuronal", "contextual_secondary"]:
        mask = tests_df["test_family"] == family
        tests_df.loc[mask, "gene_matched_bh_q"] = bh_adjust(
            tests_df.loc[mask, "gene_matched_one_sided_p"]
        )
        tests_df.loc[mask, "donor_exact_bh_q"] = bh_adjust(
            tests_df.loc[mask, "donor_sex_stratified_exact_one_sided_p"]
        )
    tests_df["frozen_gate_supportive"] = (
        tests_df["contrast_id"].isin(PRIMARY)
        & (tests_df["mapped_up_n"] >= 10)
        & (tests_df["mapped_down_n"] >= 10)
        & (tests_df["gene_concordance"] > 0)
        & (tests_df["gene_matched_bh_q"] < 0.10)
        & (tests_df["donor_score_difference"] > 0)
        & (tests_df["donor_exact_bh_q"] < 0.10)
        & tests_df["loo_all_positive"]
    )
    terminal_gate = bool(tests_df.loc[tests_df["contrast_id"].isin(PRIMARY), "frozen_gate_supportive"].any())

    scores_df = pd.DataFrame(score_rows)
    mapped_all.to_csv(TABLES / f"hDRG_signatures_to_independent_human_bulk_mapping_{DATE}.tsv", sep="\t", index=False)
    tests_df.sort_values(["test_family", "contrast_id"]).to_csv(
        TABLES / f"independent_human_DPN_bulk_signature_tests_{DATE}.tsv", sep="\t", index=False
    )
    scores_df.to_csv(TABLES / f"independent_human_DPN_bulk_donor_scores_{DATE}.tsv", sep="\t", index=False)
    metadata.to_csv(METADATA / f"independent_human_DPN_bulk_donor_metadata_{DATE}.tsv", sep="\t", index=False)
    universe[["human_gene_id", "current_human_symbol", "source_gene", "baseMean", "log2FoldChange", "stat", "padj", "expression_decile"]].to_csv(
        TABLES / f"independent_human_DPN_bulk_target_universe_{DATE}.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    plot_results(scores_df, tests_df)

    qc = {
        "status": "PASS",
        "frozen_spec": "INDEPENDENT_HUMAN_BULK_VALIDATION_FROZEN_SPEC_2026-08-27.md",
        "source_article": {
            "title": "Transcriptomic analysis of human sensory neurons in painful diabetic neuropathy reveals inflammation and neuronal loss",
            "pmid": "35304484",
            "pmcid": "PMC8933403",
            "doi": "10.1038/s41598-022-08100-8",
        },
        "target_audit": target_audit,
        "metadata_audit": metadata_audit,
        "signature_source_sha256": sha256(SIGNATURES),
        "ncbi_gene_info_sha256": sha256(gene_info_path),
        "random_seed": RANDOM_SEED,
        "matched_null_iterations": N_NULL,
        "null_summaries": null_summaries,
        "terminal_human_dpn_gate_pass": terminal_gate,
        "gate_interpretation": (
            "A pass supports a reproducible terminal human DPN program but does not validate the specific diabetes-to-DPN transition."
        ),
    }
    with (TABLES / f"independent_human_DPN_bulk_validation_qc_{DATE}.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)

    print(tests_df.sort_values(["test_family", "contrast_id"]).to_string(index=False))
    print(json.dumps({"terminal_human_dpn_gate_pass": terminal_gate, **target_audit, **metadata_audit}, indent=2))


if __name__ == "__main__":
    main()
