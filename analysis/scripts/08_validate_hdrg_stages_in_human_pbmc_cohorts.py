#!/usr/bin/env python3
"""Validate frozen human hDRG stage signatures in two human PBMC cohorts."""

from __future__ import annotations

import gzip
import hashlib
import io
import itertools
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PHASE = Path(__file__).resolve().parents[2]
RAW = PHASE / "data" / "raw" / "human_PBMC_stage_cohorts"
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
METADATA = PHASE / "metadata"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
SIGNATURES = TABLES / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"
DATE = "2026-08-27"
SEED = 20260827
N_NULL = 10_000

EARLY = "early_allcell_diabetes_vs_control"
LATE = "late_allcell_DPN_vs_diabetes"
LATE_NEURON = "late_neuron_DPN_vs_diabetes"
SEVERITY = "severity_neuron_modhigh_vs_low_nageotte"
XENIUM = "xenium_DPN_vs_control"


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


def build_lookup(info: pd.DataFrame) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    exact: dict[str, str] = {}
    folded_candidates: defaultdict[str, set[str]] = defaultdict(set)
    synonym_candidates: defaultdict[str, set[str]] = defaultdict(set)
    id_to_symbol: dict[str, str] = {}
    for row in info.itertuples(index=False):
        gene_id = str(row.GeneID)
        symbol = str(row.Symbol)
        exact[symbol] = gene_id
        folded_candidates[symbol.upper()].add(gene_id)
        id_to_symbol[gene_id] = symbol
        if row.Synonyms and row.Synonyms != "-":
            for synonym in str(row.Synonyms).split("|"):
                if synonym and synonym != "-":
                    synonym_candidates[synonym.upper()].add(gene_id)
    folded = {key: next(iter(ids)) for key, ids in folded_candidates.items() if len(ids) == 1}
    synonyms = {
        key: next(iter(ids))
        for key, ids in synonym_candidates.items()
        if len(ids) == 1 and key not in folded
    }
    return exact, folded, synonyms, id_to_symbol


def resolve_single(symbol: object, lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[str | None, str]:
    text = str(symbol).strip()
    exact, folded, synonyms, _ = lookup
    if text in exact:
        return exact[text], "official_exact"
    key = text.upper()
    if key in folded:
        return folded[key], "official_casefold"
    if key in synonyms:
        return synonyms[key], "unique_synonym"
    return None, "unresolved_or_ambiguous"


def resolve_annotation(value: object, lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[str | None, str]:
    text = str(value).strip()
    candidates = [item.strip() for item in re.split(r"\s*(?:///|[;,])\s*", text) if item.strip()]
    if not candidates:
        return None, "empty"
    resolved = [resolve_single(candidate, lookup) for candidate in candidates]
    gene_ids = {item[0] for item in resolved if item[0] is not None}
    if len(gene_ids) == 1:
        methods = {item[1] for item in resolved if item[0] is not None}
        return next(iter(gene_ids)), "+".join(sorted(methods))
    return None, "multi_or_unresolved_annotation"


def read_geo_table(path: Path, begin_marker: str, end_marker: str) -> pd.DataFrame:
    rows = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped == begin_marker:
                in_table = True
                continue
            if stripped == end_marker:
                break
            if in_table:
                rows.append(stripped)
    if not rows:
        raise RuntimeError(f"No table found in {path}")
    return pd.read_csv(io.StringIO("\n".join(rows)), sep="\t")


def read_geo_sample_metadata(path: Path) -> pd.DataFrame:
    selected: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession") or line.startswith("!Sample_title") or line.startswith("!Sample_characteristics_ch1"):
                parts = [part.strip().strip('"') for part in line.rstrip("\n").split("\t")]
                selected.setdefault(parts[0], parts[1:])
    accessions = selected["!Sample_geo_accession"]
    titles = selected["!Sample_title"]
    frame = pd.DataFrame({"sample": accessions, "title": titles})
    frame["group"] = frame["title"].str.extract(r"^(DM|CN|DPN)", expand=False).map({"DM": "DM", "CN": "HC", "DPN": "DPN"})
    if frame["group"].isna().any():
        raise RuntimeError("Could not derive all GSE95849 groups")
    return frame


def prepare_gse95849(lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    matrix_path = RAW / "GSE95849" / "GSE95849_series_matrix.txt.gz"
    platform_path = RAW / "GSE95849" / "GPL22448_family.soft.gz"
    matrix = read_geo_table(matrix_path, "!series_matrix_table_begin", "!series_matrix_table_end")
    matrix = matrix.rename(columns={matrix.columns[0]: "probe_id"}).set_index("probe_id")
    metadata = read_geo_sample_metadata(matrix_path)
    matrix = matrix[metadata["sample"].tolist()].apply(pd.to_numeric, errors="coerce")
    platform = read_geo_table(platform_path, "!platform_table_begin", "!platform_table_end")
    platform = platform[["ID", "Gene_symbol", "Transcript_type"]].rename(columns={"ID": "probe_id"})
    resolved = platform["Gene_symbol"].map(lambda value: resolve_annotation(value, lookup))
    platform["human_gene_id"] = [item[0] for item in resolved]
    platform["mapping_method"] = [item[1] for item in resolved]
    platform = platform[platform["human_gene_id"].notna() & platform["probe_id"].isin(matrix.index)].copy()
    log_matrix = np.log2(matrix.clip(lower=0) + 1)
    platform["across_sample_median"] = platform["probe_id"].map(log_matrix.median(axis=1))
    platform = platform.sort_values(
        ["human_gene_id", "across_sample_median", "probe_id"], ascending=[True, False, True]
    ).drop_duplicates("human_gene_id", keep="first")
    expression = log_matrix.loc[platform["probe_id"]].copy()
    expression.index = platform["human_gene_id"].to_numpy()
    expression.index.name = "human_gene_id"
    expression = expression[np.isfinite(expression).all(axis=1)]
    audit = {
        "deposited_probe_rows": int(len(matrix)),
        "annotated_resolved_probe_rows": int(platform["human_gene_id"].notna().sum()),
        "unique_gene_ids_after_probe_reduction": int(len(expression)),
        "sample_count": int(expression.shape[1]),
        "group_counts": metadata["group"].value_counts().to_dict(),
        "series_matrix_sha256": sha256(matrix_path),
        "platform_soft_sha256": sha256(platform_path),
    }
    return expression, metadata, audit


def prepare_gse185011(lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    path = RAW / "GSE185011" / "GSE185011_AllSample_matrix_FPKM.txt.gz"
    raw = pd.read_csv(path, sep="\t", compression="gzip")
    raw = raw.rename(columns={raw.columns[0]: "source_gene"})
    raw["source_row"] = np.arange(2, len(raw) + 2)
    resolved = raw["source_gene"].map(lambda value: resolve_single(value, lookup))
    raw["human_gene_id"] = [item[0] for item in resolved]
    raw["mapping_method"] = [item[1] for item in resolved]
    samples = [column for column in raw.columns if re.match(r"^(HC|T2DM|DR|DPN|DN)\d+$", str(column))]
    values = raw[samples].apply(pd.to_numeric, errors="coerce")
    raw["across_sample_mean_raw"] = values.mean(axis=1)
    raw["across_sample_median_log"] = np.log2(values.clip(lower=0) + 1).median(axis=1)
    resolved_raw = raw[raw["human_gene_id"].notna() & (raw["across_sample_mean_raw"] > 0)].sort_values(
        ["human_gene_id", "across_sample_median_log", "source_row"], ascending=[True, False, True]
    ).drop_duplicates("human_gene_id", keep="first")
    expression = np.log2(values.loc[resolved_raw.index].clip(lower=0) + 1)
    expression.index = resolved_raw["human_gene_id"].to_numpy()
    expression.index.name = "human_gene_id"
    expression = expression[np.isfinite(expression).all(axis=1)]
    metadata = pd.DataFrame({"sample": samples})
    metadata["group"] = metadata["sample"].str.extract(r"^(HC|T2DM|DR|DPN|DN)", expand=False).replace({"T2DM": "DM"})
    audit = {
        "deposited_gene_rows": int(len(raw)),
        "resolved_positive_rows": int(len(resolved_raw)),
        "unique_gene_ids": int(len(expression)),
        "sample_count": int(expression.shape[1]),
        "group_counts": metadata["group"].value_counts().to_dict(),
        "fpkm_matrix_sha256": sha256(path),
    }
    return expression, metadata, audit


def load_signatures(lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> pd.DataFrame:
    signatures = pd.read_csv(SIGNATURES, sep="\t")
    signatures = signatures[signatures["contrast_id"].isin([EARLY, LATE, LATE_NEURON, SEVERITY, XENIUM])].copy()
    resolved = signatures["gene"].map(lambda value: resolve_single(value, lookup))
    signatures["human_gene_id"] = [item[0] for item in resolved]
    signatures["mapping_method"] = [item[1] for item in resolved]
    signatures = signatures[signatures["human_gene_id"].notna()].sort_values(
        ["contrast_id", "direction", "p_val_adj", "gene"]
    ).drop_duplicates(["contrast_id", "direction", "human_gene_id"], keep="first")
    return signatures


def expression_matched_null(up: pd.DataFrame, down: pd.DataFrame, universe: pd.DataFrame, rng: np.random.Generator) -> tuple[np.ndarray, float]:
    effects = {int(decile): frame["effect"].to_numpy(dtype=float) for decile, frame in universe.groupby("expression_decile")}
    up_counts = up["expression_decile"].value_counts().to_dict()
    down_counts = down["expression_decile"].value_counts().to_dict()
    null = np.empty(N_NULL, dtype=float)
    for iteration in range(N_NULL):
        up_sum = sum(float(rng.choice(effects[int(decile)], int(count), replace=False).sum()) for decile, count in up_counts.items())
        down_sum = sum(float(rng.choice(effects[int(decile)], int(count), replace=False).sum()) for decile, count in down_counts.items())
        null[iteration] = up_sum / len(up) - down_sum / len(down)
    observed = float(up["effect"].mean() - down["effect"].mean())
    p = (1 + int(np.sum(null >= observed - 1e-12))) / (N_NULL + 1)
    return null, float(p)


def exact_permutation(scores: pd.Series, groups: pd.Series, earlier: str, later: str) -> tuple[float, int]:
    keep = groups.isin([earlier, later])
    current_scores = scores[keep]
    current_groups = groups[keep]
    n_later = int((current_groups == later).sum())
    observed = float(current_scores[current_groups == later].mean() - current_scores[current_groups == earlier].mean())
    values = current_scores.to_numpy(dtype=float)
    null = []
    for later_indices in itertools.combinations(range(len(values)), n_later):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(later_indices)] = True
        null.append(float(values[mask].mean() - values[~mask].mean()))
    return float(np.mean(np.asarray(null) >= observed - 1e-12)), len(null)


def hedges_g(scores: pd.Series, groups: pd.Series, earlier: str, later: str) -> tuple[float, float]:
    x = scores[groups == later].to_numpy(dtype=float)
    y = scores[groups == earlier].to_numpy(dtype=float)
    df = len(x) + len(y) - 2
    pooled = ((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / df
    if pooled <= 0:
        return float("nan"), float("nan")
    correction = 1 - 3 / (4 * df - 1)
    g = float(((x.mean() - y.mean()) / math.sqrt(pooled)) * correction)
    variance = float((len(x) + len(y)) / (len(x) * len(y)) + (g * g) / (2 * df))
    return g, variance


def leave_one_out(scores: pd.Series, groups: pd.Series, earlier: str, later: str) -> tuple[bool, float, float]:
    keep = groups.isin([earlier, later])
    scores = scores[keep]
    groups = groups[keep]
    differences = []
    for sample in scores.index:
        retained = scores.index != sample
        current_scores = scores[retained]
        current_groups = groups[retained]
        differences.append(float(current_scores[current_groups == later].mean() - current_scores[current_groups == earlier].mean()))
    return bool(all(item > 0 for item in differences)), min(differences), max(differences)


def fixed_effect_meta(rows: pd.DataFrame, transition: str) -> dict[str, object]:
    current = rows[rows["transition"] == transition].copy()
    weights = 1 / current["hedges_g_variance"].to_numpy(dtype=float)
    effects = current["hedges_g"].to_numpy(dtype=float)
    pooled = float(np.sum(weights * effects) / np.sum(weights))
    se = float(math.sqrt(1 / np.sum(weights)))
    z = pooled / se
    p = float(2 * stats.norm.sf(abs(z)))
    q_stat = float(np.sum(weights * (effects - pooled) ** 2))
    df = len(effects) - 1
    i2 = float(max(0, (q_stat - df) / q_stat) * 100) if q_stat > 0 else 0.0
    return {
        "transition": transition,
        "cohort_count": int(len(current)),
        "fixed_effect_hedges_g": pooled,
        "standard_error": se,
        "z": z,
        "two_sided_p": p,
        "cochran_Q": q_stat,
        "heterogeneity_df": df,
        "I2_percent": i2,
        "all_cohort_score_differences_positive": bool((current["score_difference"] > 0).all()),
        "at_least_one_full_support": bool(current["projection_supportive"].any()),
    }


def plot_stage(scores: pd.DataFrame, tests: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), sharey=False)
    rng = np.random.default_rng(SEED)
    colors = {"HC": "#4C78A8", "DM": "#F2A541", "DPN": "#D1495B"}
    for row, cohort in enumerate(["GSE95849", "GSE185011"]):
        for column, signature in enumerate([EARLY, LATE]):
            axis = axes[row, column]
            subset = scores[(scores["cohort"] == cohort) & (scores["contrast_id"] == signature) & scores["group"].isin(["HC", "DM", "DPN"])]
            for x, group in enumerate(["HC", "DM", "DPN"]):
                values = subset.loc[subset["group"] == group, "score"].to_numpy(dtype=float)
                jitter = rng.uniform(-0.08, 0.08, size=len(values))
                axis.scatter(np.full(len(values), x) + jitter, values, color=colors[group], s=48, edgecolor="white", linewidth=0.6)
                axis.errorbar(x, values.mean(), yerr=values.std(ddof=1) / math.sqrt(len(values)), color="black", marker="_", markersize=18, capsize=4)
            axis.axhline(0, color="#BBBBBB", lw=0.7)
            axis.set_xticks([0, 1, 2], ["Healthy", "T2DM", "DPN"])
            axis.set_ylabel("Within-sample rank score")
            axis.set_title(f"{cohort}: {'early all-cell' if signature == EARLY else 'late all-cell'}", fontweight="bold")
            transition = "early_HC_to_DM" if signature == EARLY else "late_DM_to_DPN"
            result = tests[(tests["cohort"] == cohort) & (tests["transition"] == transition)].iloc[0]
            axis.text(0.03, 0.97, f"g={result['hedges_g']:.2f}; exact Q={result['exact_bh_q']:.3g}\ngene Q={result['gene_matched_bh_q']:.3g}", transform=axis.transAxes, va="top", fontsize=9)
            axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Human hDRG disease-stage signatures in two PBMC cohorts", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"human_PBMC_stage_signature_validation_{DATE}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.2, 5.2))
    subset = scores[(scores["cohort"] == "GSE185011") & (scores["contrast_id"] == LATE)]
    order = ["DM", "DPN", "DR", "DN"]
    colors_specific = {"DM": "#F2A541", "DPN": "#D1495B", "DR": "#7A5195", "DN": "#5F9E6E"}
    for x, group in enumerate(order):
        values = subset.loc[subset["group"] == group, "score"].to_numpy(dtype=float)
        jitter = rng.uniform(-0.08, 0.08, size=len(values))
        axis.scatter(np.full(len(values), x) + jitter, values, color=colors_specific[group], s=55, edgecolor="white", linewidth=0.6)
        axis.errorbar(x, values.mean(), yerr=values.std(ddof=1) / math.sqrt(len(values)), color="black", marker="_", markersize=18, capsize=4)
    axis.set_xticks(range(len(order)), ["T2DM", "DPN", "Retinopathy", "Nephropathy"])
    axis.set_ylabel("Late all-cell hDRG rank score")
    axis.set_title("GSE185011 complication-specificity boundary", fontweight="bold")
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"human_PBMC_DPN_complication_specificity_{DATE}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    info_path = NCBI / "Homo_sapiens.gene_info.gz"
    lookup = build_lookup(load_gene_info(info_path))
    signatures = load_signatures(lookup)
    cohorts = {
        "GSE95849": prepare_gse95849(lookup),
        "GSE185011": prepare_gse185011(lookup),
    }
    rng = np.random.default_rng(SEED)
    test_rows = []
    score_rows = []
    mapping_rows = []
    audits = {}
    null_summaries = {}

    test_definitions = [
        (EARLY, "HC", "DM", "early_HC_to_DM", "primary_transition"),
        (LATE, "DM", "DPN", "late_DM_to_DPN", "primary_transition"),
        (LATE_NEURON, "DM", "DPN", "late_neuron_DM_to_DPN", "secondary"),
        (SEVERITY, "DM", "DPN", "severity_neuron_DM_to_DPN", "secondary"),
        (XENIUM, "HC", "DPN", "xenium_HC_to_DPN", "secondary"),
    ]

    for cohort, (expression, metadata, audit) in cohorts.items():
        audits[cohort] = audit
        group_by_sample = metadata.set_index("sample")["group"].loc[expression.columns]
        ranks = expression.rank(axis=0, method="average", pct=True) - 0.5
        mean_expression = expression.mean(axis=1)
        for signature, earlier, later, transition, family in test_definitions:
            source = signatures[signatures["contrast_id"] == signature].copy()
            source["cohort"] = cohort
            source["mapped_to_target"] = source["human_gene_id"].isin(expression.index)
            source["target_current_symbol"] = source["human_gene_id"].map(lookup[3])
            mapping_rows.append(source)
            mapped = source[source["mapped_to_target"]]
            up_ids = mapped.loc[mapped["direction"] == "up", "human_gene_id"].drop_duplicates().tolist()
            down_ids = mapped.loc[mapped["direction"] == "down", "human_gene_id"].drop_duplicates().tolist()
            minimum = 5 if signature == XENIUM else 10
            if len(up_ids) < minimum or len(down_ids) < minimum:
                raise RuntimeError(f"Insufficient mapping {cohort}/{signature}: up={len(up_ids)}, down={len(down_ids)}")
            samples = group_by_sample[group_by_sample.isin([earlier, later])].index
            target_effect = expression[samples[group_by_sample.loc[samples] == later]].mean(axis=1) - expression[samples[group_by_sample.loc[samples] == earlier]].mean(axis=1)
            universe = pd.DataFrame({"effect": target_effect, "mean_expression": mean_expression.loc[target_effect.index]})
            universe = universe[np.isfinite(universe).all(axis=1)].copy()
            universe["expression_decile"] = pd.qcut(universe["mean_expression"].rank(method="first"), 10, labels=False).astype(int)
            up = universe.loc[up_ids].reset_index()
            down = universe.loc[down_ids].reset_index()
            null, gene_p = expression_matched_null(up, down, universe, rng)
            gene_concordance = float(up["effect"].mean() - down["effect"].mean())
            scores = ranks.loc[up_ids].mean(axis=0) - ranks.loc[down_ids].mean(axis=0)
            difference = float(scores[group_by_sample == later].mean() - scores[group_by_sample == earlier].mean())
            exact_p, permutation_count = exact_permutation(scores, group_by_sample, earlier, later)
            g, g_variance = hedges_g(scores, group_by_sample, earlier, later)
            loo_positive, loo_min, loo_max = leave_one_out(scores, group_by_sample, earlier, later)
            test_rows.append(
                {
                    "cohort": cohort,
                    "contrast_id": signature,
                    "transition": transition,
                    "test_family": family,
                    "earlier_group": earlier,
                    "later_group": later,
                    "earlier_n": int((group_by_sample == earlier).sum()),
                    "later_n": int((group_by_sample == later).sum()),
                    "mapped_up_n": len(up_ids),
                    "mapped_down_n": len(down_ids),
                    "gene_concordance": gene_concordance,
                    "gene_matched_p": gene_p,
                    "matched_null_mean": float(null.mean()),
                    "matched_null_sd": float(null.std(ddof=1)),
                    "earlier_score_mean": float(scores[group_by_sample == earlier].mean()),
                    "later_score_mean": float(scores[group_by_sample == later].mean()),
                    "score_difference": difference,
                    "hedges_g": g,
                    "hedges_g_variance": g_variance,
                    "exact_one_sided_p": exact_p,
                    "exact_permutation_count": permutation_count,
                    "loo_all_positive": loo_positive,
                    "loo_min_difference": loo_min,
                    "loo_max_difference": loo_max,
                }
            )
            null_summaries[f"{cohort}__{transition}"] = {
                "gene_null_quantiles": {str(q): float(np.quantile(null, q)) for q in [0, 0.025, 0.5, 0.975, 1]}
            }
            for sample, score in scores.items():
                score_rows.append({"cohort": cohort, "contrast_id": signature, "sample": sample, "group": group_by_sample.loc[sample], "score": float(score)})

        if cohort == "GSE185011":
            signature = LATE
            source = signatures[signatures["contrast_id"] == signature]
            up_ids = source.loc[(source["direction"] == "up") & source["human_gene_id"].isin(expression.index), "human_gene_id"].drop_duplicates().tolist()
            down_ids = source.loc[(source["direction"] == "down") & source["human_gene_id"].isin(expression.index), "human_gene_id"].drop_duplicates().tolist()
            scores = ranks.loc[up_ids].mean(axis=0) - ranks.loc[down_ids].mean(axis=0)
            for other in ["DR", "DN"]:
                earlier, later = other, "DPN"
                target_effect = expression.loc[:, group_by_sample.isin([earlier, later])].groupby(group_by_sample[group_by_sample.isin([earlier, later])], axis=1).mean()
                effect = target_effect[later] - target_effect[earlier]
                universe = pd.DataFrame({"effect": effect, "mean_expression": mean_expression.loc[effect.index]})
                universe["expression_decile"] = pd.qcut(universe["mean_expression"].rank(method="first"), 10, labels=False).astype(int)
                up = universe.loc[up_ids].reset_index()
                down = universe.loc[down_ids].reset_index()
                null, gene_p = expression_matched_null(up, down, universe, rng)
                exact_p, permutation_count = exact_permutation(scores, group_by_sample, earlier, later)
                g, g_variance = hedges_g(scores, group_by_sample, earlier, later)
                loo_positive, loo_min, loo_max = leave_one_out(scores, group_by_sample, earlier, later)
                difference = float(scores[group_by_sample == later].mean() - scores[group_by_sample == earlier].mean())
                test_rows.append(
                    {
                        "cohort": cohort,
                        "contrast_id": signature,
                        "transition": f"DPN_vs_{other}",
                        "test_family": "specificity",
                        "earlier_group": earlier,
                        "later_group": later,
                        "earlier_n": 5,
                        "later_n": 5,
                        "mapped_up_n": len(up_ids),
                        "mapped_down_n": len(down_ids),
                        "gene_concordance": float(up["effect"].mean() - down["effect"].mean()),
                        "gene_matched_p": gene_p,
                        "matched_null_mean": float(null.mean()),
                        "matched_null_sd": float(null.std(ddof=1)),
                        "earlier_score_mean": float(scores[group_by_sample == earlier].mean()),
                        "later_score_mean": float(scores[group_by_sample == later].mean()),
                        "score_difference": difference,
                        "hedges_g": g,
                        "hedges_g_variance": g_variance,
                        "exact_one_sided_p": exact_p,
                        "exact_permutation_count": permutation_count,
                        "loo_all_positive": loo_positive,
                        "loo_min_difference": loo_min,
                        "loo_max_difference": loo_max,
                    }
                )
                null_summaries[f"{cohort}__DPN_vs_{other}"] = {
                    "gene_null_quantiles": {str(q): float(np.quantile(null, q)) for q in [0, 0.025, 0.5, 0.975, 1]}
                }

    tests = pd.DataFrame(test_rows)
    tests["gene_matched_bh_q"] = np.nan
    tests["exact_bh_q"] = np.nan
    for family in tests["test_family"].unique():
        mask = tests["test_family"] == family
        tests.loc[mask, "gene_matched_bh_q"] = bh_adjust(tests.loc[mask, "gene_matched_p"])
        tests.loc[mask, "exact_bh_q"] = bh_adjust(tests.loc[mask, "exact_one_sided_p"])
    tests["projection_supportive"] = (
        (tests["test_family"] == "primary_transition")
        & (tests["score_difference"] > 0)
        & (tests["hedges_g"] >= 0.5)
        & (tests["exact_bh_q"] < 0.10)
        & tests["loo_all_positive"]
        & (tests["gene_matched_bh_q"] < 0.10)
    )

    meta_rows = [fixed_effect_meta(tests[tests["test_family"] == "primary_transition"], transition) for transition in ["early_HC_to_DM", "late_DM_to_DPN"]]
    meta = pd.DataFrame(meta_rows)
    early_meta = meta.set_index("transition").loc["early_HC_to_DM"]
    late_meta = meta.set_index("transition").loc["late_DM_to_DPN"]
    specificity = tests[tests["test_family"] == "specificity"]
    specificity_pass = bool(
        (specificity["score_difference"] > 0).all()
        and (specificity["exact_one_sided_p"] < 0.10).any()
    )
    early_pass = bool(
        early_meta["all_cohort_score_differences_positive"]
        and early_meta["fixed_effect_hedges_g"] >= 0.5
        and early_meta["at_least_one_full_support"]
    )
    late_pass = bool(
        late_meta["all_cohort_score_differences_positive"]
        and late_meta["fixed_effect_hedges_g"] >= 0.5
        and late_meta["at_least_one_full_support"]
    )
    stage_gate = bool(early_pass and late_pass and specificity_pass)

    scores = pd.DataFrame(score_rows).drop_duplicates(["cohort", "contrast_id", "sample"])
    sample_metadata = pd.concat(
        [metadata.assign(cohort=cohort) for cohort, (_, metadata, _) in cohorts.items()], ignore_index=True
    )
    tests.to_csv(TABLES / f"human_PBMC_stage_projection_tests_{DATE}.tsv", sep="\t", index=False)
    meta.to_csv(TABLES / f"human_PBMC_stage_projection_meta_analysis_{DATE}.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / f"human_PBMC_stage_projection_scores_{DATE}.tsv", sep="\t", index=False)
    pd.concat(mapping_rows, ignore_index=True).to_csv(
        TABLES / f"hDRG_to_human_PBMC_signature_mapping_{DATE}.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    sample_metadata.to_csv(METADATA / f"human_PBMC_stage_cohort_metadata_{DATE}.tsv", sep="\t", index=False)
    plot_stage(scores, tests)

    qc = {
        "status": "PASS",
        "frozen_spec": "HUMAN_PBMC_STAGE_VALIDATION_FROZEN_SPEC_2026-08-27.md",
        "seed": SEED,
        "matched_null_iterations": N_NULL,
        "cohort_audits": audits,
        "signature_sha256": sha256(SIGNATURES),
        "human_gene_info_sha256": sha256(info_path),
        "null_summaries": null_summaries,
        "gate_components": {
            "early_transition_pass": early_pass,
            "late_transition_pass": late_pass,
            "GSE185011_DPN_complication_specificity_pass": specificity_pass,
        },
        "human_PBMC_stage_gate_pass": stage_gate,
        "interpretation_boundary": (
            "PBMC validation supports a systemic correlate, not neuronal localization or blood-to-ganglion causality."
        ),
    }
    with (TABLES / f"human_PBMC_stage_projection_qc_{DATE}.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)

    print(tests[["cohort", "transition", "contrast_id", "score_difference", "hedges_g", "exact_one_sided_p", "exact_bh_q", "gene_matched_bh_q", "loo_all_positive", "projection_supportive"]].to_string(index=False))
    print(meta.to_string(index=False))
    print(json.dumps({"human_PBMC_stage_gate_pass": stage_gate, **qc["gate_components"]}, indent=2))


if __name__ == "__main__":
    main()
