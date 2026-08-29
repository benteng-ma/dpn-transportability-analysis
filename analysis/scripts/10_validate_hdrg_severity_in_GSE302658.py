#!/usr/bin/env python3
"""Run the frozen GSE302658 clinical validation of human hDRG signatures."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
import statsmodels.api as sm


PHASE = Path(__file__).resolve().parents[2]
RAW = PHASE / "data" / "raw" / "human_PDN_trial_GSE302658"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
TABLES = PHASE / "results" / "tables"
FIGURES = PHASE / "results" / "figures"
METADATA = PHASE / "metadata"
DATE = "2026-08-27"
SEED = 20260827
N_PERMUTATIONS = 10_000
N_BOOTSTRAPS = 10_000
CPM_THRESHOLD = 1.0
BASELINE_PREVALENCE = 0.20

COUNTS = RAW / "GSE302658_salmon.merged.transcript_counts.tsv.gz"
GENE_INFO = NCBI / "Homo_sapiens.gene_info.gz"
MAPPING = TABLES / f"GSE302658_ensembl_to_ncbi_gene_mapping_{DATE}.tsv.gz"
SAMPLE_METADATA = METADATA / f"GSE302658_clinical_sample_metadata_{DATE}.tsv"
SIGNATURES = TABLES / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"

SEVERITY = "severity_neuron_modhigh_vs_low_nageotte"
LATE_NEURON = "late_neuron_DPN_vs_diabetes"
LATE_ALL = "late_allcell_DPN_vs_diabetes"
EARLY = "early_allcell_diabetes_vs_control"
TESTED_SIGNATURES = [SEVERITY, LATE_NEURON, LATE_ALL, EARLY]

NPSI_SUBSCORES = [
    "burning_superfic_spont_pain_sub_score",
    "pressing_deep_spont_pain_sub_score",
    "paroxysmal_pain_sub_score",
    "evoked_pain_sub_score",
    "paresthesia_dysesthesia_sub_score",
]
SECONDARY_OUTCOMES = [
    "average_pain_nrs_last_12_hours",
    "worst_pain_nrs_last_12_hours",
    *NPSI_SUBSCORES,
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    array = values.to_numpy(dtype=float)
    result = np.full(len(array), np.nan)
    valid = np.isfinite(array)
    if not valid.any():
        return pd.Series(result, index=values.index)
    subset = array[valid]
    order = np.argsort(subset)
    ranked = subset[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    result[np.where(valid)[0]] = restored
    return pd.Series(result, index=values.index)


def load_gene_info(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip", dtype=str, na_filter=False)


def build_symbol_lookup(info: pd.DataFrame) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
    exact: dict[str, str] = {}
    folded_candidates: defaultdict[str, set[str]] = defaultdict(set)
    synonym_candidates: defaultdict[str, set[str]] = defaultdict(set)
    id_to_symbol: dict[str, str] = {}
    for row in info[["GeneID", "Symbol", "Synonyms"]].itertuples(index=False):
        exact[row.Symbol] = row.GeneID
        folded_candidates[row.Symbol.upper()].add(row.GeneID)
        id_to_symbol[row.GeneID] = row.Symbol
        if row.Synonyms and row.Synonyms != "-":
            for synonym in row.Synonyms.split("|"):
                if synonym and synonym != "-":
                    synonym_candidates[synonym.upper()].add(row.GeneID)
    folded = {key: next(iter(ids)) for key, ids in folded_candidates.items() if len(ids) == 1}
    synonyms = {
        key: next(iter(ids))
        for key, ids in synonym_candidates.items()
        if len(ids) == 1 and key not in folded
    }
    return exact, folded, synonyms, id_to_symbol


def resolve_symbol(symbol: object, lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[str | None, str]:
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


def aggregate_transcripts_to_genes(
    counts_path: Path,
    mapping_path: Path,
) -> tuple[pd.DataFrame, pd.Series, dict[str, object]]:
    mapping = pd.read_csv(mapping_path, sep="\t", compression="gzip", dtype=str)
    one_to_one = mapping[mapping["ensembl_mapping_is_one_to_one"].str.lower().eq("true")].copy()
    one_to_one = one_to_one.drop_duplicates("ensembl_gene_id")
    ensembl_to_gene = one_to_one.set_index("ensembl_gene_id")["human_gene_id"].to_dict()
    gene_ids = sorted(set(ensembl_to_gene.values()), key=lambda value: int(value))
    gene_to_row = {gene_id: index for index, gene_id in enumerate(gene_ids)}

    with gzip.open(counts_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    samples = header[2:]
    gene_counts = np.zeros((len(gene_ids), len(samples)), dtype=np.float64)
    library_sizes = np.zeros(len(samples), dtype=np.float64)
    mapped_transcript_rows = 0
    deposited_rows = 0

    for chunk in pd.read_csv(counts_path, sep="\t", compression="gzip", chunksize=5_000):
        deposited_rows += len(chunk)
        values = chunk[samples].to_numpy(dtype=np.float64, copy=False)
        library_sizes += values.sum(axis=0)
        deposited_ensembl = chunk["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        mapped_gene = deposited_ensembl.map(ensembl_to_gene)
        keep = mapped_gene.notna().to_numpy()
        mapped_transcript_rows += int(keep.sum())
        if not keep.any():
            continue
        local = pd.DataFrame(values[keep], columns=samples)
        local.insert(0, "human_gene_id", mapped_gene[keep].to_numpy())
        grouped = local.groupby("human_gene_id", sort=False)[samples].sum()
        indices = np.fromiter((gene_to_row[str(value)] for value in grouped.index), dtype=int)
        gene_counts[indices] += grouped.to_numpy(dtype=np.float64)

    counts = pd.DataFrame(gene_counts, index=pd.Index(gene_ids, name="human_gene_id"), columns=samples)
    nonzero = counts.sum(axis=1) > 0
    counts = counts.loc[nonzero]
    audit = {
        "deposited_transcript_rows": deposited_rows,
        "mapped_transcript_rows": mapped_transcript_rows,
        "mapped_positive_gene_count": int(nonzero.sum()),
        "sample_count": len(samples),
        "library_size_min": float(library_sizes.min()),
        "library_size_median": float(np.median(library_sizes)),
        "library_size_max": float(library_sizes.max()),
    }
    return counts, pd.Series(library_sizes, index=samples, name="salmon_library_size"), audit


def load_signatures(
    lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> pd.DataFrame:
    signatures = pd.read_csv(SIGNATURES, sep="\t")
    signatures = signatures[signatures["contrast_id"].isin(TESTED_SIGNATURES)].copy()
    resolved = signatures["gene"].map(lambda value: resolve_symbol(value, lookup))
    signatures["human_gene_id"] = [item[0] for item in resolved]
    signatures["symbol_mapping_method"] = [item[1] for item in resolved]
    signatures = signatures.sort_values(["contrast_id", "direction", "p_val_adj", "gene"])
    signatures["duplicate_resolved_gene_direction"] = signatures.duplicated(
        ["contrast_id", "direction", "human_gene_id"], keep="first"
    ) & signatures["human_gene_id"].notna()
    return signatures


def standardized_rank(values: np.ndarray) -> np.ndarray:
    ranks = stats.rankdata(values, method="average")
    centered = ranks - ranks.mean()
    scale = centered.std(ddof=1)
    if scale <= 0:
        return np.zeros_like(centered, dtype=float)
    return centered / scale


def correlation_test(
    x: pd.Series,
    y: pd.Series,
    rng: np.random.Generator,
    n_permutations: int = N_PERMUTATIONS,
    n_bootstraps: int = N_BOOTSTRAPS,
) -> dict[str, object]:
    joined = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    xv = joined["x"].to_numpy(dtype=float)
    yv = joined["y"].to_numpy(dtype=float)
    if len(xv) < 10 or np.unique(xv).size < 2 or np.unique(yv).size < 2:
        return {
            "n": len(xv),
            "rho": np.nan,
            "asymptotic_two_sided_p": np.nan,
            "permutation_positive_p": np.nan,
            "permutation_two_sided_p": np.nan,
            "bootstrap_ci_low": np.nan,
            "bootstrap_ci_high": np.nan,
            "loo_all_positive": False,
            "loo_min_rho": np.nan,
            "loo_max_rho": np.nan,
        }
    rho, asymptotic_p = stats.spearmanr(xv, yv)
    xrank = standardized_rank(xv)
    yrank = standardized_rank(yv)
    denominator = len(xv) - 1
    null = np.empty(n_permutations, dtype=float)
    for iteration in range(n_permutations):
        null[iteration] = float(np.dot(xrank, rng.permutation(yrank)) / denominator)
    positive_p = (1 + int(np.sum(null >= rho - 1e-12))) / (n_permutations + 1)
    two_sided_p = (1 + int(np.sum(np.abs(null) >= abs(rho) - 1e-12))) / (n_permutations + 1)

    bootstrap = np.empty(n_bootstraps, dtype=float)
    for iteration in range(n_bootstraps):
        indices = rng.integers(0, len(xv), size=len(xv))
        estimate = stats.spearmanr(xv[indices], yv[indices]).statistic
        bootstrap[iteration] = estimate if np.isfinite(estimate) else 0.0
    loo = np.asarray(
        [stats.spearmanr(np.delete(xv, index), np.delete(yv, index)).statistic for index in range(len(xv))],
        dtype=float,
    )
    return {
        "n": len(xv),
        "rho": float(rho),
        "asymptotic_two_sided_p": float(asymptotic_p),
        "permutation_positive_p": float(positive_p),
        "permutation_two_sided_p": float(two_sided_p),
        "bootstrap_ci_low": float(np.quantile(bootstrap, 0.025)),
        "bootstrap_ci_high": float(np.quantile(bootstrap, 0.975)),
        "loo_all_positive": bool(np.all(loo > 0)),
        "loo_min_rho": float(np.nanmin(loo)),
        "loo_max_rho": float(np.nanmax(loo)),
    }


def fit_hc3(
    data: pd.DataFrame,
    outcome: str,
    exposure: str,
    covariates: list[str],
    rank_outcome: bool,
) -> tuple[dict[str, object], pd.DataFrame]:
    columns = [outcome, exposure, *covariates]
    frame = data[columns].dropna().copy()
    if rank_outcome:
        frame[outcome] = stats.rankdata(frame[outcome].to_numpy(dtype=float), method="average")
    continuous = [column for column in [exposure, "age_years", "bmi", "baseline_npsi", "baseline_score", "actual_visit8_day"] if column in frame]
    for column in continuous:
        std = frame[column].std(ddof=1)
        if std > 0:
            frame[column] = (frame[column] - frame[column].mean()) / std
    categorical = [column for column in covariates if column in {"sex", "randomized_treatment"}]
    numeric_covariates = [column for column in covariates if column not in categorical]
    design = pd.concat(
        [
            frame[[exposure, *numeric_covariates]].astype(float),
            pd.get_dummies(frame[categorical], drop_first=True, dtype=float) if categorical else pd.DataFrame(index=frame.index),
        ],
        axis=1,
    )
    design = sm.add_constant(design, has_constant="add")
    model = sm.OLS(frame[outcome].astype(float), design).fit(cov_type="HC3")
    coefficients = pd.DataFrame(
        {
            "term": model.params.index,
            "coefficient": model.params.to_numpy(dtype=float),
            "standard_error_hc3": model.bse.to_numpy(dtype=float),
            "two_sided_p_hc3": model.pvalues.to_numpy(dtype=float),
            "ci_low_hc3": model.conf_int()[0].to_numpy(dtype=float),
            "ci_high_hc3": model.conf_int()[1].to_numpy(dtype=float),
        }
    )
    row = coefficients.set_index("term").loc[exposure]
    summary = {
        "n": int(model.nobs),
        "coefficient": float(row["coefficient"]),
        "standard_error_hc3": float(row["standard_error_hc3"]),
        "two_sided_p_hc3": float(row["two_sided_p_hc3"]),
        "ci_low_hc3": float(row["ci_low_hc3"]),
        "ci_high_hc3": float(row["ci_high_hc3"]),
        "r_squared": float(model.rsquared),
    }
    return summary, coefficients


def matched_signature_null(
    rank_expression: pd.DataFrame,
    mean_baseline_cpm: pd.Series,
    up_ids: list[str],
    down_ids: list[str],
    outcome: pd.Series,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float, float, dict[str, object]]:
    baseline_samples = outcome.dropna().index
    expression = rank_expression.loc[:, baseline_samples]
    y = outcome.loc[baseline_samples].to_numpy(dtype=float)
    yrank = standardized_rank(y)
    denominator = len(y) - 1
    universe = pd.DataFrame({"mean_cpm": mean_baseline_cpm.loc[rank_expression.index]})
    universe["decile"] = pd.qcut(universe["mean_cpm"].rank(method="first"), 10, labels=False).astype(int)
    signature_ids = set(up_ids) | set(down_ids)
    eligible = universe.loc[~universe.index.isin(signature_ids)]
    row_lookup = {gene_id: index for index, gene_id in enumerate(rank_expression.index)}
    pools = {
        decile: np.asarray([row_lookup[gene_id] for gene_id in frame.index], dtype=int)
        for decile, frame in eligible.groupby("decile")
    }
    up_counts = universe.loc[up_ids, "decile"].value_counts().sort_index().to_dict()
    down_counts = universe.loc[down_ids, "decile"].value_counts().sort_index().to_dict()
    matrix = rank_expression.to_numpy(dtype=float)[:, rank_expression.columns.get_indexer(baseline_samples)]
    null = np.empty(N_PERMUTATIONS, dtype=float)
    for iteration in range(N_PERMUTATIONS):
        sampled_up: list[int] = []
        sampled_down: list[int] = []
        for decile in range(10):
            n_up = int(up_counts.get(decile, 0))
            n_down = int(down_counts.get(decile, 0))
            needed = n_up + n_down
            if needed == 0:
                continue
            pool = pools[decile]
            if needed > len(pool):
                raise RuntimeError(f"Insufficient matched-null pool in decile {decile}: need {needed}, have {len(pool)}")
            selected = rng.choice(pool, size=needed, replace=False)
            sampled_up.extend(selected[:n_up])
            sampled_down.extend(selected[n_up:])
        random_score = matrix[np.asarray(sampled_up)].mean(axis=0) - matrix[np.asarray(sampled_down)].mean(axis=0)
        null[iteration] = float(np.dot(standardized_rank(random_score), yrank) / denominator)
    observed_score = expression.loc[up_ids].mean(axis=0) - expression.loc[down_ids].mean(axis=0)
    observed_rho = float(stats.spearmanr(observed_score, outcome.loc[baseline_samples]).statistic)
    empirical_p = (1 + int(np.sum(null >= observed_rho - 1e-12))) / (N_PERMUTATIONS + 1)
    audit = {
        "eligible_null_gene_count": int(len(eligible)),
        "up_counts_by_expression_decile": {str(key): int(value) for key, value in up_counts.items()},
        "down_counts_by_expression_decile": {str(key): int(value) for key, value in down_counts.items()},
        "null_quantiles": {str(q): float(np.quantile(null, q)) for q in [0, 0.025, 0.5, 0.975, 1]},
    }
    return null, observed_rho, float(empirical_p), audit


def make_paired_frame(scores: pd.DataFrame, metadata: pd.DataFrame, signature: str, outcome: str) -> pd.DataFrame:
    source = metadata.merge(
        scores.loc[scores["contrast_id"] == signature, ["matrix_sample", "score"]],
        on="matrix_sample",
        how="inner",
    )
    score_wide = source.pivot(index="subject_id", columns="visit", values="score")
    outcome_wide = source.pivot(index="subject_id", columns="visit", values=outcome)
    base = source[source["visit"] == "Visit 3"].drop_duplicates("subject_id").set_index("subject_id")
    follow = source[source["visit"] == "Visit 8"].drop_duplicates("subject_id").set_index("subject_id")
    frame = pd.DataFrame(index=score_wide.index.intersection(outcome_wide.index))
    frame["baseline_score"] = score_wide.get("Visit 3")
    frame["visit8_score"] = score_wide.get("Visit 8")
    frame["delta_score"] = frame["visit8_score"] - frame["baseline_score"]
    frame["baseline_outcome"] = outcome_wide.get("Visit 3")
    frame["visit8_outcome"] = outcome_wide.get("Visit 8")
    frame["delta_outcome"] = frame["visit8_outcome"] - frame["baseline_outcome"]
    for column in ["randomized_treatment", "age_years", "sex", "bmi"]:
        frame[column] = base[column]
    frame["actual_visit8_day"] = follow["study_day"]
    return frame.reset_index()


def plot_results(
    metadata_scores: pd.DataFrame,
    paired_primary: pd.DataFrame,
    tests: pd.DataFrame,
    coverage: pd.DataFrame,
    primary_summary: dict[str, object],
    longitudinal_summary: dict[str, object],
) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    colors = {"Placebo": "#777777", "AZD2423 20 mg": "#3B82B4", "AZD2423 150 mg": "#D95F59"}
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 9.2))

    baseline = metadata_scores[
        (metadata_scores["contrast_id"] == SEVERITY) & (metadata_scores["visit"] == "Visit 3")
    ].dropna(subset=["score", "npsi_total_score"])
    ax = axes[0, 0]
    for treatment, frame in baseline.groupby("randomized_treatment"):
        ax.scatter(frame["score"], frame["npsi_total_score"], s=34, alpha=0.78, color=colors[treatment], label=treatment)
    if len(baseline) > 1:
        xline = np.linspace(baseline["score"].min(), baseline["score"].max(), 100)
        slope, intercept = np.polyfit(baseline["score"], baseline["npsi_total_score"], 1)
        ax.plot(xline, intercept + slope * xline, color="black", linewidth=1.4)
    ax.set_xlabel("Frozen hDRG Nageotte severity score in blood")
    ax.set_ylabel("Baseline NPSI total")
    ax.set_title("A  Baseline clinical transfer", loc="left", fontweight="bold")
    ax.text(
        0.03,
        0.97,
        f"n={primary_summary['n']}; rho={primary_summary['rho']:.3f}\n"
        f"permutation P={primary_summary['permutation_positive_p']:.3g}\n"
        f"matched-set P={primary_summary['matched_signature_p']:.3g}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )

    ax = axes[0, 1]
    complete = paired_primary.dropna(subset=["delta_score", "delta_outcome"])
    for treatment, frame in complete.groupby("randomized_treatment"):
        ax.scatter(frame["delta_score"], frame["delta_outcome"], s=38, alpha=0.8, color=colors[treatment], label=treatment)
    if len(complete) > 1:
        xline = np.linspace(complete["delta_score"].min(), complete["delta_score"].max(), 100)
        slope, intercept = np.polyfit(complete["delta_score"], complete["delta_outcome"], 1)
        ax.plot(xline, intercept + slope * xline, color="black", linewidth=1.4)
    ax.axhline(0, color="#BBBBBB", linewidth=0.8)
    ax.axvline(0, color="#BBBBBB", linewidth=0.8)
    ax.set_xlabel("Change in hDRG severity score")
    ax.set_ylabel("Change in NPSI total")
    ax.set_title("B  Within-person clinical covariation", loc="left", fontweight="bold")
    ax.text(
        0.03,
        0.97,
        f"n={longitudinal_summary['n']}; rho={longitudinal_summary['rho']:.3f}\n"
        f"permutation P={longitudinal_summary['permutation_positive_p']:.3g}",
        transform=ax.transAxes,
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    ax.legend(frameon=False, fontsize=8, loc="best")

    ax = axes[1, 0]
    display = tests[
        (tests["contrast_id"] == SEVERITY)
        & tests["endpoint"].isin(["npsi_total_score", "average_pain_nrs_last_12_hours", "worst_pain_nrs_last_12_hours", *NPSI_SUBSCORES])
    ].copy()
    display["label"] = display["window"].map({"baseline": "B", "change": "Delta"}) + ": " + display["endpoint"].replace(
        {
            "npsi_total_score": "NPSI total",
            "average_pain_nrs_last_12_hours": "Average pain",
            "worst_pain_nrs_last_12_hours": "Worst pain",
            "burning_superfic_spont_pain_sub_score": "Burning",
            "pressing_deep_spont_pain_sub_score": "Pressing/deep",
            "paroxysmal_pain_sub_score": "Paroxysmal",
            "evoked_pain_sub_score": "Evoked",
            "paresthesia_dysesthesia_sub_score": "Paresthesia",
        }
    )
    display = display.sort_values(["window", "endpoint"], ascending=[False, True]).reset_index(drop=True)
    y_positions = np.arange(len(display))
    point_colors = display["window"].map({"baseline": "#3B82B4", "change": "#D95F59"})
    ax.hlines(y_positions, display["bootstrap_ci_low"], display["bootstrap_ci_high"], color=point_colors, linewidth=1.5)
    ax.scatter(display["rho"], y_positions, color=point_colors, s=32, zorder=3)
    ax.axvline(0, color="#777777", linewidth=0.8)
    ax.set_yticks(y_positions, display["label"], fontsize=8)
    ax.set_xlabel("Spearman rho (bootstrap 95% CI)")
    ax.set_title("C  Clinical outcome profile", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 1]
    cov = coverage[coverage["contrast_id"].isin(TESTED_SIGNATURES)].copy()
    cov["short"] = cov["contrast_id"].map(
        {SEVERITY: "Nageotte severity", LATE_NEURON: "Late neuronal", LATE_ALL: "Late all-cell", EARLY: "Early diabetes"}
    )
    pivot = cov.pivot(index="short", columns="direction", values="expressed_fraction_of_original").reindex(
        ["Nageotte severity", "Late neuronal", "Late all-cell", "Early diabetes"]
    )
    x = np.arange(len(pivot))
    ax.bar(x - 0.18, pivot["up"] * 100, width=0.36, color="#C44E52", label="Up")
    ax.bar(x + 0.18, pivot["down"] * 100, width=0.36, color="#4C72B0", label="Down")
    ax.axhline(20, color="#444444", linestyle="--", linewidth=0.9, label="Frozen minimum")
    ax.set_xticks(x, pivot.index, rotation=20, ha="right")
    ax.set_ylabel("Original signature retained in blood (%)")
    ax.set_title("D  Prespecified blood-expression coverage", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    for axis in axes.flat:
        axis.grid(color="#EEEEEE", linewidth=0.6, zorder=0)
    fig.suptitle("Clinical test of a human sensory-ganglion degeneration program in GSE302658", fontweight="bold", y=0.995)
    fig.tight_layout()
    for extension in ["png", "pdf"]:
        fig.savefig(FIGURES / f"GSE302658_hDRG_severity_clinical_validation_{DATE}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(SAMPLE_METADATA, sep="\t", dtype={"subject_id": str, "matrix_sample": str})
    info = load_gene_info(GENE_INFO)
    lookup = build_symbol_lookup(info)
    counts, library_sizes, aggregation_audit = aggregate_transcripts_to_genes(COUNTS, MAPPING)
    if list(counts.columns) != metadata["matrix_sample"].tolist():
        raise RuntimeError("Aggregated count columns no longer match audited metadata order")

    baseline_samples = metadata.loc[metadata["visit"] == "Visit 3", "matrix_sample"].tolist()
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    minimum_baseline_samples = math.ceil(BASELINE_PREVALENCE * len(baseline_samples))
    expressed = (cpm[baseline_samples] >= CPM_THRESHOLD).sum(axis=1) >= minimum_baseline_samples
    filtered_counts = counts.loc[expressed].copy()
    filtered_cpm = cpm.loc[expressed].copy()
    rank_expression = filtered_counts.rank(axis=0, method="average", pct=True) - 0.5
    mean_baseline_cpm = filtered_cpm[baseline_samples].mean(axis=1)

    signatures = load_signatures(lookup)
    signatures["resolved_current_symbol"] = signatures["human_gene_id"].map(lookup[3])
    signatures["mapped_to_deposited_gene"] = signatures["human_gene_id"].isin(counts.index)
    signatures["passes_frozen_blood_expression_filter"] = signatures["human_gene_id"].isin(rank_expression.index)
    unique_signatures = signatures[
        signatures["human_gene_id"].notna() & ~signatures["duplicate_resolved_gene_direction"]
    ].copy()

    coverage_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    score_vectors: dict[str, pd.Series] = {}
    for contrast in TESTED_SIGNATURES:
        source_all = signatures[signatures["contrast_id"] == contrast]
        source = unique_signatures[unique_signatures["contrast_id"] == contrast]
        for direction in ["up", "down"]:
            original_n = int((source_all["direction"] == direction).sum())
            resolved = source[source["direction"] == direction]
            expressed_ids = resolved.loc[
                resolved["human_gene_id"].isin(rank_expression.index), "human_gene_id"
            ].drop_duplicates().tolist()
            coverage_rows.append(
                {
                    "contrast_id": contrast,
                    "direction": direction,
                    "original_gene_n": original_n,
                    "resolved_unique_gene_n": int(resolved["human_gene_id"].nunique()),
                    "deposited_gene_n": int(resolved["human_gene_id"].isin(counts.index).sum()),
                    "expressed_gene_n": len(expressed_ids),
                    "expressed_fraction_of_original": len(expressed_ids) / original_n,
                }
            )
        up_ids = source.loc[
            (source["direction"] == "up") & source["human_gene_id"].isin(rank_expression.index), "human_gene_id"
        ].drop_duplicates().tolist()
        down_ids = source.loc[
            (source["direction"] == "down") & source["human_gene_id"].isin(rank_expression.index), "human_gene_id"
        ].drop_duplicates().tolist()
        score = rank_expression.loc[up_ids].mean(axis=0) - rank_expression.loc[down_ids].mean(axis=0)
        score_vectors[contrast] = score
        for sample, value in score.items():
            score_rows.append({"contrast_id": contrast, "matrix_sample": sample, "score": float(value)})

    coverage = pd.DataFrame(coverage_rows)
    coverage["passes_count_20"] = coverage["expressed_gene_n"] >= 20
    coverage["passes_fraction_0_20"] = coverage["expressed_fraction_of_original"] >= 0.20
    coverage["direction_coverage_pass"] = coverage["passes_count_20"] & coverage["passes_fraction_0_20"]
    scores = pd.DataFrame(score_rows)
    metadata_scores = metadata.merge(scores, on="matrix_sample", how="left")

    primary_coverage_pass = bool(
        coverage.loc[coverage["contrast_id"] == SEVERITY, "direction_coverage_pass"].all()
    )
    baseline_primary = metadata_scores[
        (metadata_scores["contrast_id"] == SEVERITY) & (metadata_scores["visit"] == "Visit 3")
    ].set_index("matrix_sample")
    rng = np.random.default_rng(SEED)
    primary_correlation = correlation_test(
        baseline_primary["score"], baseline_primary["npsi_total_score"], rng
    )

    primary_source = unique_signatures[unique_signatures["contrast_id"] == SEVERITY]
    primary_up = primary_source.loc[
        (primary_source["direction"] == "up") & primary_source["human_gene_id"].isin(rank_expression.index), "human_gene_id"
    ].drop_duplicates().tolist()
    primary_down = primary_source.loc[
        (primary_source["direction"] == "down") & primary_source["human_gene_id"].isin(rank_expression.index), "human_gene_id"
    ].drop_duplicates().tolist()
    baseline_outcome = baseline_primary["npsi_total_score"]
    null, matched_observed_rho, matched_p, matched_audit = matched_signature_null(
        rank_expression[baseline_samples], mean_baseline_cpm, primary_up, primary_down, baseline_outcome, rng
    )
    if not np.isclose(primary_correlation["rho"], matched_observed_rho, atol=1e-12):
        raise RuntimeError("Observed primary rho differs between correlation and matched-null implementations")

    baseline_model_data = baseline_primary.reset_index()
    baseline_model, baseline_coefficients = fit_hc3(
        baseline_model_data,
        outcome="npsi_total_score",
        exposure="score",
        covariates=["age_years", "sex", "bmi"],
        rank_outcome=True,
    )
    baseline_gate = bool(
        primary_coverage_pass
        and primary_correlation["rho"] > 0
        and primary_correlation["permutation_positive_p"] < 0.05
        and matched_p < 0.05
        and baseline_model["coefficient"] > 0
        and baseline_model["two_sided_p_hc3"] < 0.05
    )

    paired_primary = make_paired_frame(scores, metadata, SEVERITY, "npsi_total_score")
    paired_primary["baseline_npsi"] = paired_primary["baseline_outcome"]
    longitudinal_correlation = correlation_test(
        paired_primary.set_index("subject_id")["delta_score"],
        paired_primary.set_index("subject_id")["delta_outcome"],
        rng,
    )
    longitudinal_model, longitudinal_coefficients = fit_hc3(
        paired_primary,
        outcome="visit8_outcome",
        exposure="visit8_score",
        covariates=[
            "baseline_npsi",
            "baseline_score",
            "randomized_treatment",
            "age_years",
            "sex",
            "bmi",
            "actual_visit8_day",
        ],
        rank_outcome=False,
    )
    longitudinal_confirmation = bool(
        longitudinal_correlation["rho"] > 0
        and longitudinal_correlation["permutation_positive_p"] < 0.05
        and longitudinal_model["coefficient"] > 0
        and longitudinal_model["two_sided_p_hc3"] < 0.05
    )

    test_rows: list[dict[str, object]] = []
    primary_row = {
        "test_family": "primary",
        "contrast_id": SEVERITY,
        "window": "baseline",
        "endpoint": "npsi_total_score",
        **primary_correlation,
        "matched_signature_p": matched_p,
        "hc3_exposure_coefficient": baseline_model["coefficient"],
        "hc3_exposure_p": baseline_model["two_sided_p_hc3"],
        "frozen_gate_pass": baseline_gate,
    }
    test_rows.append(primary_row)
    longitudinal_row = {
        "test_family": "longitudinal_primary",
        "contrast_id": SEVERITY,
        "window": "change",
        "endpoint": "npsi_total_score",
        **longitudinal_correlation,
        "matched_signature_p": np.nan,
        "hc3_exposure_coefficient": longitudinal_model["coefficient"],
        "hc3_exposure_p": longitudinal_model["two_sided_p_hc3"],
        "frozen_gate_pass": longitudinal_confirmation,
    }
    test_rows.append(longitudinal_row)

    for outcome in SECONDARY_OUTCOMES:
        baseline_result = correlation_test(
            baseline_primary["score"], baseline_primary[outcome], rng, n_bootstraps=2_000
        )
        test_rows.append(
            {
                "test_family": "baseline_secondary",
                "contrast_id": SEVERITY,
                "window": "baseline",
                "endpoint": outcome,
                **baseline_result,
                "matched_signature_p": np.nan,
                "hc3_exposure_coefficient": np.nan,
                "hc3_exposure_p": np.nan,
                "frozen_gate_pass": False,
            }
        )
        paired = make_paired_frame(scores, metadata, SEVERITY, outcome).set_index("subject_id")
        change_result = correlation_test(
            paired["delta_score"], paired["delta_outcome"], rng, n_bootstraps=2_000
        )
        test_rows.append(
            {
                "test_family": "longitudinal_secondary",
                "contrast_id": SEVERITY,
                "window": "change",
                "endpoint": outcome,
                **change_result,
                "matched_signature_p": np.nan,
                "hc3_exposure_coefficient": np.nan,
                "hc3_exposure_p": np.nan,
                "frozen_gate_pass": False,
            }
        )

    for contrast in [LATE_NEURON, LATE_ALL]:
        baseline = metadata_scores[
            (metadata_scores["contrast_id"] == contrast) & (metadata_scores["visit"] == "Visit 3")
        ].set_index("matrix_sample")
        result = correlation_test(baseline["score"], baseline["npsi_total_score"], rng, n_bootstraps=2_000)
        test_rows.append(
            {
                "test_family": "baseline_secondary",
                "contrast_id": contrast,
                "window": "baseline",
                "endpoint": "npsi_total_score",
                **result,
                "matched_signature_p": np.nan,
                "hc3_exposure_coefficient": np.nan,
                "hc3_exposure_p": np.nan,
                "frozen_gate_pass": False,
            }
        )
        paired = make_paired_frame(scores, metadata, contrast, "npsi_total_score").set_index("subject_id")
        result = correlation_test(paired["delta_score"], paired["delta_outcome"], rng, n_bootstraps=2_000)
        test_rows.append(
            {
                "test_family": "longitudinal_secondary",
                "contrast_id": contrast,
                "window": "change",
                "endpoint": "npsi_total_score",
                **result,
                "matched_signature_p": np.nan,
                "hc3_exposure_coefficient": np.nan,
                "hc3_exposure_p": np.nan,
                "frozen_gate_pass": False,
            }
        )

    negative_baseline = metadata_scores[
        (metadata_scores["contrast_id"] == EARLY) & (metadata_scores["visit"] == "Visit 3")
    ].set_index("matrix_sample")
    result = correlation_test(negative_baseline["score"], negative_baseline["npsi_total_score"], rng, n_bootstraps=2_000)
    test_rows.append(
        {
            "test_family": "negative_control",
            "contrast_id": EARLY,
            "window": "baseline",
            "endpoint": "npsi_total_score",
            **result,
            "matched_signature_p": np.nan,
            "hc3_exposure_coefficient": np.nan,
            "hc3_exposure_p": np.nan,
            "frozen_gate_pass": False,
        }
    )
    negative_paired = make_paired_frame(scores, metadata, EARLY, "npsi_total_score").set_index("subject_id")
    result = correlation_test(negative_paired["delta_score"], negative_paired["delta_outcome"], rng, n_bootstraps=2_000)
    test_rows.append(
        {
            "test_family": "negative_control",
            "contrast_id": EARLY,
            "window": "change",
            "endpoint": "npsi_total_score",
            **result,
            "matched_signature_p": np.nan,
            "hc3_exposure_coefficient": np.nan,
            "hc3_exposure_p": np.nan,
            "frozen_gate_pass": False,
        }
    )

    tests = pd.DataFrame(test_rows)
    tests["secondary_bh_q"] = np.nan
    for family in ["baseline_secondary", "longitudinal_secondary"]:
        mask = tests["test_family"] == family
        tests.loc[mask, "secondary_bh_q"] = bh_adjust(tests.loc[mask, "permutation_two_sided_p"])

    baseline_rank_matrix = rank_expression[baseline_samples].T
    pca = PCA(n_components=3, svd_solver="randomized", random_state=SEED)
    pc_values = pca.fit_transform(baseline_rank_matrix)
    pc_frame = pd.DataFrame(pc_values, index=baseline_samples, columns=["PC1", "PC2", "PC3"])
    pc_frame["library_size"] = library_sizes.loc[baseline_samples]
    confound = baseline_primary[["score", "age_years", "bmi", "sex"]].join(pc_frame)
    confound_rows = []
    for variable in ["age_years", "bmi", "library_size", "PC1", "PC2", "PC3"]:
        rho, p = stats.spearmanr(confound["score"], confound[variable], nan_policy="omit")
        confound_rows.append({"variable": variable, "association": "spearman", "estimate": rho, "two_sided_p": p})
    female = confound["sex"].eq("F").astype(int)
    estimate, p = stats.pointbiserialr(female, confound["score"])
    confound_rows.append({"variable": "sex_F", "association": "point_biserial", "estimate": estimate, "two_sided_p": p})
    confound_tests = pd.DataFrame(confound_rows)

    primary_summary = {**primary_correlation, "matched_signature_p": matched_p}
    longitudinal_summary = dict(longitudinal_correlation)
    if baseline_gate:
        interpretation_tier = "validated_clinical_transfer"
    elif longitudinal_confirmation:
        interpretation_tier = "supportive_longitudinal_signal"
    else:
        interpretation_tier = "compartment_boundary"

    counts_output = filtered_counts.copy()
    counts_output.insert(0, "current_symbol", counts_output.index.map(lookup[3]))
    counts_output.to_csv(
        TABLES / f"GSE302658_gene_level_estimated_counts_filtered_{DATE}.tsv.gz",
        sep="\t",
        compression="gzip",
        float_format="%.8g",
    )
    signatures.to_csv(TABLES / f"GSE302658_hDRG_signature_mapping_{DATE}.tsv", sep="\t", index=False)
    coverage.to_csv(TABLES / f"GSE302658_hDRG_signature_coverage_{DATE}.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / f"GSE302658_hDRG_signature_scores_{DATE}.tsv", sep="\t", index=False)
    tests.to_csv(TABLES / f"GSE302658_clinical_signature_tests_{DATE}.tsv", sep="\t", index=False)
    confound_tests.to_csv(TABLES / f"GSE302658_primary_score_confound_checks_{DATE}.tsv", sep="\t", index=False)
    baseline_coefficients.assign(model="baseline_primary").to_csv(
        TABLES / f"GSE302658_baseline_primary_HC3_model_{DATE}.tsv", sep="\t", index=False
    )
    longitudinal_coefficients.assign(model="longitudinal_primary").to_csv(
        TABLES / f"GSE302658_longitudinal_primary_HC3_model_{DATE}.tsv", sep="\t", index=False
    )
    pd.DataFrame({"matched_random_signature_rho": null}).to_csv(
        TABLES / f"GSE302658_primary_matched_signature_null_{DATE}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    qc = {
        "dataset": "GSE302658",
        "analysis_date": DATE,
        "random_seed": SEED,
        "frozen_spec": str(PHASE / f"CLINICAL_PDN_GSE302658_FROZEN_SPEC_{DATE}.md"),
        "input_sha256": {
            "counts": sha256(COUNTS),
            "sample_metadata": sha256(SAMPLE_METADATA),
            "signatures": sha256(SIGNATURES),
            "human_gene_info": sha256(GENE_INFO),
        },
        "aggregation": aggregation_audit,
        "filter": {
            "cpm_threshold": CPM_THRESHOLD,
            "baseline_prevalence": BASELINE_PREVALENCE,
            "minimum_baseline_sample_count": minimum_baseline_samples,
            "retained_gene_count": int(len(rank_expression)),
        },
        "pca_baseline_variance_explained": {
            f"PC{index + 1}": float(value) for index, value in enumerate(pca.explained_variance_ratio_)
        },
        "matched_signature_null": matched_audit,
        "primary": {
            "coverage_pass": primary_coverage_pass,
            "correlation": primary_correlation,
            "matched_signature_p": matched_p,
            "hc3_model": baseline_model,
            "frozen_baseline_gate_pass": baseline_gate,
        },
        "longitudinal": {
            "correlation": longitudinal_correlation,
            "hc3_model": longitudinal_model,
            "confirmation_pass": longitudinal_confirmation,
        },
        "interpretation_tier": interpretation_tier,
    }
    with (TABLES / f"GSE302658_clinical_validation_qc_{DATE}.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, ensure_ascii=False)

    plot_results(metadata_scores, paired_primary, tests, coverage, primary_summary, longitudinal_summary)
    print(tests[["test_family", "contrast_id", "window", "endpoint", "n", "rho", "permutation_positive_p", "permutation_two_sided_p", "secondary_bh_q", "frozen_gate_pass"]].to_string(index=False))
    print(json.dumps({"baseline_gate": baseline_gate, "longitudinal_confirmation": longitudinal_confirmation, "interpretation_tier": interpretation_tier}, indent=2))


if __name__ == "__main__":
    main()
