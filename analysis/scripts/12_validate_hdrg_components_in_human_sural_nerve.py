#!/usr/bin/env python3
"""Test frozen hDRG components in independent human sural-nerve DPN cohorts."""

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
RAW = PHASE / "data" / "raw" / "human_sural_nerve_JCI184075" / "supplementary"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
SIGNATURES = TABLES / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"
DATE = "2026-08-27"
RANDOM_SEED = 20260827
N_NULL = 10_000

FILES = {
    "sample_info": RAW / "jci-135-184075-s171.xlsx",
    "quantile_tpm": RAW / "jci-135-184075-s175.xlsx",
    "dpn_control": RAW / "jci-135-184075-s176.xlsx",
    "axon_severity": RAW / "jci-135-184075-s172.xlsx",
}

EXPECTED_COMPONENT_COUNTS = {
    "late_shared_concordant_neuronal_core": {"all": 196, "up": 77, "down": 119},
    "late_allcell_residual": {"all": 2951, "up": 1339, "down": 1612},
    "late_neuron_residual": {"all": 522, "up": 215, "down": 307},
    "late_directionally_opposed_overlap": {"all": 26, "up": 0, "down": 0},
    "severity_neuron_shared_concordant_core": {"all": 137, "up": 104, "down": 33},
    "severity_neuron_residual": {"all": 348, "up": 230, "down": 118},
    "severity_directionally_opposed_overlap": {"all": 24, "up": 0, "down": 0},
}

DISPLAY = {
    "late_shared_concordant_neuronal_core": "Late shared neuronal core",
    "late_allcell_residual": "Late all-cell residual",
    "late_neuron_residual": "Late neuron residual",
    "severity_neuron_shared_concordant_core": "Severity shared neuronal core",
    "severity_neuron_residual": "Severity residual",
    "original_late_allcell": "Original late all-cell",
    "original_late_neuron": "Original late neuron",
    "original_severity": "Original neuron severity",
    "original_early_allcell": "Original early all-cell",
    "original_xenium": "Original Xenium",
}

DPN_MODULES = [
    "late_shared_concordant_neuronal_core",
    "late_neuron_residual",
    "late_allcell_residual",
    "severity_neuron_shared_concordant_core",
    "severity_neuron_residual",
    "original_late_allcell",
    "original_late_neuron",
    "original_severity",
    "original_early_allcell",
    "original_xenium",
]
SEVERITY_MODULES = [
    "severity_neuron_shared_concordant_core",
    "severity_neuron_residual",
    "original_severity",
    "original_late_neuron",
]
PRIMARY_DPN = {
    "late_shared_concordant_neuronal_core",
    "late_neuron_residual",
}
PRIMARY_SEVERITY = {
    "severity_neuron_shared_concordant_core",
    "severity_neuron_residual",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    valid = values.notna()
    result = pd.Series(np.nan, index=values.index, dtype=float)
    if not valid.any():
        return result
    array = values.loc[valid].to_numpy()
    order = np.argsort(array)
    ranked = array[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    result.loc[valid] = restored
    return result


def load_gene_info(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, na_filter=False)
    return frame[["GeneID", "Symbol", "Synonyms"]].copy()


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
        if row.Synonyms and row.Synonyms != "-":
            for synonym in str(row.Synonyms).split("|"):
                synonym = synonym.strip()
                if synonym and synonym != "-":
                    synonym_candidates[synonym.upper()].add(gene_id)
    official_folded = {
        symbol: next(iter(ids))
        for symbol, ids in official_folded_candidates.items()
        if len(ids) == 1
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
    frame = frame.copy()
    repaired = frame[gene_column].map(repair_excel_gene)
    frame["gene_input"] = [item[0] for item in repaired]
    frame["gene_recovered_from_excel_date"] = [item[1] for item in repaired]
    resolved = frame["gene_input"].map(
        lambda symbol: resolve_symbol(symbol, official_exact, official_folded, synonym_unique)
    )
    frame["human_gene_id"] = [item[0] for item in resolved]
    frame["mapping_method"] = [item[1] for item in resolved]
    frame["current_human_symbol"] = frame["human_gene_id"].map(id_to_symbol)
    return frame


def make_module(
    frame: pd.DataFrame,
    module_id: str,
    source_contrast: str,
    effect_column: str = "avg_log2FC",
) -> pd.DataFrame:
    result = frame[["gene", "direction", effect_column]].copy()
    result = result.rename(columns={effect_column: "source_log2fc"})
    result["module_id"] = module_id
    result["source_contrast"] = source_contrast
    result["audit_only"] = False
    return result


def build_source_modules() -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict[str, object]]:
    signatures = pd.read_csv(SIGNATURES, sep="\t")
    membership = signatures["primary_signature_member"]
    if membership.dtype != bool:
        membership = membership.astype(str).str.lower().isin(["true", "1", "yes"])
    signatures = signatures[membership].copy()
    signatures["gene"] = signatures["gene"].astype(str).str.strip()
    signatures["direction"] = signatures["direction"].astype(str).str.lower()

    by_contrast = {
        contrast: block.drop_duplicates("gene").set_index("gene", drop=False)
        for contrast, block in signatures.groupby("contrast_id", sort=False)
    }
    late_all = by_contrast["late_allcell_DPN_vs_diabetes"]
    late_neuron = by_contrast["late_neuron_DPN_vs_diabetes"]
    severity = by_contrast["severity_neuron_modhigh_vs_low_nageotte"]
    early = by_contrast["early_allcell_diabetes_vs_control"]
    xenium = by_contrast["xenium_DPN_vs_control"]

    late_overlap = set(late_all.index) & set(late_neuron.index)
    late_same = {
        gene for gene in late_overlap if late_all.at[gene, "direction"] == late_neuron.at[gene, "direction"]
    }
    late_opposed = late_overlap - late_same
    severity_overlap = set(severity.index) & set(late_neuron.index)
    severity_same = {
        gene
        for gene in severity_overlap
        if severity.at[gene, "direction"] == late_neuron.at[gene, "direction"]
    }
    severity_opposed = severity_overlap - severity_same

    modules: dict[str, pd.DataFrame] = {}
    modules["late_shared_concordant_neuronal_core"] = make_module(
        late_neuron.loc[sorted(late_same)],
        "late_shared_concordant_neuronal_core",
        "late_neuron_DPN_vs_diabetes",
    )
    modules["late_allcell_residual"] = make_module(
        late_all.loc[sorted(set(late_all.index) - late_overlap)],
        "late_allcell_residual",
        "late_allcell_DPN_vs_diabetes",
    )
    modules["late_neuron_residual"] = make_module(
        late_neuron.loc[sorted(set(late_neuron.index) - late_overlap)],
        "late_neuron_residual",
        "late_neuron_DPN_vs_diabetes",
    )
    modules["severity_neuron_shared_concordant_core"] = make_module(
        severity.loc[sorted(severity_same)],
        "severity_neuron_shared_concordant_core",
        "severity_neuron_modhigh_vs_low_nageotte",
    )
    modules["severity_neuron_residual"] = make_module(
        severity.loc[sorted(set(severity.index) - severity_overlap)],
        "severity_neuron_residual",
        "severity_neuron_modhigh_vs_low_nageotte",
    )

    opposed_rows: list[pd.DataFrame] = []
    for module_id, genes, left, right in [
        ("late_directionally_opposed_overlap", late_opposed, late_all, late_neuron),
        ("severity_directionally_opposed_overlap", severity_opposed, severity, late_neuron),
    ]:
        block = pd.DataFrame(
            {
                "gene": sorted(genes),
                "direction_source_a": [left.at[g, "direction"] for g in sorted(genes)],
                "direction_source_b": [right.at[g, "direction"] for g in sorted(genes)],
                "source_a_log2fc": [left.at[g, "avg_log2FC"] for g in sorted(genes)],
                "source_b_log2fc": [right.at[g, "avg_log2FC"] for g in sorted(genes)],
                "module_id": module_id,
                "audit_only": True,
            }
        )
        opposed_rows.append(block)

    originals = {
        "original_late_allcell": (late_all, "late_allcell_DPN_vs_diabetes"),
        "original_late_neuron": (late_neuron, "late_neuron_DPN_vs_diabetes"),
        "original_severity": (severity, "severity_neuron_modhigh_vs_low_nageotte"),
        "original_early_allcell": (early, "early_allcell_diabetes_vs_control"),
        "original_xenium": (xenium, "xenium_DPN_vs_control"),
    }
    for module_id, (frame, contrast) in originals.items():
        modules[module_id] = make_module(frame, module_id, contrast)

    counts = {}
    for module_id in EXPECTED_COMPONENT_COUNTS:
        if "opposed" in module_id:
            block = next(frame for frame in opposed_rows if frame["module_id"].iloc[0] == module_id)
            observed = {"all": int(len(block)), "up": 0, "down": 0}
        else:
            block = modules[module_id]
            observed = {
                "all": int(len(block)),
                "up": int((block["direction"] == "up").sum()),
                "down": int((block["direction"] == "down").sum()),
            }
        if observed != EXPECTED_COMPONENT_COUNTS[module_id]:
            raise RuntimeError(
                f"Frozen source-component count mismatch for {module_id}: {observed}"
            )
        counts[module_id] = observed

    component_membership = pd.concat(
        [
            *[frame for key, frame in modules.items() if not key.startswith("original_")],
            *opposed_rows,
        ],
        ignore_index=True,
        sort=False,
    )
    audit = {
        "signature_file_sha256": sha256(SIGNATURES),
        "source_signature_counts": {
            contrast: int(len(frame)) for contrast, frame in by_contrast.items()
        },
        "asserted_component_counts": counts,
    }
    return modules, component_membership, audit


def load_target(
    target_id: str,
    path: Path,
    sheet: str,
    sample_groups: dict[str, str],
    target_orientation_multiplier: int,
    lookups: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    raw = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
    raw.columns = [str(column).strip() for column in raw.columns]
    raw = raw.rename(columns={raw.columns[0]: "source_gene"})
    raw["source_row"] = np.arange(2, len(raw) + 2)
    raw = add_resolution(raw, "source_gene", lookups)
    for column in ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", *sample_groups]:
        raw[column] = pd.to_numeric(raw[column], errors="coerce")
    raw["wald_stat"] = raw["log2FoldChange"] / raw["lfcSE"]
    raw["target_log2fc_oriented"] = raw["log2FoldChange"] * target_orientation_multiplier
    raw["target_wald_oriented"] = raw["wald_stat"] * target_orientation_multiplier
    raw["finite_stat"] = np.isfinite(raw["target_wald_oriented"])
    raw["padj_sort"] = raw["padj"].fillna(np.inf)
    raw["abs_stat"] = raw["target_wald_oriented"].abs().fillna(-np.inf)
    resolved = raw[raw["human_gene_id"].notna()].copy()
    resolved = resolved.sort_values(
        ["human_gene_id", "finite_stat", "padj_sort", "abs_stat", "source_row"],
        ascending=[True, False, True, False, True],
    ).drop_duplicates("human_gene_id", keep="first")
    universe = resolved[
        resolved["finite_stat"]
        & np.isfinite(resolved["baseMean"])
        & (resolved["baseMean"] > 0)
    ].copy()
    expression = universe.set_index("human_gene_id")[list(sample_groups)].copy()
    expression = expression.replace([np.inf, -np.inf], np.nan)
    complete_expression = expression.notna().all(axis=1)
    universe = universe[universe["human_gene_id"].isin(expression.index[complete_expression])].copy()
    expression = expression.loc[universe["human_gene_id"]]
    expression.index.name = "human_gene_id"
    metadata = pd.DataFrame(
        {
            "target_id": target_id,
            "sample_id": list(sample_groups),
            "group": [sample_groups[sample] for sample in sample_groups],
            "primary_expression_source": "DESeq2_normalized_counts",
        }
    )
    audit = {
        "target_id": target_id,
        "source_workbook": path.name,
        "source_workbook_sha256": sha256(path),
        "source_sheet": sheet,
        "deposited_rows": int(len(raw)),
        "resolved_rows": int(raw["human_gene_id"].notna().sum()),
        "unique_resolved_gene_ids": int(len(resolved)),
        "finite_positive_basemean_complete_expression_universe": int(len(universe)),
        "excel_date_cells_recovered": int(raw["gene_recovered_from_excel_date"].sum()),
        "sample_columns": list(sample_groups),
        "group_counts": metadata["group"].value_counts().to_dict(),
        "target_orientation_multiplier": target_orientation_multiplier,
    }
    return universe, expression, metadata, audit


def load_tpm_sensitivity(
    path: Path,
    sample_groups: dict[str, str],
    lookups: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = pd.read_excel(
        path,
        sheet_name="A - Quantile normalized TPMs",
        engine="openpyxl",
    )
    raw.columns = [str(column).strip() for column in raw.columns]
    raw = raw.rename(columns={raw.columns[0]: "source_gene"})
    raw["source_row"] = np.arange(2, len(raw) + 2)
    raw = add_resolution(raw, "source_gene", lookups)
    for sample in sample_groups:
        raw[sample] = pd.to_numeric(raw[sample], errors="coerce")
    raw["complete"] = raw[list(sample_groups)].notna().all(axis=1)
    raw["mean_expression"] = raw[list(sample_groups)].mean(axis=1)
    resolved = raw[raw["human_gene_id"].notna() & raw["complete"]].copy()
    resolved = resolved.sort_values(
        ["human_gene_id", "mean_expression", "source_row"],
        ascending=[True, False, True],
    ).drop_duplicates("human_gene_id", keep="first")
    resolved = resolved[resolved["mean_expression"] > 0]
    matrix = resolved.set_index("human_gene_id")[list(sample_groups)]
    audit = {
        "source_workbook": path.name,
        "source_workbook_sha256": sha256(path),
        "source_sheet": "A - Quantile normalized TPMs",
        "deposited_rows": int(len(raw)),
        "positive_complete_unique_gene_ids": int(len(matrix)),
        "excel_date_cells_recovered": int(raw["gene_recovered_from_excel_date"].sum()),
    }
    return matrix, audit


def map_modules(
    modules: dict[str, pd.DataFrame],
    target_universe: pd.DataFrame,
    target_id: str,
    module_ids: list[str],
    lookups: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    target_columns = [
        "human_gene_id",
        "current_human_symbol",
        "baseMean",
        "target_log2fc_oriented",
        "target_wald_oriented",
        "padj",
    ]
    target = target_universe[target_columns].copy()
    mapped_modules: dict[str, pd.DataFrame] = {}
    mapping_rows = []
    coverage_rows = []
    for module_id in module_ids:
        source = add_resolution(modules[module_id], "gene", lookups)
        source["source_row"] = np.arange(len(source))
        source = source.sort_values(
            ["human_gene_id", "source_row"], ascending=[True, True]
        ).drop_duplicates("human_gene_id", keep="first")
        merged = source.merge(
            target,
            how="left",
            on="human_gene_id",
            suffixes=("_source", "_target"),
        )
        merged["target_id"] = target_id
        merged["mapped_to_target_universe"] = merged["target_wald_oriented"].notna()
        mapping_rows.append(merged)
        mapped = merged[merged["mapped_to_target_universe"]].copy()
        mapped_modules[module_id] = mapped
        counts = {
            direction: int((mapped["direction"] == direction).sum())
            for direction in ["up", "down"]
        }
        coverage_rows.append(
            {
                "target_id": target_id,
                "module_id": module_id,
                "source_n": int(len(source)),
                "source_up_n": int((source["direction"] == "up").sum()),
                "source_down_n": int((source["direction"] == "down").sum()),
                "mapped_n": int(len(mapped)),
                "mapped_up_n": counts["up"],
                "mapped_down_n": counts["down"],
                "coverage_fraction": float(len(mapped) / len(source)) if len(source) else np.nan,
                "minimum_10_per_direction": bool(counts["up"] >= 10 and counts["down"] >= 10),
            }
        )
    return (
        mapped_modules,
        pd.concat(mapping_rows, ignore_index=True),
        pd.DataFrame(coverage_rows),
    )


def matched_null_p(
    mapped: pd.DataFrame,
    target_universe: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    universe = target_universe.copy()
    universe["expression_decile"] = pd.qcut(
        universe["baseMean"].rank(method="first"), 10, labels=False
    ).astype(int)
    decile_by_gene = universe.set_index("human_gene_id")["expression_decile"]
    values_by_decile = {
        decile: block["target_wald_oriented"].to_numpy(dtype=float)
        for decile, block in universe.groupby("expression_decile")
    }
    mapped = mapped.copy()
    mapped["expression_decile"] = mapped["human_gene_id"].map(decile_by_gene)
    up = mapped[mapped["direction"] == "up"]
    down = mapped[mapped["direction"] == "down"]
    observed = float(up["target_wald_oriented"].mean() - down["target_wald_oriented"].mean())
    up_counts = up["expression_decile"].value_counts().to_dict()
    down_counts = down["expression_decile"].value_counts().to_dict()
    null = np.empty(N_NULL, dtype=float)
    for iteration in range(N_NULL):
        up_sum = 0.0
        down_sum = 0.0
        up_n = 0
        down_n = 0
        for decile, count in up_counts.items():
            chosen = rng.choice(values_by_decile[int(decile)], size=int(count), replace=False)
            up_sum += float(chosen.sum())
            up_n += int(count)
        for decile, count in down_counts.items():
            chosen = rng.choice(values_by_decile[int(decile)], size=int(count), replace=False)
            down_sum += float(chosen.sum())
            down_n += int(count)
        null[iteration] = up_sum / up_n - down_sum / down_n
    p_value = float((1 + np.sum(null >= observed)) / (N_NULL + 1))
    return observed, p_value, float(np.mean(null))


def gene_tests(
    target_id: str,
    mapped_modules: dict[str, pd.DataFrame],
    target_universe: pd.DataFrame,
    module_ids: list[str],
) -> pd.DataFrame:
    rows = []
    for index, module_id in enumerate(module_ids):
        mapped = mapped_modules[module_id]
        up = mapped[mapped["direction"] == "up"]
        down = mapped[mapped["direction"] == "down"]
        enough = len(up) >= 10 and len(down) >= 10
        if enough:
            rng = np.random.default_rng(RANDOM_SEED + index + (0 if target_id == "DPN_vs_control" else 1000))
            concordance, empirical_p, null_mean = matched_null_p(mapped, target_universe, rng)
            rho, rho_p = stats.spearmanr(
                mapped["source_log2fc"], mapped["target_log2fc_oriented"], nan_policy="omit"
            )
        else:
            concordance = empirical_p = null_mean = rho = rho_p = np.nan
        rows.append(
            {
                "target_id": target_id,
                "module_id": module_id,
                "module_display": DISPLAY[module_id],
                "family": family_for(target_id, module_id),
                "mapped_n": int(len(mapped)),
                "mapped_up_n": int(len(up)),
                "mapped_down_n": int(len(down)),
                "minimum_10_per_direction": bool(enough),
                "gene_concordance": concordance,
                "matched_null_mean": null_mean,
                "matched_null_one_sided_p": empirical_p,
                "spearman_source_target_log2fc": float(rho) if np.isfinite(rho) else np.nan,
                "spearman_two_sided_p": float(rho_p) if np.isfinite(rho_p) else np.nan,
                "median_target_log2fc_source_up": float(up["target_log2fc_oriented"].median()) if len(up) else np.nan,
                "median_target_log2fc_source_down": float(down["target_log2fc_oriented"].median()) if len(down) else np.nan,
                "fraction_directionally_concordant_up": float((up["target_log2fc_oriented"] > 0).mean()) if len(up) else np.nan,
                "fraction_directionally_concordant_down": float((down["target_log2fc_oriented"] < 0).mean()) if len(down) else np.nan,
            }
        )
    result = pd.DataFrame(rows)
    result["matched_null_bh_q"] = result.groupby("family", group_keys=False)[
        "matched_null_one_sided_p"
    ].apply(bh_adjust)
    return result


def family_for(target_id: str, module_id: str) -> str:
    if target_id == "DPN_vs_control":
        return "DPN_neuronal_primary" if module_id in PRIMARY_DPN else "DPN_contextual"
    return "severity_primary" if module_id in PRIMARY_SEVERITY else "severity_contextual"


def hedges_g(positive: np.ndarray, negative: np.ndarray) -> float:
    n1 = len(positive)
    n0 = len(negative)
    if n1 < 2 or n0 < 2:
        return np.nan
    pooled_var = (
        (n1 - 1) * np.var(positive, ddof=1) + (n0 - 1) * np.var(negative, ddof=1)
    ) / (n1 + n0 - 2)
    if pooled_var <= 0:
        return np.nan
    d_value = (np.mean(positive) - np.mean(negative)) / math.sqrt(pooled_var)
    correction = 1 - 3 / (4 * (n1 + n0 - 2) - 1)
    return float(correction * d_value)


def exact_label_p(scores: np.ndarray, positive_n: int) -> tuple[float, int]:
    observed = float(scores[:positive_n].mean() - scores[positive_n:].mean())
    indices = np.arange(len(scores))
    total_sum = float(scores.sum())
    extreme = 0
    assignments = 0
    for positive_indices in itertools.combinations(indices, positive_n):
        positive_sum = float(scores[list(positive_indices)].sum())
        negative_sum = total_sum - positive_sum
        difference = positive_sum / positive_n - negative_sum / (len(scores) - positive_n)
        extreme += int(difference >= observed - 1e-15)
        assignments += 1
    return float(extreme / assignments), assignments


def score_modules(
    target_id: str,
    expression: pd.DataFrame,
    sample_groups: dict[str, str],
    positive_group: str,
    negative_group: str,
    mapped_modules: dict[str, pd.DataFrame],
    module_ids: list[str],
    expression_source: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = list(sample_groups)
    expression = expression[samples]
    ranks = expression.rank(axis=0, method="average", pct=True) - 0.5
    score_rows = []
    test_rows = []
    for module_id in module_ids:
        mapped = mapped_modules[module_id]
        up_ids = [gene for gene in mapped.loc[mapped["direction"] == "up", "human_gene_id"] if gene in ranks.index]
        down_ids = [gene for gene in mapped.loc[mapped["direction"] == "down", "human_gene_id"] if gene in ranks.index]
        enough = len(up_ids) >= 10 and len(down_ids) >= 10
        if enough:
            score = ranks.loc[up_ids].mean(axis=0) - ranks.loc[down_ids].mean(axis=0)
        else:
            score = pd.Series(np.nan, index=samples)
        for sample in samples:
            score_rows.append(
                {
                    "target_id": target_id,
                    "expression_source": expression_source,
                    "module_id": module_id,
                    "sample_id": sample,
                    "group": sample_groups[sample],
                    "score": float(score[sample]) if np.isfinite(score[sample]) else np.nan,
                }
            )
        positive_samples = [s for s in samples if sample_groups[s] == positive_group]
        negative_samples = [s for s in samples if sample_groups[s] == negative_group]
        ordered = positive_samples + negative_samples
        values = score.loc[ordered].to_numpy(dtype=float)
        if enough and np.isfinite(values).all():
            positive_values = score.loc[positive_samples].to_numpy(dtype=float)
            negative_values = score.loc[negative_samples].to_numpy(dtype=float)
            difference = float(positive_values.mean() - negative_values.mean())
            effect = hedges_g(positive_values, negative_values)
            exact_p, assignments = exact_label_p(values, len(positive_samples))
            loo_differences = []
            for omitted in ordered:
                retained_positive = [s for s in positive_samples if s != omitted]
                retained_negative = [s for s in negative_samples if s != omitted]
                loo_differences.append(
                    float(score.loc[retained_positive].mean() - score.loc[retained_negative].mean())
                )
            loo_min = float(min(loo_differences))
            loo_all_positive = bool(all(value > 0 for value in loo_differences))
        else:
            difference = effect = exact_p = loo_min = np.nan
            assignments = 0
            loo_all_positive = False
        test_rows.append(
            {
                "target_id": target_id,
                "expression_source": expression_source,
                "module_id": module_id,
                "module_display": DISPLAY[module_id],
                "family": family_for(target_id, module_id),
                "positive_group": positive_group,
                "negative_group": negative_group,
                "mapped_up_n": len(up_ids),
                "mapped_down_n": len(down_ids),
                "minimum_10_per_direction": bool(enough),
                "score_mean_difference": difference,
                "hedges_g": effect,
                "exact_one_sided_p": exact_p,
                "exact_assignments": assignments,
                "loo_min_difference": loo_min,
                "loo_all_positive": loo_all_positive,
            }
        )
    scores = pd.DataFrame(score_rows)
    tests = pd.DataFrame(test_rows)
    tests["exact_bh_q"] = tests.groupby("family", group_keys=False)["exact_one_sided_p"].apply(
        bh_adjust
    )
    return scores, tests


def build_gate_table(gene: pd.DataFrame, sample: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    merged = gene.merge(
        sample[sample["expression_source"] == "DESeq2_normalized_counts"],
        on=["target_id", "module_id", "module_display", "family"],
        suffixes=("_gene", "_sample"),
    )
    merged["component_pass"] = (
        merged["minimum_10_per_direction_gene"]
        & merged["minimum_10_per_direction_sample"]
        & (merged["gene_concordance"] > 0)
        & (merged["matched_null_bh_q"] < 0.10)
        & (merged["score_mean_difference"] > 0)
        & (merged["exact_bh_q"] < 0.10)
        & merged["loo_all_positive"]
    )
    distal_rows = merged[
        (merged["target_id"] == "DPN_vs_control") & merged["module_id"].isin(PRIMARY_DPN)
    ]
    severity_rows = merged[
        (merged["target_id"] == "severe_vs_moderate_axonal_loss")
        & merged["module_id"].isin(PRIMARY_SEVERITY)
    ]
    distal_pass = bool(distal_rows["component_pass"].any())
    severity_pass = bool(severity_rows["component_pass"].any())
    upgrade_pass = bool(distal_pass and severity_pass)
    summary = {
        "distal_neuronal_state_gate": "PASS" if distal_pass else "FAIL",
        "within_DPN_axonal_severity_gate": "PASS" if severity_pass else "FAIL",
        "pure_computational_manuscript_upgrade_gate": "PASS" if upgrade_pass else "FAIL",
        "distal_passing_components": distal_rows.loc[
            distal_rows["component_pass"], "module_id"
        ].tolist(),
        "severity_passing_components": severity_rows.loc[
            severity_rows["component_pass"], "module_id"
        ].tolist(),
    }
    return merged, summary


def make_figure(gene: pd.DataFrame, sample: pd.DataFrame, output: Path) -> None:
    primary_sample = sample[sample["expression_source"] == "DESeq2_normalized_counts"].copy()
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5))
    panels = [
        ("DPN_vs_control", "gene_concordance", gene, "Gene concordance", axes[0, 0]),
        ("DPN_vs_control", "hedges_g", primary_sample, "Sample score Hedges g", axes[0, 1]),
        (
            "severe_vs_moderate_axonal_loss",
            "gene_concordance",
            gene,
            "Gene concordance",
            axes[1, 0],
        ),
        (
            "severe_vs_moderate_axonal_loss",
            "hedges_g",
            primary_sample,
            "Sample score Hedges g",
            axes[1, 1],
        ),
    ]
    for panel_index, (target_id, metric, frame, x_label, axis) in enumerate(panels):
        block = frame[frame["target_id"] == target_id].copy()
        block = block.iloc[::-1].reset_index(drop=True)
        colors = [
            "#007C91"
            if family in {"DPN_neuronal_primary", "severity_primary"}
            else "#9AA0A6"
            for family in block["family"]
        ]
        axis.axvline(0, color="#333333", linewidth=0.8, linestyle="--")
        axis.scatter(block[metric], np.arange(len(block)), c=colors, s=55, zorder=3)
        q_column = "matched_null_bh_q" if metric == "gene_concordance" else "exact_bh_q"
        for y, row in block.iterrows():
            value = row[metric]
            if np.isfinite(value):
                axis.annotate(
                    f"q={row[q_column]:.3g}",
                    (value, y),
                    xytext=(6, 0),
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    fontsize=7.5,
                )
        axis.set_yticks(np.arange(len(block)))
        axis.set_yticklabels(block["module_display"], fontsize=8)
        axis.set_xlabel(x_label)
        axis.margins(x=0.24)
        axis.grid(axis="x", alpha=0.2)
        prefix = "DPN vs control" if target_id == "DPN_vs_control" else "Severe vs moderate axonal loss"
        axis.set_title(f"{'ABCD'[panel_index]}. {prefix}: {x_label.lower()}", loc="left", fontweight="bold")
    fig.suptitle(
        "Transport of frozen human hDRG components into independent human sural nerve",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.012,
        "Teal = frozen primary component; grey = contextual anchor. Positive values support source-defined direction.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.975))
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    for path in FILES.values():
        if not path.exists():
            raise FileNotFoundError(path)

    lookups = build_symbol_lookup(load_gene_info(NCBI / "Homo_sapiens.gene_info.gz"))
    modules, component_membership, source_audit = build_source_modules()

    dpn_groups = {
        **{f"C{i}": "Control" for i in range(1, 7)},
        **{sample: "DPN" for sample in ["S38", "S39", "S40", "S42", "S43", "S44"]},
    }
    moderate = {"S8", "S10", "S20", "S33"}
    severity_samples = [
        "S8", "S10", "S11", "S12", "S13", "S14", "S18", "S19", "S20",
        "S21", "S22", "S25", "S28", "S30", "S31", "S33", "S34",
    ]
    severity_groups = {
        sample: ("Moderate" if sample in moderate else "Severe") for sample in severity_samples
    }

    dpn_universe, dpn_expression, dpn_metadata, dpn_audit = load_target(
        "DPN_vs_control",
        FILES["dpn_control"],
        "deseq2_results",
        dpn_groups,
        1,
        lookups,
    )
    severity_universe, severity_expression, severity_metadata, severity_audit = load_target(
        "severe_vs_moderate_axonal_loss",
        FILES["axon_severity"],
        "Deseq_mod_vs_severe_axonLoss",
        severity_groups,
        -1,
        lookups,
    )

    dpn_mapped, dpn_mapping, dpn_coverage = map_modules(
        modules, dpn_universe, "DPN_vs_control", DPN_MODULES, lookups
    )
    severity_mapped, severity_mapping, severity_coverage = map_modules(
        modules,
        severity_universe,
        "severe_vs_moderate_axonal_loss",
        SEVERITY_MODULES,
        lookups,
    )

    dpn_gene = gene_tests("DPN_vs_control", dpn_mapped, dpn_universe, DPN_MODULES)
    severity_gene = gene_tests(
        "severe_vs_moderate_axonal_loss",
        severity_mapped,
        severity_universe,
        SEVERITY_MODULES,
    )
    gene = pd.concat([dpn_gene, severity_gene], ignore_index=True)

    dpn_scores, dpn_score_tests = score_modules(
        "DPN_vs_control",
        dpn_expression,
        dpn_groups,
        "DPN",
        "Control",
        dpn_mapped,
        DPN_MODULES,
        "DESeq2_normalized_counts",
    )
    severity_scores, severity_score_tests = score_modules(
        "severe_vs_moderate_axonal_loss",
        severity_expression,
        severity_groups,
        "Severe",
        "Moderate",
        severity_mapped,
        SEVERITY_MODULES,
        "DESeq2_normalized_counts",
    )

    tpm_expression, tpm_audit = load_tpm_sensitivity(FILES["quantile_tpm"], dpn_groups, lookups)
    common_tpm_ids = dpn_universe["human_gene_id"].isin(tpm_expression.index)
    tpm_universe = dpn_universe[common_tpm_ids].copy()
    tpm_expression = tpm_expression.loc[tpm_universe["human_gene_id"]]
    tpm_mapped, _, _ = map_modules(
        modules, tpm_universe, "DPN_vs_control", DPN_MODULES, lookups
    )
    tpm_scores, tpm_score_tests = score_modules(
        "DPN_vs_control",
        tpm_expression,
        dpn_groups,
        "DPN",
        "Control",
        tpm_mapped,
        DPN_MODULES,
        "quantile_normalized_TPM_sensitivity",
    )

    scores = pd.concat([dpn_scores, severity_scores, tpm_scores], ignore_index=True)
    score_tests = pd.concat(
        [dpn_score_tests, severity_score_tests, tpm_score_tests], ignore_index=True
    )
    gate_table, gate_summary = build_gate_table(gene, score_tests)

    component_membership.to_csv(
        TABLES / f"hDRG_source_defined_transport_components_{DATE}.tsv", sep="\t", index=False
    )
    pd.concat([dpn_metadata, severity_metadata], ignore_index=True).to_csv(
        METADATA / f"JCI184075_sural_sample_metadata_{DATE}.tsv", sep="\t", index=False
    )
    pd.concat([dpn_mapping, severity_mapping], ignore_index=True).to_csv(
        TABLES / f"JCI184075_hDRG_component_mapping_{DATE}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.concat([dpn_coverage, severity_coverage], ignore_index=True).to_csv(
        TABLES / f"JCI184075_hDRG_component_coverage_{DATE}.tsv", sep="\t", index=False
    )
    gene.to_csv(
        TABLES / f"JCI184075_hDRG_component_gene_tests_{DATE}.tsv", sep="\t", index=False
    )
    scores.to_csv(
        TABLES / f"JCI184075_hDRG_component_sample_scores_{DATE}.tsv", sep="\t", index=False
    )
    score_tests.to_csv(
        TABLES / f"JCI184075_hDRG_component_sample_tests_{DATE}.tsv", sep="\t", index=False
    )
    gate_table.to_csv(
        TABLES / f"JCI184075_hDRG_component_gate_table_{DATE}.tsv", sep="\t", index=False
    )

    qc = {
        "analysis_date": DATE,
        "random_seed": RANDOM_SEED,
        "matched_null_iterations": N_NULL,
        "source_audit": source_audit,
        "dpn_control_target_audit": dpn_audit,
        "axon_severity_target_audit": severity_audit,
        "tpm_sensitivity_audit": tpm_audit,
        "gate_summary": gate_summary,
        "interpretation": (
            "peer_reviewed_human_distal_nerve_upgrade"
            if gate_summary["pure_computational_manuscript_upgrade_gate"] == "PASS"
            else "transportability_boundary_or_partial_support"
        ),
    }
    with (TABLES / f"JCI184075_hDRG_component_validation_qc_{DATE}.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)

    make_figure(
        gene,
        score_tests,
        FIGURES / f"JCI184075_hDRG_component_transportability_{DATE}",
    )

    print(json.dumps(gate_summary, ensure_ascii=False, indent=2))
    print("\nPrimary gate rows:")
    print(
        gate_table[
            gate_table["family"].isin(["DPN_neuronal_primary", "severity_primary"])
        ][
            [
                "target_id",
                "module_id",
                "gene_concordance",
                "matched_null_bh_q",
                "score_mean_difference",
                "hedges_g",
                "exact_bh_q",
                "loo_all_positive",
                "component_pass",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
