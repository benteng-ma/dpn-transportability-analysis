#!/usr/bin/env python3
"""Audit hDRG-signature coverage in PXD062366 without testing group outcomes."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
RAW = PHASE / "data" / "raw" / "PXD062366" / "Proteins_Report.csv"
GENE_INFO = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27" / "Homo_sapiens.gene_info.gz"
SIGNATURES = PHASE / "results" / "tables" / "hDRG_frozen_primary_stage_signatures_2026-08-27.tsv"
TABLES = PHASE / "results" / "tables"
DATE = "2026-08-27"


def build_lookup(info: pd.DataFrame) -> tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]:
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


def resolve(symbol: str, lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[str | None, str]:
    exact, folded, synonyms, _ = lookup
    text = symbol.strip()
    if text in exact:
        return exact[text], "official_exact"
    key = text.upper()
    if key in folded:
        return folded[key], "official_casefold"
    if key in synonyms:
        return synonyms[key], "unique_synonym"
    return None, "unresolved"


def resolve_protein_group(value: object, lookup: tuple[dict[str, str], dict[str, str], dict[str, str], dict[str, str]]) -> tuple[str | None, str, str]:
    if pd.isna(value):
        return None, "missing_gene_annotation", ""
    symbols = [item.strip() for item in str(value).split(";") if item.strip()]
    resolved = [resolve(symbol, lookup) for symbol in symbols]
    gene_ids = sorted({item[0] for item in resolved if item[0] is not None}, key=lambda value: int(value))
    if len(gene_ids) == 1:
        methods = sorted({item[1] for item in resolved if item[0] is not None})
        return gene_ids[0], "+".join(methods), ";".join(symbols)
    if len(gene_ids) > 1:
        return None, "ambiguous_multi_gene_protein_group", ";".join(symbols)
    return None, "unresolved_protein_group", ";".join(symbols)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(RAW, na_values=["NaN", ""])
    sample_columns = list(raw.columns[3:])
    values = raw[sample_columns].apply(pd.to_numeric, errors="coerce")
    raw["detected_sample_n"] = values.notna().sum(axis=1)
    raw["detection_fraction"] = raw["detected_sample_n"] / len(sample_columns)
    raw["median_log2_intensity_observed"] = np.log2(values.where(values > 0)).median(axis=1)

    info = pd.read_csv(GENE_INFO, sep="\t", compression="gzip", dtype=str, na_filter=False)
    lookup = build_lookup(info)
    resolved = raw["PG.Genes"].map(lambda value: resolve_protein_group(value, lookup))
    raw["human_gene_id"] = [item[0] for item in resolved]
    raw["gene_mapping_method"] = [item[1] for item in resolved]
    raw["parsed_gene_symbols"] = [item[2] for item in resolved]
    raw["current_symbol"] = raw["human_gene_id"].map(lookup[3])
    raw["source_row"] = np.arange(2, len(raw) + 2)

    resolved_rows = raw[raw["human_gene_id"].notna()].copy()
    resolved_rows = resolved_rows.sort_values(
        ["human_gene_id", "detected_sample_n", "median_log2_intensity_observed", "PG.UniProtIds", "source_row"],
        ascending=[True, False, False, True, True],
    )
    resolved_rows["selected_gene_representative"] = ~resolved_rows.duplicated("human_gene_id", keep="first")
    representatives = resolved_rows[resolved_rows["selected_gene_representative"]].copy()
    detection_minimum = math.ceil(0.50 * len(sample_columns))
    representatives["passes_detection_50pct"] = representatives["detected_sample_n"] >= detection_minimum

    signatures = pd.read_csv(SIGNATURES, sep="\t")
    signature_resolution = signatures["gene"].map(lambda value: resolve(str(value), lookup))
    signatures["human_gene_id"] = [item[0] for item in signature_resolution]
    signatures["signature_mapping_method"] = [item[1] for item in signature_resolution]
    signatures = signatures.sort_values(["contrast_id", "direction", "p_val_adj", "gene"])
    signatures["duplicate_resolved_signature_gene"] = signatures.duplicated(
        ["contrast_id", "direction", "human_gene_id"], keep="first"
    ) & signatures["human_gene_id"].notna()
    unique_signatures = signatures[
        signatures["human_gene_id"].notna() & ~signatures["duplicate_resolved_signature_gene"]
    ]
    available_all = set(representatives["human_gene_id"])
    available_50 = set(representatives.loc[representatives["passes_detection_50pct"], "human_gene_id"])
    coverage_rows = []
    for (contrast, direction), frame in signatures.groupby(["contrast_id", "direction"]):
        unique = unique_signatures[
            (unique_signatures["contrast_id"] == contrast) & (unique_signatures["direction"] == direction)
        ]
        coverage_rows.append(
            {
                "contrast_id": contrast,
                "direction": direction,
                "original_gene_n": int(len(frame)),
                "resolved_unique_gene_n": int(unique["human_gene_id"].nunique()),
                "tear_proteome_gene_n_any_detection": int(unique["human_gene_id"].isin(available_all).sum()),
                "tear_proteome_gene_n_detection_ge_50pct": int(unique["human_gene_id"].isin(available_50).sum()),
                "fraction_original_detection_ge_50pct": float(unique["human_gene_id"].isin(available_50).sum() / len(frame)),
            }
        )
    coverage = pd.DataFrame(coverage_rows)

    output_mapping = TABLES / f"PXD062366_protein_to_human_gene_mapping_{DATE}.tsv.gz"
    raw.to_csv(output_mapping, sep="\t", index=False, compression="gzip")
    representatives[[
        "PG.Genes",
        "PG.ProteinDescriptions",
        "PG.UniProtIds",
        "human_gene_id",
        "current_symbol",
        "gene_mapping_method",
        "detected_sample_n",
        "detection_fraction",
        "median_log2_intensity_observed",
        "selected_gene_representative",
        "passes_detection_50pct",
        "source_row",
    ]].to_csv(TABLES / f"PXD062366_gene_representatives_{DATE}.tsv", sep="\t", index=False)
    signatures.to_csv(TABLES / f"PXD062366_hDRG_signature_mapping_{DATE}.tsv", sep="\t", index=False)
    coverage.to_csv(TABLES / f"PXD062366_hDRG_signature_coverage_{DATE}.tsv", sep="\t", index=False)

    audit = {
        "sample_count": len(sample_columns),
        "protein_group_rows": int(len(raw)),
        "unambiguous_gene_mapped_rows": int(raw["human_gene_id"].notna().sum()),
        "unique_gene_representatives": int(len(representatives)),
        "detection_50pct_minimum_samples": detection_minimum,
        "gene_representatives_passing_detection_50pct": int(representatives["passes_detection_50pct"].sum()),
        "mapping_method_counts": raw["gene_mapping_method"].value_counts(dropna=False).to_dict(),
        "analysis_boundary": "No Healthy/preDM/T2DM group association was calculated in this coverage audit.",
    }
    with (TABLES / f"PXD062366_hDRG_signature_coverage_audit_{DATE}.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    print(coverage.to_string(index=False))


if __name__ == "__main__":
    main()
