#!/usr/bin/env python3
"""Map subject IDs in the raw PXD062366 protein report to groups in Table S1.

The supplement relabels subjects as within-group ordinals. This script recovers the
mapping by exact matching of abundance vectors over shared UniProt identifiers.
Raw inputs are read only; all derived outputs are written under the phase results
and metadata directories.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
RAW = PHASE / "data" / "raw" / "PXD062366" / "Proteins_Report.csv"
SUPPLEMENT = PHASE / "data" / "raw" / "PXD062366" / "supplement" / "Table S1.xlsx"
RESULTS = PHASE / "results" / "tables"
METADATA = PHASE / "metadata"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def best_match(raw_values: pd.Series, supplement: pd.DataFrame) -> list[dict]:
    a = pd.to_numeric(raw_values, errors="coerce").to_numpy(dtype=float)
    candidates: list[dict] = []
    for column in supplement.columns:
        b = pd.to_numeric(supplement[column], errors="coerce").to_numpy(dtype=float)
        both = np.isfinite(a) & np.isfinite(b)
        union = np.isfinite(a) | np.isfinite(b)
        exact = np.isclose(a[both], b[both], rtol=1e-7, atol=1e-5)
        correlation = float(np.corrcoef(a[both], b[both])[0, 1]) if both.sum() >= 3 else np.nan
        candidates.append(
            {
                "supplement_column": column,
                "n_both_observed": int(both.sum()),
                "n_union_observed": int(union.sum()),
                "n_exact": int(exact.sum()),
                "exact_fraction_among_both": float(exact.mean()) if both.any() else 0.0,
                "missingness_agreement_fraction": float((np.isfinite(a) == np.isfinite(b)).mean()),
                "pearson_r": correlation,
            }
        )
    return sorted(
        candidates,
        key=lambda item: (
            item["exact_fraction_among_both"],
            item["missingness_agreement_fraction"],
            -1 if np.isnan(item["pearson_r"]) else item["pearson_r"],
            item["n_both_observed"],
        ),
        reverse=True,
    )


def group_for_index(index: int) -> tuple[str, int]:
    if 0 <= index <= 12:
        return "Healthy", index + 1
    if 13 <= index <= 28:
        return "preDM", index - 12
    if 29 <= index <= 46:
        return "T2DM", index - 28
    raise ValueError(f"Unexpected supplement sample index: {index}")


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(RAW)
    sample_columns = list(raw.columns[3:])
    supplement_raw = pd.read_excel(SUPPLEMENT, header=None)

    expected_shape = (1373, 50)
    if supplement_raw.shape != expected_shape:
        raise RuntimeError(f"Unexpected Table S1 shape {supplement_raw.shape}; expected {expected_shape}")

    supplement = supplement_raw.iloc[3:].copy()
    supplement.columns = ["Genes", "ProteinDescriptions", "UniProt"] + [f"S{i:02d}" for i in range(47)]
    merged = raw.merge(
        supplement,
        left_on="PG.UniProtIds",
        right_on="UniProt",
        how="inner",
        validate="one_to_one",
    )
    if len(merged) < 1000:
        raise RuntimeError(f"Only {len(merged)} shared one-to-one proteins; expected at least 1000")

    supplement_columns = [f"S{i:02d}" for i in range(47)]
    rows: list[dict] = []
    all_ranked: dict[str, list[dict]] = {}
    for sample in sample_columns:
        ranked = best_match(merged[sample], merged[supplement_columns])
        all_ranked[sample] = ranked
        best = ranked[0]
        second = ranked[1]
        supplement_index = int(best["supplement_column"][1:])
        group, within_group_ordinal = group_for_index(supplement_index)
        rows.append(
            {
                "sample_id": sample,
                "group": group,
                "supplement_within_group_ordinal": within_group_ordinal,
                "supplement_column_index_zero_based": supplement_index,
                **best,
                "second_best_exact_fraction": second["exact_fraction_among_both"],
                "second_best_missingness_agreement_fraction": second["missingness_agreement_fraction"],
                "exact_fraction_margin": best["exact_fraction_among_both"]
                - second["exact_fraction_among_both"],
            }
        )

    mapping = pd.DataFrame(rows).sort_values(
        ["group", "supplement_within_group_ordinal"], kind="stable"
    )
    unique_assignments = mapping["supplement_column"].nunique() == len(mapping)
    perfect_exact = bool((mapping["exact_fraction_among_both"] == 1.0).all())
    strong_missingness = bool((mapping["missingness_agreement_fraction"] >= 0.999).all())
    expected_groups = mapping["group"].value_counts().to_dict() == {
        "T2DM": 18,
        "preDM": 16,
        "Healthy": 13,
    }
    passed = bool(unique_assignments and perfect_exact and strong_missingness and expected_groups)

    output_mapping = RESULTS / "PXD062366_sample_group_mapping_2026-08-27.tsv"
    output_metadata = METADATA / "PXD062366_sample_metadata_2026-08-27.tsv"
    output_qc = RESULTS / "PXD062366_mapping_qc_2026-08-27.json"
    mapping.to_csv(output_mapping, sep="\t", index=False)
    mapping[["sample_id", "group", "supplement_within_group_ordinal"]].to_csv(
        output_metadata, sep="\t", index=False
    )

    qc = {
        "status": "PASS" if passed else "FAIL",
        "raw_file": str(RAW),
        "raw_sha256": sha256(RAW),
        "supplement_file": str(SUPPLEMENT),
        "supplement_sha256": sha256(SUPPLEMENT),
        "raw_shape": list(raw.shape),
        "supplement_shape": list(supplement_raw.shape),
        "shared_one_to_one_uniprot_rows": int(len(merged)),
        "n_sample_columns": len(sample_columns),
        "group_counts": mapping["group"].value_counts().sort_index().to_dict(),
        "unique_assignments": unique_assignments,
        "all_best_matches_exact": perfect_exact,
        "all_missingness_agreement_at_least_0_999": strong_missingness,
        "minimum_exact_fraction": float(mapping["exact_fraction_among_both"].min()),
        "minimum_missingness_agreement_fraction": float(mapping["missingness_agreement_fraction"].min()),
        "minimum_exact_match_margin": float(mapping["exact_fraction_margin"].min()),
        "outputs": [str(output_mapping), str(output_metadata)],
        "interpretation": (
            "Subject IDs in Proteins_Report.csv are deterministically assigned to the Healthy, preDM, "
            "and T2DM groups by exact abundance-vector matches to the author supplement."
        ),
    }
    output_qc.write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

