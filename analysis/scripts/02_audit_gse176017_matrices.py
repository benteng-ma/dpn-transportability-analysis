#!/usr/bin/env python3
"""Audit GSE176017 per-animal single-cell count matrices and emit pseudobulk counts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
RAW_DIR = PHASE / "data" / "raw" / "GSE176017"
MATRIX_DIR = RAW_DIR / "extracted"
SOFT = RAW_DIR / "GSE176017_family.soft.gz"
RESULTS = PHASE / "results" / "tables"
METADATA = PHASE / "metadata"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_soft() -> dict[str, dict]:
    samples: dict[str, dict] = {}
    current: str | None = None
    with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line.startswith("^SAMPLE = "):
                current = line.split("=", 1)[1].strip()
                samples[current] = {"gsm": current, "characteristics": {}}
            elif current and line.startswith("!Sample_title = "):
                samples[current]["title"] = line.split("=", 1)[1].strip()
            elif current and line.startswith("!Sample_source_name_ch1 = "):
                samples[current]["source"] = line.split("=", 1)[1].strip()
            elif current and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split("=", 1)[1].strip()
                if ":" in value:
                    key, entry = value.split(":", 1)
                    samples[current]["characteristics"][key.strip().lower()] = entry.strip()
    return samples


def group_from_title(title: str, disease_state: str) -> str:
    if title.startswith("N") or disease_state.lower() == "normal":
        return "Normal"
    if title.startswith("DM") or disease_state.lower() == "dm":
        return "Diabetes_no_allodynia"
    if title.startswith("PDPN") or disease_state.lower() == "pdpn":
        return "Painful_DPN"
    raise ValueError(f"Cannot assign group for title={title!r}, disease_state={disease_state!r}")


def audit_matrix(path: Path) -> tuple[dict, dict[str, int]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
        if not header or header == [""]:
            raise RuntimeError(f"Empty header in {path}")
        n_cells = len(header)
        total_counts = np.zeros(n_cells, dtype=np.int64)
        detected_genes = np.zeros(n_cells, dtype=np.int64)
        gene_totals: dict[str, int] = {}
        gene_seen: Counter[str] = Counter()
        malformed_rows = 0
        negative_rows = 0
        noninteger_rows = 0
        n_gene_rows = 0

        for line_number, raw_line in enumerate(handle, start=2):
            line = raw_line.rstrip("\n\r")
            first_tab = line.find("\t")
            if first_tab < 1:
                malformed_rows += 1
                continue
            gene = line[:first_tab]
            values = np.fromstring(line[first_tab + 1 :], sep="\t", dtype=np.float64)
            if values.size != n_cells:
                malformed_rows += 1
                continue
            if np.any(values < 0):
                negative_rows += 1
            if not np.all(values == np.floor(values)):
                noninteger_rows += 1
            counts = values.astype(np.int64)
            total_counts += counts
            detected_genes += counts > 0
            gene_seen[gene] += 1
            gene_totals[gene] = gene_totals.get(gene, 0) + int(counts.sum())
            n_gene_rows += 1

    duplicate_gene_symbols = sum(count - 1 for count in gene_seen.values() if count > 1)
    record = {
        "file": path.name,
        "file_size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "n_cells": n_cells,
        "n_gene_rows": n_gene_rows,
        "n_unique_gene_symbols": len(gene_seen),
        "duplicate_gene_symbol_rows": duplicate_gene_symbols,
        "malformed_rows": malformed_rows,
        "negative_value_rows": negative_rows,
        "noninteger_value_rows": noninteger_rows,
        "total_umi": int(total_counts.sum()),
        "median_umi_per_cell": float(np.median(total_counts)),
        "minimum_umi_per_cell": int(total_counts.min()),
        "maximum_umi_per_cell": int(total_counts.max()),
        "median_detected_genes_per_cell": float(np.median(detected_genes)),
        "minimum_detected_genes_per_cell": int(detected_genes.min()),
        "maximum_detected_genes_per_cell": int(detected_genes.max()),
        "zero_library_cells": int((total_counts == 0).sum()),
        "duplicate_cell_barcodes": n_cells - len(set(header)),
    }
    return record, gene_totals


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    soft_samples = parse_soft()
    matrices = sorted(MATRIX_DIR.glob("*.expression_matrix.txt.gz"))
    if len(matrices) != 8:
        raise RuntimeError(f"Expected 8 per-animal matrices, found {len(matrices)}")

    audit_rows: list[dict] = []
    metadata_rows: list[dict] = []
    pseudobulk_by_sample: dict[str, dict[str, int]] = {}
    for path in matrices:
        match = re.match(r"(GSM\d+)_([^.]*)\.expression_matrix\.txt\.gz$", path.name)
        if not match:
            raise RuntimeError(f"Unexpected matrix filename: {path.name}")
        gsm, filename_title = match.groups()
        meta = soft_samples.get(gsm)
        if meta is None:
            raise RuntimeError(f"{gsm} absent from SOFT metadata")
        title = meta.get("title", filename_title)
        if title != filename_title:
            raise RuntimeError(f"Title mismatch for {gsm}: filename={filename_title}, SOFT={title}")
        disease_state = meta["characteristics"].get("disease state", title.rstrip("1234567890"))
        group = group_from_title(title, disease_state)
        record, gene_totals = audit_matrix(path)
        record.update({"gsm": gsm, "sample_id": title, "group": group})
        audit_rows.append(record)
        pseudobulk_by_sample[title] = gene_totals
        metadata_rows.append(
            {
                "gsm": gsm,
                "sample_id": title,
                "group": group,
                "disease_state_geo": disease_state,
                "treatment": meta["characteristics"].get("treatment", ""),
                "strain": meta["characteristics"].get("strain", ""),
                "source": meta.get("source", ""),
                "biological_unit": "one rat",
                "inferential_role": "animal-level progression replication",
            }
        )

    audit = pd.DataFrame(audit_rows).sort_values("sample_id", kind="stable")
    metadata = pd.DataFrame(metadata_rows).sort_values("sample_id", kind="stable")
    group_counts = metadata["group"].value_counts().to_dict()
    expected_counts = {
        "Normal": 2,
        "Diabetes_no_allodynia": 2,
        "Painful_DPN": 4,
    }

    all_genes = sorted({gene for sample in pseudobulk_by_sample.values() for gene in sample})
    pseudobulk = pd.DataFrame(
        {
            sample: [counts.get(gene, 0) for gene in all_genes]
            for sample, counts in pseudobulk_by_sample.items()
        },
        index=all_genes,
    )
    pseudobulk.index.name = "gene_symbol"

    output_audit = RESULTS / "GSE176017_per_animal_matrix_qc_2026-08-27.tsv"
    output_meta = METADATA / "GSE176017_animal_metadata_2026-08-27.tsv"
    output_pseudobulk = RESULTS / "GSE176017_animal_pseudobulk_raw_counts_2026-08-27.tsv.gz"
    output_qc = RESULTS / "GSE176017_audit_summary_2026-08-27.json"
    audit.to_csv(output_audit, sep="\t", index=False)
    metadata.to_csv(output_meta, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    pseudobulk.to_csv(output_pseudobulk, sep="\t", compression="gzip")

    passed = bool(
        group_counts == expected_counts
        and (audit["malformed_rows"] == 0).all()
        and (audit["negative_value_rows"] == 0).all()
        and (audit["noninteger_value_rows"] == 0).all()
        and (audit["zero_library_cells"] == 0).all()
        and (audit["duplicate_cell_barcodes"] == 0).all()
    )
    qc = {
        "status": "PASS" if passed else "FAIL",
        "soft_file": str(SOFT),
        "soft_sha256": sha256(SOFT),
        "n_animal_matrices": len(matrices),
        "group_counts": group_counts,
        "total_cells": int(audit["n_cells"].sum()),
        "union_gene_symbols": len(all_genes),
        "total_umi": int(audit["total_umi"].sum()),
        "minimum_cells_in_an_animal": int(audit["n_cells"].min()),
        "maximum_cells_in_an_animal": int(audit["n_cells"].max()),
        "all_matrices_structurally_valid": bool(
            (audit[["malformed_rows", "negative_value_rows", "noninteger_value_rows"]] == 0).all().all()
        ),
        "animal_is_inferential_unit": True,
        "small_group_warning": (
            "Normal and diabetes-without-allodynia each have n=2 animals. Effect size, exact uncertainty, "
            "and leave-one-animal-out behavior must be shown; this dataset cannot provide definitive inference alone."
        ),
        "outputs": [str(output_audit), str(output_meta), str(output_pseudobulk)],
    }
    output_qc.write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(qc, indent=2, ensure_ascii=False))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

