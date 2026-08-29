#!/usr/bin/env python3
"""Extract audited major-cell-type proportion tests from the source hDRG supplement."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
RAW = PHASE / "data" / "raw" / "human_hDRG_preprint"
TABLES = PHASE / "results" / "tables"
DATE = "2026-08-27"
WORKBOOK = RAW / "700028_file12.xlsx"
XML = RAW / "PMC12871725_fullText.xml"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    raw = pd.read_excel(WORKBOOK, sheet_name="Sup_Data1", header=None, engine="openpyxl")
    header = [str(value).strip() if pd.notna(value) else f"unnamed_{index}" for index, value in enumerate(raw.iloc[0])]
    frame = raw.iloc[3:].copy()
    frame.columns = header
    frame = frame.rename(columns={header[0]: "row_id", header[1]: "cell_type"})
    keep = [
        "row_id", "cell_type", "group1", "group2", "n1", "n2", "statistic", "p", "p.adj", "p.adj.signif"
    ]
    frame = frame[keep].copy()
    for column in ["row_id", "n1", "n2", "statistic", "p", "p.adj"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame = frame[frame["cell_type"].notna()].copy()
    frame["pairwise_q_below_0_05"] = frame["p.adj"] < 0.05
    narrative = {
        "Neurons": "decreased_in_DPN",
        "Fibroblast": "increased_in_DPN",
        "SMC": "increased_in_DPN",
    }
    frame["source_article_narrative_for_DPN"] = frame["cell_type"].map(narrative)
    frame["direction_source"] = frame["source_article_narrative_for_DPN"].notna().map(
        {True: "article_results_text", False: "not_assigned"}
    )
    frame.to_csv(
        TABLES / f"hDRG_major_celltype_pairwise_proportion_tests_{DATE}.tsv",
        sep="\t",
        index=False,
    )
    diabetes_to_dpn = frame[(frame["group1"] == "Diabetic") & (frame["group2"] == "DPN")]
    summary = {
        "analysis_date": DATE,
        "workbook_sha256": sha256(WORKBOOK),
        "xml_sha256": sha256(XML),
        "sheet": "Sup_Data1",
        "extracted_pairwise_rows": int(len(frame)),
        "diabetes_to_DPN_significant_cell_types_q_below_0_05": diabetes_to_dpn.loc[
            diabetes_to_dpn["pairwise_q_below_0_05"], "cell_type"
        ].tolist(),
        "article_results_text": (
            "The source preprint states that DPN showed increased fibroblast and smooth-muscle-cell "
            "proportions and decreased neuronal proportion. The direction labels in the output are "
            "taken from that narrative, not inferred from the sign of the deposited pairwise statistic."
        ),
        "evidence_status": "source_preprint_context_not_independent_validation",
    }
    with (TABLES / f"hDRG_major_celltype_proportion_context_qc_{DATE}.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(diabetes_to_dpn.to_string(index=False))


if __name__ == "__main__":
    main()
