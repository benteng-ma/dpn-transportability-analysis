#!/usr/bin/env python3
"""Audit the GSE302658 painful diabetic neuropathy blood RNA-seq trial deposit.

This script intentionally performs no signature-outcome testing.  It resolves the
sample/subject/visit structure, checks the deposited Salmon matrix, and writes the
metadata and Ensembl-to-NCBI mapping needed for the subsequently frozen analysis.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
RAW = PHASE / "data" / "raw" / "human_PDN_trial_GSE302658"
NCBI = PHASE / "data" / "raw" / "NCBI_orthology_2026-08-27"
TABLES = PHASE / "results" / "tables"
METADATA = PHASE / "metadata"
DATE = "2026-08-27"

SERIES = RAW / "GSE302658_series_matrix.txt.gz"
COUNTS = RAW / "GSE302658_salmon.merged.transcript_counts.tsv.gz"
GENE_INFO = NCBI / "Homo_sapiens.gene_info.gz"

FIXED_DESCRIPTION_FIELDS = [
    "library_name",
    "subject_id",
    "visit",
    "randomized_treatment",
    "age_years",
    "sex",
    "weight_kg",
    "height_cm",
    "bmi",
    "race",
    "matrix_sample",
    "extraction_id",
    "study_day",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def snake_case(text: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return value or "unnamed"


def read_geo_metadata_rows(path: Path) -> dict[str, list[list[str]]]:
    rows: defaultdict[str, list[list[str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for raw_line in handle:
            if not raw_line.startswith("!"):
                continue
            parsed = next(csv.reader([raw_line.rstrip("\r\n")], delimiter="\t", quotechar='"'))
            rows[parsed[0]].append(parsed[1:])
    return dict(rows)


def parse_clinical_description(text: str) -> dict[str, object]:
    fields = next(csv.reader([text], delimiter=",", quotechar='"'))
    if len(fields) < len(FIXED_DESCRIPTION_FIELDS):
        raise RuntimeError(f"Clinical description has only {len(fields)} fields: {text[:160]}")
    if len(fields) < 55:
        fields.extend([""] * (55 - len(fields)))
    record: dict[str, object] = dict(zip(FIXED_DESCRIPTION_FIELDS, fields[:13]))
    # The deposit stores two named NRS pairs, followed by one unlabeled NPSI
    # assessment-day field, and then the named NPSI item/value pairs.
    leading_outcomes = fields[13:17]
    for index in range(0, len(leading_outcomes), 2):
        label = leading_outcomes[index].strip()
        value = leading_outcomes[index + 1].strip()
        if label:
            record[snake_case(label)] = value
    record["npsi_assessment_study_day"] = fields[17].strip()
    remainder = fields[18:]
    if len(remainder) % 2:
        remainder.append("")
    for index in range(0, len(remainder), 2):
        label = remainder[index].strip()
        value = remainder[index + 1].strip()
        if label:
            record[snake_case(label)] = value
    return record


def build_metadata(rows: dict[str, list[list[str]]]) -> pd.DataFrame:
    titles = rows["!Sample_title"][0]
    accessions = rows["!Sample_geo_accession"][0]
    description_rows = rows["!Sample_description"]
    if len(description_rows) != 2:
        raise RuntimeError(f"Expected two !Sample_description rows, found {len(description_rows)}")
    library_descriptions, clinical_descriptions = description_rows
    lengths = {len(titles), len(accessions), len(library_descriptions), len(clinical_descriptions)}
    if len(lengths) != 1:
        raise RuntimeError(f"GEO metadata row lengths disagree: {sorted(lengths)}")

    records: list[dict[str, object]] = []
    for title, accession, library_text, clinical_text in zip(
        titles, accessions, library_descriptions, clinical_descriptions
    ):
        record = parse_clinical_description(clinical_text)
        record["geo_accession"] = accession
        record["sample_title"] = title
        record["library_description"] = library_text.removeprefix("Library name: ")
        match = re.fullmatch(r"(E\d+)-Visit ([38])-(.+)", title)
        if not match:
            raise RuntimeError(f"Unrecognized sample title: {title}")
        record["title_subject_id"] = match.group(1)
        record["title_visit"] = f"Visit {match.group(2)}"
        record["title_treatment"] = match.group(3)
        records.append(record)

    frame = pd.DataFrame(records)
    numeric_columns = [
        "age_years",
        "weight_kg",
        "height_cm",
        "bmi",
        "study_day",
        "average_pain_nrs_last_12_hours",
        "worst_pain_nrs_last_12_hours",
        "npsi_total_score",
        "burning_pain",
        "burning_superfic_spont_pain_sub_score",
        "electric_shock_pain",
        "evoked_pain_sub_score",
        "pain_attacks",
        "pain_brushing",
        "pain_pins_and_needle",
        "pain_present",
        "pain_pressure",
        "pain_something_cold",
        "paresthesia_dysesthesia_sub_score",
        "paroxysmal_pain_sub_score",
        "pressing_deep_spont_pain_sub_score",
        "pressure_pain",
        "squeezing_pain",
        "stabbing_pain",
        "tingling",
    ]
    for column in numeric_columns:
        if column not in frame:
            frame[column] = np.nan
        frame[column] = pd.to_numeric(frame[column], errors="coerce")

    frame["npsi_assessment_study_day"] = pd.to_numeric(
        frame["npsi_assessment_study_day"], errors="coerce"
    )
    frame["race"] = frame["race"].str.strip()

    frame["metadata_title_agreement"] = (
        (frame["subject_id"] == frame["title_subject_id"])
        & (frame["visit"] == frame["title_visit"])
        & (frame["randomized_treatment"] == frame["title_treatment"])
    )
    frame["library_name_agreement"] = frame["library_name"] == frame["library_description"]
    frame["nominal_study_day"] = frame["visit"].map({"Visit 3": -1, "Visit 8": 28})
    frame["visit_day_compatible"] = (
        ((frame["visit"] == "Visit 3") & frame["study_day"].eq(-1))
        | ((frame["visit"] == "Visit 8") & frame["study_day"].between(26, 35))
    )
    preferred = [
        "geo_accession",
        "matrix_sample",
        "library_name",
        "subject_id",
        "visit",
        "study_day",
        "randomized_treatment",
        "age_years",
        "sex",
        "weight_kg",
        "height_cm",
        "bmi",
        "race",
        "average_pain_nrs_last_12_hours",
        "worst_pain_nrs_last_12_hours",
        "npsi_total_score",
    ]
    remaining = [column for column in frame.columns if column not in preferred]
    return frame[preferred + remaining]


def build_ensembl_mapping(path: Path) -> pd.DataFrame:
    info = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, na_filter=False)
    records: list[dict[str, str]] = []
    for row in info[["GeneID", "Symbol", "dbXrefs", "type_of_gene"]].itertuples(index=False):
        ensembl_ids = re.findall(r"(?:^|\|)Ensembl:(ENSG\d+)", row.dbXrefs)
        for ensembl_id in ensembl_ids:
            records.append(
                {
                    "ensembl_gene_id": ensembl_id,
                    "human_gene_id": row.GeneID,
                    "current_symbol": row.Symbol,
                    "type_of_gene": row.type_of_gene,
                }
            )
    mapping = pd.DataFrame(records).drop_duplicates()
    ambiguity = mapping.groupby("ensembl_gene_id")["human_gene_id"].nunique()
    mapping["ensembl_mapping_is_one_to_one"] = mapping["ensembl_gene_id"].map(ambiguity).eq(1)
    return mapping.sort_values(["ensembl_gene_id", "human_gene_id"])


def audit_count_matrix(path: Path, metadata: pd.DataFrame, mapping: pd.DataFrame) -> tuple[dict[str, object], set[str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    if header[:2] != ["tx", "gene_id"]:
        raise RuntimeError(f"Unexpected count matrix leading columns: {header[:2]}")
    matrix_samples = header[2:]
    deposited_gene_ids: set[str] = set()
    deposited_tx_ids: set[str] = set()
    row_count = 0
    duplicate_transcript_rows = 0
    fractional_value_seen = False
    first_numeric_columns = matrix_samples[: min(6, len(matrix_samples))]
    usecols = ["tx", "gene_id", *first_numeric_columns]
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", usecols=usecols, chunksize=25_000):
        row_count += len(chunk)
        duplicate_transcript_rows += int(chunk["tx"].isin(deposited_tx_ids).sum())
        deposited_tx_ids.update(chunk["tx"].astype(str))
        deposited_gene_ids.update(chunk["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True))
        if not fractional_value_seen:
            numeric = chunk[first_numeric_columns].to_numpy(dtype=float)
            fractional_value_seen = bool(np.any(np.abs(numeric - np.round(numeric)) > 1e-8))

    mapped_ids = set(mapping.loc[mapping["ensembl_mapping_is_one_to_one"], "ensembl_gene_id"])
    metadata_samples = metadata["matrix_sample"].tolist()
    audit = {
        "matrix_transcript_rows": row_count,
        "matrix_unique_transcript_ids": len(deposited_tx_ids),
        "duplicate_transcript_rows_detected_within_chunks_or_prior_chunks": duplicate_transcript_rows,
        "matrix_unique_ensembl_gene_ids": len(deposited_gene_ids),
        "matrix_sample_count": len(matrix_samples),
        "metadata_sample_count": len(metadata_samples),
        "matrix_metadata_sample_set_equal": set(matrix_samples) == set(metadata_samples),
        "matrix_metadata_sample_order_equal": matrix_samples == metadata_samples,
        "unmatched_matrix_samples": sorted(set(matrix_samples) - set(metadata_samples)),
        "unmatched_metadata_samples": sorted(set(metadata_samples) - set(matrix_samples)),
        "deposited_gene_ids_one_to_one_mapped_to_ncbi": len(deposited_gene_ids & mapped_ids),
        "deposited_gene_id_one_to_one_mapping_fraction": len(deposited_gene_ids & mapped_ids) / len(deposited_gene_ids),
        "fractional_salmon_estimated_counts_observed": fractional_value_seen,
        "matrix_sha256": sha256(path),
    }
    return audit, deposited_gene_ids


def summarize_metadata(frame: pd.DataFrame) -> dict[str, object]:
    duplicate_subject_visit = frame.duplicated(["subject_id", "visit"], keep=False)
    subjects_by_visit = frame.groupby("subject_id")["visit"].agg(lambda values: sorted(set(values)))
    paired_subjects = subjects_by_visit[subjects_by_visit.map(lambda values: values == ["Visit 3", "Visit 8"])].index
    paired = frame[frame["subject_id"].isin(paired_subjects)]
    paired_outcome_counts: dict[str, int] = {}
    for outcome in [
        "average_pain_nrs_last_12_hours",
        "worst_pain_nrs_last_12_hours",
        "npsi_total_score",
    ]:
        complete = paired.pivot(index="subject_id", columns="visit", values=outcome)
        paired_outcome_counts[outcome] = int(complete[["Visit 3", "Visit 8"]].notna().all(axis=1).sum())

    outcome_completeness = {}
    clinical_columns = [
        column
        for column in frame.columns
        if column.endswith("score") or column.endswith("pain") or column in {
            "average_pain_nrs_last_12_hours",
            "worst_pain_nrs_last_12_hours",
            "pain_attacks",
            "pain_brushing",
            "pain_pins_and_needle",
            "pain_present",
            "pain_pressure",
            "pain_something_cold",
            "pressure_pain",
            "squeezing_pain",
            "stabbing_pain",
            "tingling",
        }
    ]
    for outcome in sorted(set(clinical_columns)):
        outcome_completeness[outcome] = {
            str(visit): int(values.notna().sum())
            for visit, values in frame.groupby("visit")[outcome]
        }

    return {
        "sample_count": int(len(frame)),
        "unique_subject_count": int(frame["subject_id"].nunique()),
        "visit_counts": frame["visit"].value_counts().sort_index().to_dict(),
        "treatment_counts_by_sample": frame["randomized_treatment"].value_counts().to_dict(),
        "unique_subject_counts_by_treatment": frame.groupby("randomized_treatment")["subject_id"].nunique().to_dict(),
        "sex_counts_by_subject": frame.drop_duplicates("subject_id")["sex"].value_counts(dropna=False).to_dict(),
        "race_counts_by_subject": frame.drop_duplicates("subject_id")["race"].value_counts(dropna=False).to_dict(),
        "duplicate_subject_visit_row_count": int(duplicate_subject_visit.sum()),
        "complete_visit3_visit8_subject_count": int(len(paired_subjects)),
        "paired_outcome_counts": paired_outcome_counts,
        "all_title_metadata_fields_agree": bool(frame["metadata_title_agreement"].all()),
        "all_library_description_fields_agree": bool(frame["library_name_agreement"].all()),
        "all_visit_day_fields_compatible": bool(frame["visit_day_compatible"].all()),
        "actual_study_day_counts_by_visit": {
            str(visit): {str(day): int(count) for day, count in values.value_counts().sort_index().items()}
            for visit, values in frame.groupby("visit")["study_day"]
        },
        "baseline_age_range": [
            float(frame.loc[frame["visit"] == "Visit 3", "age_years"].min()),
            float(frame.loc[frame["visit"] == "Visit 3", "age_years"].max()),
        ],
        "baseline_bmi_range": [
            float(frame.loc[frame["visit"] == "Visit 3", "bmi"].min()),
            float(frame.loc[frame["visit"] == "Visit 3", "bmi"].max()),
        ],
        "outcome_nonmissing_counts_by_visit": outcome_completeness,
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    METADATA.mkdir(parents=True, exist_ok=True)
    rows = read_geo_metadata_rows(SERIES)
    metadata = build_metadata(rows)
    mapping = build_ensembl_mapping(GENE_INFO)
    matrix_audit, deposited_gene_ids = audit_count_matrix(COUNTS, metadata, mapping)
    mapping["present_in_GSE302658"] = mapping["ensembl_gene_id"].isin(deposited_gene_ids)

    metadata.to_csv(METADATA / f"GSE302658_clinical_sample_metadata_{DATE}.tsv", sep="\t", index=False)
    mapping.to_csv(
        TABLES / f"GSE302658_ensembl_to_ncbi_gene_mapping_{DATE}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    qc = {
        "dataset": "GSE302658",
        "audit_date": DATE,
        "design_interpretation": {
            "population": "painful diabetic polyneuropathy",
            "visit_3": "pretreatment; study day -1",
            "visit_8": "post-randomized treatment; study day 28",
            "randomized_arms": ["Placebo", "AZD2423 20 mg", "AZD2423 150 mg"],
            "matrix_unit": "Salmon transcript-level estimated counts; fractional values retained",
        },
        "metadata": summarize_metadata(metadata),
        "matrix": matrix_audit,
        "files": {
            "series_matrix": str(SERIES),
            "series_matrix_sha256": sha256(SERIES),
            "salmon_matrix": str(COUNTS),
            "human_gene_info": str(GENE_INFO),
            "human_gene_info_sha256": sha256(GENE_INFO),
        },
        "analysis_boundary": "No signature-outcome association was calculated by this audit script.",
    }
    with (TABLES / f"GSE302658_deposit_audit_{DATE}.json").open("w", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, ensure_ascii=False)

    print(json.dumps(qc["metadata"], indent=2, ensure_ascii=False))
    print(json.dumps(qc["matrix"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
