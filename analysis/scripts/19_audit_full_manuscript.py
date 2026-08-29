#!/usr/bin/env python3
"""Audit full-manuscript structure, citations, boundaries, and critical numbers."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
MANUSCRIPT = PHASE / "manuscript"
DATE = "2026-08-27"
DRAFT = MANUSCRIPT / f"FULL_MANUSCRIPT_DRAFT_V1_{DATE}.md"
REGISTRY = MANUSCRIPT / f"MANUSCRIPT_NUMERICAL_CLAIM_REGISTRY_V1_{DATE}.tsv"
BIBLIOGRAPHY = MANUSCRIPT / f"BIBLIOGRAPHIC_AND_DATASET_PROVENANCE_AUDIT_V1_{DATE}.tsv"
OUTPUT = MANUSCRIPT / f"FULL_MANUSCRIPT_AUDIT_V1_{DATE}.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def section(text: str, start: str, end: str) -> str:
    return text.split(start, 1)[1].split(end, 1)[0]


def main() -> None:
    text = DRAFT.read_text(encoding="utf-8")
    registry = pd.read_csv(REGISTRY, sep="\t")
    bibliography = pd.read_csv(BIBLIOGRAPHY, sep="\t")

    abstract = section(text, "## Abstract\n", "**Keywords:**")
    abstract_words = re.findall(r"\b[\w'-]+\b", abstract)
    body_before_references = text.split("## References", 1)[0]
    reference_section = text.split("## References", 1)[1]

    required_tokens = {
        "N001": ["656", "3,173", "744", "509", "60-gene"],
        "N002": ["196-gene", "522-gene", "2,951-gene", "Twenty-six"],
        "N004": ["concordance was 0.516", "Q=0.000400", "g=0.863"],
        "N006": ["Q=0.00520", "g=1.334", "Q value was 0.0281"],
        "N009": ["g=1.674", "Q value was 0.00504"],
        "N016": ["rho=0.070", "P=0.244", "P=0.280", "P=0.995"],
        "N017": ["rho=0.198", "P=0.0273", "P=0.221"],
        "N018": ["only seven source-up and 33 source-down"],
    }
    numerical_checks = {
        claim_id: {
            "missing_tokens": [token for token in tokens if token not in text],
            "pass": all(token in text for token in tokens),
        }
        for claim_id, tokens in required_tokens.items()
    }

    citation_checks = {}
    for ref_id in range(1, 20):
        body_pattern = re.compile(rf"\[(?:\d+[,-]?\s*)*\b{ref_id}\b(?:\s*[,-]\s*\d+)*\]")
        reference_pattern = re.compile(rf"(?m)^{ref_id}\. ")
        citation_checks[str(ref_id)] = {
            "cited_in_body": bool(body_pattern.search(body_before_references)),
            "listed_in_references": bool(reference_pattern.search(reference_section)),
        }

    required_sections = [
        "## Abstract",
        "## Introduction",
        "## Results",
        "## Discussion",
        "## Conclusions",
        "## Methods",
        "## Draft figure legends",
        "## Data availability",
        "## Code availability",
        "## Ethics statement",
        "## References",
    ]
    structural_checks = {
        "required_sections_present": all(item in text for item in required_sections),
        "abstract_word_count_150_to_300": 150 <= len(abstract_words) <= 300,
        "no_reference_placeholders": "[REF-" not in text,
        "nineteen_bibliography_rows": len(bibliography) == 19,
        "nineteen_references_listed": all(item["listed_in_references"] for item in citation_checks.values()),
        "all_references_cited": all(item["cited_in_body"] for item in citation_checks.values()),
        "numerical_registry_unchanged_and_passed": bool(len(registry) == 18 and (registry["audit_status"] == "PASS").all()),
        "preprint_boundary": "non-peer-reviewed" in text and "not peer reviewed" in text,
        "causal_boundary": "does not establish causal propagation" in text,
        "bulk_cell_origin_boundary": "cannot be assigned to purified neurons" in text,
        "ocular_axis_boundary": "do not support a unified TG-cornea-sural axis" in text,
        "biomarker_boundary": "Tissue-program transportability and biomarker portability are therefore distinct hypotheses" in text,
        "tear_coverage_boundary": "could not directly test the neural programs" in text,
        "funding_not_asserted_without_confirmation": "no project number is asserted in this draft" in text,
    }

    all_pass = (
        all(item["pass"] for item in numerical_checks.values())
        and all(structural_checks.values())
        and all(item["cited_in_body"] and item["listed_in_references"] for item in citation_checks.values())
    )
    audit = {
        "date": DATE,
        "draft": DRAFT.name,
        "draft_sha256": sha256(DRAFT),
        "abstract_word_count": len(abstract_words),
        "registry_sha256": sha256(REGISTRY),
        "bibliography_audit_sha256": sha256(BIBLIOGRAPHY),
        "numerical_checks": numerical_checks,
        "citation_checks": citation_checks,
        "structural_and_boundary_checks": structural_checks,
        "all_checks_pass": all_pass,
        "audit_scope": "Full-draft structure, citation completeness, critical numerical tokens, and required evidence-boundary language. Journal formatting and final language editing remain separate tasks.",
    }
    OUTPUT.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
