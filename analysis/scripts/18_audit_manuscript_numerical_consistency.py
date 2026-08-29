#!/usr/bin/env python3
"""Audit critical manuscript numbers and evidence-boundary wording."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
MANUSCRIPT = PHASE / "manuscript"
DATE = "2026-08-27"
DRAFT = MANUSCRIPT / f"RESULTS_METHODS_DRAFT_V1_{DATE}.md"
REGISTRY = MANUSCRIPT / f"MANUSCRIPT_NUMERICAL_CLAIM_REGISTRY_V1_{DATE}.tsv"
OUTPUT = MANUSCRIPT / f"RESULTS_METHODS_NUMERICAL_AUDIT_V1_{DATE}.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    text = DRAFT.read_text(encoding="utf-8")
    registry = pd.read_csv(REGISTRY, sep="\t")
    assert len(registry) == 18 and (registry["audit_status"] == "PASS").all()

    required_tokens = {
        "N001": ["656", "3,173", "744", "509", "60-gene"],
        "N002": ["196-gene", "522-gene", "2,951-gene", "Twenty-six"],
        "N003": ["Q=0.0103", "Q=0.0492", "Q=0.00475"],
        "N004": ["concordance was 0.516", "Q=0.000400", "0.00809", "g=0.863", "Q value was 0.0933"],
        "N005": ["Q=0.000150", "g=1.164", "0.932", "g=-0.393"],
        "N006": ["137 mapped", "252 mapped", "Q=0.00520", "0.02142", "g=1.334", "Q value was 0.0281"],
        "N007": ["0.03164", "g=1.243", "Q=0.0476"],
        "N008": ["Q value of 0.922", "g=-0.063", "g=-1.209"],
        "N009": ["Eighty-five", "22 source-down", "Q value of 0.0236", "0.01767", "g=1.674", "Q value was 0.00504"],
        "N010": ["Q=0.0530", "g=0.573", "Q=0.155"],
        "N011": ["g=0.797", "g=1.277", "g=-1.180"],
        "N012": ["g=2.273", "Q=0.0115", "g=0.820", "Q=0.000300"],
        "N013": ["-0.531", "-0.248", "-0.069"],
        "N014": ["Q=0.0444", "g=0.743", "Q=0.0255", "g=0.967"],
        "N015": ["g=-0.861", "I-squared=82.8%"],
        "N016": ["100 participants", "rho=0.070", "P=0.244", "P=0.280", "P=0.995"],
        "N017": ["95 paired", "rho=0.198", "P=0.0273", "P=0.0554", "-0.00013", "P=0.221"],
        "N018": ["only seven source-up and 33 source-down"],
    }

    checks = {}
    for claim_id, tokens in required_tokens.items():
        missing = [token for token in tokens if token not in text]
        checks[claim_id] = {"required_tokens": tokens, "missing_tokens": missing, "pass": not missing}

    structural_checks = {
        "six_result_sections": text.count("### ") >= 6,
        "six_figure_legends": all(f"### Figure {index}." in text for index in range(1, 7)),
        "preprint_label_present": "non-peer-reviewed" in text,
        "post_primary_label_present": "post-primary" in text.lower(),
        "provisional_label_present": "provisionally supportive" in text,
        "causal_boundary_present": "do not establish longitudinal progression or causal propagation" in text,
        "biomarker_boundary_present": "does not imply a stable blood or tear biomarker" in text,
        "not_testable_tear_present": "classified as not testable" in text,
        "results_present": "## Results" in text,
        "methods_present": "## Methods" in text,
    }
    all_pass = all(item["pass"] for item in checks.values()) and all(structural_checks.values())
    audit = {
        "date": DATE,
        "draft": DRAFT.name,
        "draft_sha256": sha256(DRAFT),
        "registry": REGISTRY.name,
        "registry_sha256": sha256(REGISTRY),
        "registry_rows": int(len(registry)),
        "claim_token_checks": checks,
        "structural_and_boundary_checks": structural_checks,
        "all_checks_pass": all_pass,
        "audit_scope": "Critical number tokens and required evidence-boundary language; not a substitute for final copyediting or citation verification.",
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
    print(json.dumps(audit, ensure_ascii=False, indent=2))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
