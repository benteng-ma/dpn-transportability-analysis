#!/usr/bin/env python3
"""Build a checksum manifest for phase 0.6 specifications, code and deliverables."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PHASE = Path(__file__).resolve().parents[2]
TABLES = PHASE / "results" / "tables"
DATE = "2026-08-27"
OUTPUT_TSV = TABLES / f"PHASE0_6_RUN_MANIFEST_{DATE}.tsv"
OUTPUT_JSON = TABLES / f"PHASE0_6_RUN_MANIFEST_SUMMARY_{DATE}.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def category(path: Path) -> str:
    relative = path.relative_to(PHASE)
    first = relative.parts[0]
    if len(relative.parts) == 1:
        return "frozen_spec_or_status"
    if first == "analysis":
        return "analysis_script"
    if first == "metadata":
        return "metadata"
    if first == "reports":
        return "report"
    if first == "manuscript":
        return "manuscript"
    if first == "results" and len(relative.parts) > 1:
        return f"result_{relative.parts[1]}"
    return first


def main() -> None:
    candidates = []
    candidates.extend(PHASE.glob("*.md"))
    for folder in [
        PHASE / "analysis" / "scripts",
        PHASE / "metadata",
        PHASE / "reports",
        PHASE / "results",
        PHASE / "manuscript",
    ]:
        candidates.extend(path for path in folder.rglob("*") if path.is_file())
    excluded = {OUTPUT_TSV.resolve(), OUTPUT_JSON.resolve()}
    candidates = sorted(
        {path.resolve() for path in candidates if path.resolve() not in excluded and "__pycache__" not in path.parts},
        key=lambda path: str(path).lower(),
    )
    rows = []
    for path in candidates:
        rows.append(
            {
                "relative_path": path.relative_to(PHASE).as_posix(),
                "category": category(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUTPUT_TSV, sep="\t", index=False)
    summary = {
        "analysis_date": DATE,
        "manifest_entry_count": int(len(manifest)),
        "total_manifested_bytes": int(manifest["size_bytes"].sum()),
        "category_counts": manifest["category"].value_counts().to_dict(),
        "manifest_sha256": sha256(OUTPUT_TSV),
        "scope_note": (
            "Specifications, scripts, metadata, reports, manuscript files and result deliverables are manifested. "
            "Raw-input checksums are recorded in dataset-specific QC JSON files rather than duplicated here."
        ),
    }
    with OUTPUT_JSON.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
