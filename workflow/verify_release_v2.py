#!/usr/bin/env python3
"""Validate the portable v2.0.0 release candidate without raw inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import py_compile
import re
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
V5 = REPO / "analysis_v5" / "public_clinical_extensions_2026-09-14"
REQUIRED = [
    REPO / "README.md",
    REPO / "CITATION.cff",
    REPO / ".zenodo.json",
    REPO / "VERSION_2_SOURCE_REGISTER.tsv",
    REPO / "VERSION_2_PROVENANCE.md",
    REPO / "analysis_v3/gene_definition_repair_v1/strict_background_addendum/consolidated_results/ASSOCIATION_IMPACT_ALL.tsv",
    REPO / "analysis_v3/gene_definition_repair_v1/strict_background_addendum/inputs/REPAIRED_MEMBERS.tsv",
    REPO / "analysis_v4/targeted_repair_2026-09-12/03_module_block_repair/results/M01_CORRECTED_BLOCK_INFLUENCE.tsv",
    REPO / "analysis_v4/closeout_20260912T062015Z/04_M03_spatial/SPATIAL_DONOR_EFFECTS.tsv",
    V5 / "06_integration/CLINICAL_EVIDENCE_MATRIX.tsv",
    V5 / "06_integration/MANUSCRIPT_CLINICAL_INTEGRATION_DRAFT.md",
    V5 / "tests/MODULE_SEMANTIC_TESTS.tsv",
    V5 / "MANIFEST_SHA256.tsv",
]
LOCAL_PATH_PATTERNS = [
    re.compile(r"[A-Za-z]:\\(?:Users\\[^\\]+|1 Codex project)"),
    re.compile("/mnt/" + "data/"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


candidate_output = subprocess.check_output(
    ["git", "-C", str(REPO), "ls-files", "--cached", "--others", "--exclude-standard"],
    text=True,
    encoding="utf-8",
)
candidate_paths = [REPO / line for line in candidate_output.splitlines() if line]

checks: list[tuple[str, bool, str]] = []
for path in REQUIRED:
    checks.append((f"required:{path.relative_to(REPO).as_posix()}", path.is_file(), "required file"))

oversize = [p.relative_to(REPO).as_posix() for p in candidate_paths if p.is_file() and p.stat().st_size >= 95_000_000]
checks.append(("no_file_at_or_above_95MB", not oversize, ";".join(oversize)))

local_hits = []
for path in candidate_paths:
    if not path.is_file() or path == Path(__file__).resolve():
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        continue
    for pattern in LOCAL_PATH_PATTERNS:
        if pattern.search(text):
            local_hits.append(f"{path.relative_to(REPO).as_posix()}::{pattern.pattern}")
checks.append(("no_local_absolute_paths", not local_hits, ";".join(local_hits[:20])))

compile_errors = []
for path in candidate_paths:
    if path.suffix == ".py" and path.is_file():
        try:
            py_compile.compile(str(path), doraise=True)
        except Exception as exc:  # pragma: no cover - diagnostic path
            compile_errors.append(f"{path.relative_to(REPO).as_posix()}::{exc}")
checks.append(("python_sources_compile", not compile_errors, ";".join(compile_errors)))

manifest_errors = []
with (V5 / "MANIFEST_SHA256.tsv").open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = V5 / row["relative_path"]
        if not path.is_file():
            manifest_errors.append(f"missing:{row['relative_path']}")
        elif path.stat().st_size != int(row["size_bytes"]):
            manifest_errors.append(f"size:{row['relative_path']}")
        elif sha256(path) != row["sha256"]:
            manifest_errors.append(f"hash:{row['relative_path']}")
checks.append(("v5_manifest", not manifest_errors, ";".join(manifest_errors[:20])))

release_manifest_errors = []
manifest_records = 0
for line in (REPO / "SHA256SUMS.txt").read_text(encoding="ascii").splitlines():
    if not line.strip():
        continue
    manifest_records += 1
    digest, rel = line.split("  ", 1)
    path = REPO / rel
    if not path.is_file():
        release_manifest_errors.append(f"missing:{rel}")
    elif sha256(path) != digest:
        release_manifest_errors.append(f"hash:{rel}")
checks.append(("release_manifest", not release_manifest_errors, f"records={manifest_records};" + ";".join(release_manifest_errors[:20])))

with (V5 / "tests/MODULE_SEMANTIC_TESTS.tsv").open(encoding="utf-8", newline="") as handle:
    semantic = list(csv.DictReader(handle, delimiter="\t"))
status_field = next((x for x in ("status", "result", "passed") if semantic and x in semantic[0]), None)
semantic_ok = bool(status_field) and all(str(row[status_field]).upper() in {"PASS", "TRUE"} for row in semantic)
checks.append(("clinical_semantic_tests", semantic_ok, f"records={len(semantic)};status_field={status_field}"))

with (REPO / "VERSION_2_SOURCE_REGISTER.tsv").open(encoding="utf-8", newline="") as handle:
    sources = list(csv.DictReader(handle, delimiter="\t"))
checks.append(("source_register", len(sources) == 13, f"records={len(sources)}"))

try:
    metadata = json.loads((REPO / ".zenodo.json").read_text(encoding="utf-8"))
    metadata_ok = metadata.get("version") == "2.0.0" and len(metadata.get("creators", [])) == 2
except Exception as exc:
    metadata_ok = False
    metadata = {"error": str(exc)}
checks.append(("zenodo_metadata", metadata_ok, str(metadata.get("version", metadata))))

failed = [name for name, passed, _ in checks if not passed]
print("check\tstatus\tdetail")
for name, passed, detail in checks:
    print(f"{name}\t{'PASS' if passed else 'FAIL'}\t{detail}")
print(f"SUMMARY\t{'PASS' if not failed else 'FAIL'}\tchecks={len(checks)};failed={len(failed)}")
sys.exit(1 if failed else 0)
