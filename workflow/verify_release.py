#!/usr/bin/env python3
"""Validate the public release structure, portability, and recorded checksums."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".cff", ".csv", ".json", ".md", ".py", ".txt", ".tsv", ".yaml", ".yml"}
WINDOWS_DRIVE = re.compile(r"(?<![A-Za-z0-9_])[A-Za-z]:[\\/]")
SENSITIVE_PATTERNS = {
    "email": re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", re.I),
    "private_key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "generic_secret": re.compile(r"(?i)(?:api[_-]?key|access[_-]?token|client[_-]?secret)\s*[:=]\s*['\"][^'\"]{8,}"),
}
EXPECTED_ROOT = {
    ".gitattributes",
    ".gitignore",
    ".zenodo.json",
    "README.md",
    "CITATION.cff",
    "LICENSE",
    "LICENSE-CONTENT.md",
    "RELEASE_NOTES.md",
    "requirements.txt",
    "environment.yml",
    "SOURCE_DATA_MANIFEST.tsv",
    "ZENODO_METADATA_DRAFT.md",
    "SHA256SUMS.txt",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-raw", action="store_true", help="Require and verify every non-redistributed source input.")
    parser.add_argument("--require-license", action="store_true", help="Require an author-approved LICENSE file.")
    return parser.parse_args()


def add(checks: list[dict], name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "pass": bool(passed), "detail": detail})


def main() -> None:
    args = parse_args()
    checks: list[dict] = []

    missing_root = sorted(name for name in EXPECTED_ROOT if not (ROOT / name).is_file())
    add(checks, "required_release_files", not missing_root, "missing=" + repr(missing_root))

    scripts = sorted((ROOT / "analysis" / "scripts").glob("[0-9][0-9]_*.py"))
    add(checks, "analysis_script_count", len(scripts) == 21, f"observed={len(scripts)} expected=21")
    syntax_errors = []
    for path in scripts + sorted((ROOT / "workflow").glob("*.py")):
        try:
            compile(path.read_text(encoding="utf-8"), str(path), "exec")
        except SyntaxError as error:
            syntax_errors.append(f"{path.relative_to(ROOT)}:{error.lineno}:{error.msg}")
    add(checks, "python_syntax", not syntax_errors, repr(syntax_errors))

    path_hits = []
    sensitive_hits = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for line_number, line in enumerate(text.splitlines(), start=1):
            if WINDOWS_DRIVE.search(line) or "/Users/" in line or "\\Users\\" in line:
                path_hits.append(f"{path.relative_to(ROOT)}:{line_number}")
            for label, pattern in SENSITIVE_PATTERNS.items():
                if pattern.search(line):
                    sensitive_hits.append(f"{label}:{path.relative_to(ROOT)}:{line_number}")
    add(checks, "no_local_absolute_paths", not path_hits, repr(path_hits[:20]))
    add(checks, "no_email_or_secret_patterns", not sensitive_hits, repr(sensitive_hits[:20]))

    manifest_path = ROOT / "SOURCE_DATA_MANIFEST.tsv"
    rows = []
    if manifest_path.is_file():
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    add(checks, "source_manifest_nonempty", len(rows) >= 20, f"rows={len(rows)}")
    pending = []
    missing_inputs = []
    mismatched_inputs = []
    for row in rows:
        relative = row.get("required_relative_path", "")
        path = ROOT / relative
        expected_hash = row.get("local_source_sha256", "")
        redistributed = row.get("redistributed", "").lower() == "yes"
        if expected_hash == "PENDING_GENERATION" or not expected_hash:
            pending.append(relative)
        should_check = redistributed or args.check_raw
        if should_check and not path.is_file():
            missing_inputs.append(relative)
        elif should_check and expected_hash not in {"", "PENDING_GENERATION"}:
            observed = sha256(path)
            if observed != expected_hash:
                mismatched_inputs.append(relative)
    add(checks, "source_manifest_checksums_finalized", not pending, f"pending={len(pending)}")
    add(checks, "required_selected_inputs_present", not missing_inputs, repr(missing_inputs[:20]))
    add(checks, "required_selected_inputs_match", not mismatched_inputs, repr(mismatched_inputs[:20]))

    provenance_path = ROOT / "data" / "processed" / "ocular" / "PROVENANCE.json"
    ocular_errors = []
    if provenance_path.is_file():
        for record in json.loads(provenance_path.read_text(encoding="utf-8")):
            path = ROOT / record["file"]
            if not path.is_file() or sha256(path) != record["sha256"]:
                ocular_errors.append(record["file"])
    else:
        ocular_errors.append(str(provenance_path.relative_to(ROOT)))
    add(checks, "ocular_derivative_provenance", not ocular_errors, repr(ocular_errors))

    figure_errors = []
    figure_dir = ROOT / "results" / "figures" / "main_composites"
    for number in range(1, 7):
        png_matches = list(figure_dir.glob(f"Figure{number}_*_2026-08-27.png"))
        pdf_matches = list(figure_dir.glob(f"Figure{number}_*_2026-08-27.pdf"))
        if len(png_matches) != 1 or len(pdf_matches) != 1:
            figure_errors.append(f"Figure{number}:file_count")
            continue
        with Image.open(png_matches[0]) as image:
            if image.width < 2500 or image.height < 2000:
                figure_errors.append(f"Figure{number}:dimensions={image.size}")
        if pdf_matches[0].stat().st_size < 20_000:
            figure_errors.append(f"Figure{number}:pdf_size")
    add(checks, "six_main_figure_pairs", not figure_errors, repr(figure_errors))

    results = list((ROOT / "results" / "tables").glob("*"))
    metadata = list((ROOT / "metadata").glob("*"))
    add(checks, "frozen_result_inventory", len(results) >= 70 and len(metadata) >= 8, f"tables={len(results)} metadata={len(metadata)}")

    doi_text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in [ROOT / "README.md", ROOT / "CITATION.cff", ROOT / "ZENODO_METADATA_DRAFT.md"]
        if path.is_file()
    )
    fake_doi = re.search(r"10\.5281/zenodo\.(?:XXXX|0000|TBD|PENDING)", doi_text, re.I)
    add(checks, "no_fabricated_repository_doi", fake_doi is None, "no fake Zenodo DOI pattern")

    release_metadata_text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace")
        for path in [ROOT / "README.md", ROOT / "CITATION.cff", ROOT / ".zenodo.json"]
        if path.is_file()
    )
    add(checks, "final_version_metadata", "1.0.0-rc1" not in release_metadata_text, "no rc1 marker in release metadata")

    approved_licenses = [path for path in ROOT.glob("LICENSE*") if path.name != "LICENSE_DECISION_REQUIRED.md"]
    license_pass = bool(approved_licenses) if args.require_license else True
    license_detail = "approved=" + repr([path.name for path in approved_licenses])
    if not args.require_license:
        license_detail += "; not required for local RC validation"
    add(checks, "author_approved_license", license_pass, license_detail)

    checksum_path = ROOT / "SHA256SUMS.txt"
    checksum_errors = []
    checksum_entries = {}
    if checksum_path.is_file():
        for line_number, line in enumerate(checksum_path.read_text(encoding="utf-8").splitlines(), start=1):
            if "  " not in line:
                checksum_errors.append(f"line_{line_number}:format")
                continue
            expected, relative = line.split("  ", 1)
            checksum_entries[relative] = expected
            path = ROOT / relative
            if not path.is_file():
                checksum_errors.append(f"{relative}:missing")
            elif sha256(path) != expected:
                checksum_errors.append(f"{relative}:hash")
        actual_files = {
            path.relative_to(ROOT).as_posix()
            for path in ROOT.rglob("*")
            if path.is_file() and path.resolve() != checksum_path.resolve()
        }
        if set(checksum_entries) != actual_files:
            checksum_errors.append(
                f"inventory_difference:manifest={len(checksum_entries)} actual={len(actual_files)}"
            )
    else:
        checksum_errors.append("SHA256SUMS.txt:missing")
    add(checks, "release_checksum_manifest", not checksum_errors, repr(checksum_errors[:20]))

    all_pass = all(item["pass"] for item in checks)
    report = {
        "release": ROOT.name,
        "check_raw": args.check_raw,
        "require_license": args.require_license,
        "checks": checks,
        "all_checks_pass": all_pass,
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
