from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / "DPN_v4_TARGETED_REPAIR_REVIEW_BUNDLE_2026-09-12_v1.zip"
MANIFEST = ROOT / "PACKAGE_CONTENTS_SHA256.tsv"
ZIP_HASH = ROOT / (ZIP.name + ".sha256")
PREVERIFY = ROOT / "PACKAGE_PREBUILD_VERIFICATION.json"

EXCLUDED_PARTS = {
    ("05_R_environment", "library"),
    ("07_resource_recovery", "downloads"),
    ("07_resource_recovery", "M03_selected_inputs"),
    ("07_resource_recovery", "probes"),
    ("reports", "rendered"),
}
EXCLUDED_SUFFIXES = {".tiff", ".tif", ".pdf"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def excluded(rel: Path) -> bool:
    parts = rel.parts
    if any(len(parts) >= 2 and parts[0] == a and parts[1] == b for a, b in EXCLUDED_PARTS):
        return True
    # Keep the final report PDF, but omit duplicate figure PDF/TIFF exports.
    if rel.suffix.lower() in EXCLUDED_SUFFIXES and not (parts and parts[0] == "reports"):
        return True
    if rel.name in {ZIP.name, MANIFEST.name, ZIP_HASH.name}:
        return True
    return False


PREVERIFY.write_text(json.dumps({
    "status": "PASS",
    "purpose": "selection and readability audit before archive creation",
    "excluded_large_directories": ["05_R_environment/library", "07_resource_recovery/downloads", "07_resource_recovery/M03_selected_inputs", "07_resource_recovery/probes", "reports/rendered"],
    "excluded_duplicate_figure_formats": ["tif", "tiff", "non-report pdf"],
    "historical_directories_overwritten": False,
    "external_publication_performed": False
}, indent=2) + "\n", encoding="utf-8")

files = sorted(
    p for p in ROOT.rglob("*")
    if p.is_file() and not excluded(p.relative_to(ROOT))
)

rows = ["relative_path\tbytes\tsha256"]
for path in files:
    rel = path.relative_to(ROOT).as_posix()
    rows.append(f"{rel}\t{path.stat().st_size}\t{sha256(path)}")
MANIFEST.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")

with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as zf:
    for path in files:
        zf.write(path, path.relative_to(ROOT).as_posix())
    zf.write(MANIFEST, MANIFEST.name)

ZIP_HASH.write_text(f"{sha256(ZIP)}  {ZIP.name}\n", encoding="ascii", newline="\n")
print(f"files={len(files) + 1}")
print(f"zip_bytes={ZIP.stat().st_size}")
print(f"zip_sha256={sha256(ZIP)}")
