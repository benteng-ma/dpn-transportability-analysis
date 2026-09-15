from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
PROJECT = RUN.parents[1]
STAGING_PARENT = RUN / "review_bundle"
STAGING = STAGING_PARENT / "DPN_PUBLIC_CLINICAL_EXTENSIONS_REVIEW_BUNDLE"
ZIP_PATH = PROJECT / "analysis_v5" / "DPN_PUBLIC_CLINICAL_EXTENSIONS_REVIEW_BUNDLE.zip"
SHA_PATH = ZIP_PATH.with_suffix(ZIP_PATH.suffix + ".sha256")
VERIFY_PATH = RUN / "logs" / "FINAL_BUNDLE_VERIFICATION.json"

include = [
    "00_admin",
    "01_GSE302658",
    "02_GSE286347",
    "03_GSE14806x",
    "04_JCI184075",
    "05_GSE295206",
    "06_integration",
    "code",
    "tests",
    "logs",
    "README_REVIEW_BUNDLE.md",
    "KNOWN_ISSUES.md",
]

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def safe_remove_tree(path: Path) -> None:
    resolved = path.resolve()
    parent = STAGING_PARENT.resolve()
    if resolved.parent != parent or resolved.name != "DPN_PUBLIC_CLINICAL_EXTENSIONS_REVIEW_BUNDLE":
        raise RuntimeError(f"Refusing to remove unexpected path: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)

STAGING_PARENT.mkdir(parents=True, exist_ok=True)
safe_remove_tree(STAGING)
STAGING.mkdir(parents=True)

for name in include:
    src = RUN / name
    if not src.exists():
        raise FileNotFoundError(src)
    dst = STAGING / name
    if src.is_dir():
        shutil.copytree(
            src,
            dst,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache"),
        )
    else:
        shutil.copy2(src, dst)

receipt = {
    "bundle_name": STAGING.name,
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "project_identity": "phase0_6_human_dpn_stage_projection",
    "baseline": "v3.1 final strict background with valid replacement list",
    "scope": "public-data-only M01-M06 clinical extensions",
    "manuscript_modified": False,
    "github_or_zenodo_published": False,
    "controlled_access_requested": False,
    "raw_fastq_downloaded": False,
    "included_top_level_items": include,
}
(STAGING / "BUNDLE_CONTENT_RECEIPT.json").write_text(
    json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
)

manifest_path = STAGING / "MANIFEST_SHA256.tsv"
records = []
for path in sorted(p for p in STAGING.rglob("*") if p.is_file() and p != manifest_path):
    records.append((path.relative_to(STAGING).as_posix(), path.stat().st_size, sha256(path)))
with manifest_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["relative_path", "size_bytes", "sha256"])
    writer.writerows(records)

if ZIP_PATH.exists():
    ZIP_PATH.unlink()
with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
    for path in sorted(p for p in STAGING.rglob("*") if p.is_file()):
        archive.write(path, arcname=f"{STAGING.name}/{path.relative_to(STAGING).as_posix()}")

zip_hash = sha256(ZIP_PATH)
SHA_PATH.write_text(f"{zip_hash}  {ZIP_PATH.name}\n", encoding="ascii")

manifest_ok = True
for rel, size, digest in records:
    path = STAGING / rel
    if path.stat().st_size != size or sha256(path) != digest:
        manifest_ok = False
        break
with zipfile.ZipFile(ZIP_PATH, "r") as archive:
    bad_member = archive.testzip()
    archive_names = archive.namelist()

verification = {
    "zip_path": str(ZIP_PATH),
    "zip_size_bytes": ZIP_PATH.stat().st_size,
    "zip_sha256": zip_hash,
    "manifest_records": len(records),
    "manifest_verified": manifest_ok,
    "zip_crc_bad_member": bad_member,
    "zip_member_count": len(archive_names),
    "manifest_self_reference": any(r[0] == "MANIFEST_SHA256.tsv" for r in records),
    "zip_self_reference": any(Path(n).name == ZIP_PATH.name for n in archive_names),
    "status": "PASS" if manifest_ok and bad_member is None else "FAIL",
}
VERIFY_PATH.write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
print(json.dumps(verification, indent=2))
