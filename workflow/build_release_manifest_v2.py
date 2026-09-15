#!/usr/bin/env python3
"""Build SHA256SUMS.txt for tracked and release-candidate files."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "SHA256SUMS.txt"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


output = subprocess.check_output(
    ["git", "-C", str(REPO), "ls-files", "--cached", "--others", "--exclude-standard"],
    text=True,
    encoding="utf-8",
)
paths = sorted(
    REPO / line
    for line in output.splitlines()
    if line and line != "SHA256SUMS.txt" and (REPO / line).is_file()
)
MANIFEST.write_text(
    "".join(f"{sha256(path)}  {path.relative_to(REPO).as_posix()}\n" for path in paths),
    encoding="ascii",
    newline="\n",
)
print(f"manifest_records={len(paths)}")
