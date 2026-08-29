#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for every release file except itself."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "SHA256SUMS.txt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = sorted(
        (path for path in ROOT.rglob("*") if path.is_file() and path.resolve() != OUTPUT.resolve()),
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )
    lines = [f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in files]
    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"manifested_files={len(files)}")
    print(f"manifest={OUTPUT}")


if __name__ == "__main__":
    main()
