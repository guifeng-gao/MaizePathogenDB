#!/usr/bin/env python3
"""Assemble a public MPDB release without figures or figure-generation code."""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
from pathlib import Path


COPY_DIRS = (
    "blast_db",
    "code",
    "curation",
    "data",
    "sequences",
    "taxonomy",
    "validation",
    "web",
)
COPY_FILES = ("LICENSE", "README.md", "RESULTS.md")
IGNORE = shutil.ignore_patterns(
    ".DS_Store",
    "__pycache__",
    "*.pyc",
    "*.png",
    "*.jpg",
    "*.jpeg",
    "*.pdf",
    "*.R",
)
LOCAL_PATH = re.compile(r"/Users/[A-Za-z0-9_.-]+|/Volumes/[A-Za-z0-9_.-]+")
TEXT_EXTENSIONS = {
    ".py",
    ".sh",
    ".md",
    ".yaml",
    ".yml",
    ".txt",
    ".cfg",
    ".ini",
    ".json",
    ".html",
    ".js",
    ".tsv",
}


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def scan_tree(root: Path) -> None:
    for path in root.rglob("*"):
        if not path.is_file() or path.name == "CHECKSUMS.sha256":
            continue
        if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".pdf"}:
            raise SystemExit(f"generated image found: {path}")
        if path.suffix.lower() == ".r":
            raise SystemExit(f"figure-generation code found: {path}")
        if path.suffix.lower() in TEXT_EXTENSIONS:
            text = path.read_text(encoding="utf-8", errors="replace")
            if LOCAL_PATH.search(text):
                raise SystemExit(f"local path found in {path}")


def write_checksums(root: Path) -> None:
    lines = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.name == "CHECKSUMS.sha256":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(root)}")
    (root / "CHECKSUMS.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()

    source = args.source.resolve()
    destination = args.destination.resolve()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)

    for directory in COPY_DIRS:
        shutil.copytree(source / directory, destination / directory, ignore=IGNORE)
    for filename in COPY_FILES:
        copy_file(source / filename, destination / filename)

    scan_tree(destination)
    write_checksums(destination)
    print(f"public release written to {destination}")


if __name__ == "__main__":
    main()
