#!/usr/bin/env python3
"""Build a marker-consistent SINTAX FASTA using taxonomy labels from the old file."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker-fasta", type=Path, required=True)
    parser.add_argument("--old-sintax", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    marker = read_fasta(args.marker_fasta)
    old = {
        header.split(";", 1)[0].split("|", 1)[0]: header
        for header, _ in read_fasta(args.old_sintax)
    }
    missing = [
        header.split("|", 1)[0]
        for header, _ in marker
        if header.split("|", 1)[0] not in old
    ]
    if missing:
        raise SystemExit(f"Missing taxonomy labels: {missing[:10]}")

    with args.output.open("w", encoding="utf-8") as handle:
        for header, sequence in marker:
            taxonomy_header = old[header.split("|", 1)[0]]
            handle.write(f">{taxonomy_header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    print(f"wrote {args.output} with {len(marker)} records")


if __name__ == "__main__":
    main()
