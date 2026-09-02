#!/usr/bin/env python3
"""Rebuild SINTAX so every entry matches a sequence in the release FASTA.

Reads the release all-sequences FASTA and the existing SINTAX file, keeps one
SINTAX entry per FASTA sequence in FASTA order, and writes a consistent file.
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.environ.get(
    "FASTA",
    os.path.join(ROOT, "release", "sequences", "maize_pathogens_all.fasta"),
)
SINTAX = os.environ.get(
    "SINTAX",
    os.path.join(ROOT, "release", "maize_pathogens_taxonomy_sintax.fasta"),
)
OUT = os.environ.get(
    "OUT",
    os.path.join(ROOT, "release", "maize_pathogens_taxonomy_sintax.fasta"),
)


def read_fasta_ids(path):
    ids = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split("|", 1)[0].strip())
    return ids


def read_sintax(path):
    entries = {}
    current = None
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    entries[current] = ("".join(seq) + "\n")
                current = line[1:].split(";", 1)[0].strip()
                header = line.strip()
                seq = [header + "\n"]
            else:
                seq.append(line.rstrip("\n"))
        if current is not None:
            entries[current] = ("\n".join(seq) + "\n")
    return entries


def main():
    ids = read_fasta_ids(FASTA)
    entries = read_sintax(SINTAX)
    missing = [sid for sid in ids if sid not in entries]
    if missing:
        raise SystemExit(f"missing SINTAX entries: {missing[:10]}")
    extra = sorted(set(entries) - set(ids))
    with open(OUT, "w", encoding="utf-8") as fh:
        for sid in ids:
            fh.write(entries[sid])
    print(f"fasta_ids={len(ids)} sintax_entries={len(entries)} "
          f"removed_stale={len(extra)}")
    if extra:
        print("stale removed:", extra[:10])


if __name__ == "__main__":
    main()
