#!/usr/bin/env python3
"""Remove the ten bacterial false positives from the frozen MPDB release."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


FALSE_ACCESSIONS = {
    "MZ092733.1",
    "NG_047476.1",
    "NG_048062.1",
    "OK482771.1",
    "OK504483.1",
    "ON286761.1",
    "AB116388.1",
    "FJ917355.1",
    "KY623684.1",
    "NG_048058.1",
}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
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


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    records = read_fasta(args.source_fasta)
    if len(records) != 6133:
        raise SystemExit(f"Expected 6,133 source records, found {len(records)}")

    with args.manifest.open(encoding="utf-8", newline="") as handle:
        manifest = {row["seqid"]: row for row in csv.DictReader(handle, delimiter="\t")}

    kept: list[tuple[str, str]] = []
    removed: list[dict[str, str]] = []
    for header, sequence in records:
        seqid = header.split("|", 1)[0]
        row = manifest[seqid]
        if row["accession"] in FALSE_ACCESSIONS:
            removed.append(
                {
                    "seqid": seqid,
                    "accession": row["accession"],
                    "species": row["species"],
                    "length": str(len(sequence)),
                    "reason": "16S rRNA methylase/methyltransferase or plasmid locus, not 16S rRNA gene",
                }
            )
        else:
            kept.append((header, sequence))

    removed_accessions = Counter(row["accession"] for row in removed)
    expected = Counter(FALSE_ACCESSIONS)
    if removed_accessions != expected:
        raise SystemExit(f"False-positive removal mismatch: {removed_accessions} != {expected}")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(kept, output_dir / "maize_pathogens_all.fasta")

    by_category: dict[str, list[tuple[str, str]]] = {}
    for header, sequence in kept:
        category = header.split("|")[3]
        by_category.setdefault(category, []).append((header, sequence))
    for category, category_records in sorted(by_category.items()):
        write_fasta(category_records, output_dir / f"maize_pathogens_{category}.fasta")

    with (output_dir / "removed_false_positives.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = ["seqid", "accession", "species", "length", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(removed)

    counts = Counter(header.split("|")[3] for header, _ in kept)
    print(f"source records: {len(records)}")
    print(f"removed records: {len(removed)}")
    print(f"retained records: {len(kept)}")
    print("retained categories:", dict(sorted(counts.items())))


if __name__ == "__main__":
    main()
