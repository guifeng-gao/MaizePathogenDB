#!/usr/bin/env python3
"""Summarize marker length and GC content after extraction."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header.split("|", 1)[0]] = (header, "".join(sequence))
                header = line[1:].strip()
                sequence = []
            else:
                sequence.append(line.strip())
    if header is not None:
        records[header.split("|", 1)[0]] = (header, "".join(sequence))
    return records


def gc_content(sequence: str) -> float:
    canonical = [base for base in sequence.upper() if base in "ACGT"]
    if not canonical:
        return float("nan")
    return (canonical.count("G") + canonical.count("C")) / len(canonical) * 100


def quantile(values: list[int], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    fasta = read_fasta(args.sequences)
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        meta = {row["seqid"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(fasta) != set(meta):
        raise SystemExit("FASTA and marker manifest contain different sequence IDs")

    per_sequence: list[dict[str, str]] = []
    grouped: dict[str, list[dict[str, float]]] = defaultdict(list)
    for sid, (header, sequence) in fasta.items():
        row = meta[sid]
        length = len(sequence)
        gc = gc_content(sequence)
        item = {
            "seqid": sid,
            "taxid": row["taxid"],
            "species": row["species"],
            "category": row["category"],
            "accession": row["accession"],
            "source_length": row["source_length"],
            "marker_length": row["marker_length"],
            "length": length,
            "gc_percent": f"{gc:.4f}",
            "marker_region": row["marker_region"],
            "extraction_method": row["extraction_method"],
        }
        per_sequence.append(item)
        grouped[row["category"]].append({"length": length, "gc": gc})

    summary: list[dict[str, str]] = []
    for category, values in sorted(grouped.items()):
        lengths = [int(item["length"]) for item in values]
        gcs = [float(item["gc"]) for item in values]
        summary.append(
            {
                "category": category,
                "n": str(len(values)),
                "length_min": str(min(lengths)),
                "length_q1": f"{quantile(lengths, 0.25):.1f}",
                "length_median": f"{statistics.median(lengths):.1f}",
                "length_q3": f"{quantile(lengths, 0.75):.1f}",
                "length_max": str(max(lengths)),
                "length_mean": f"{statistics.mean(lengths):.1f}",
                "gc_min": f"{min(gcs):.1f}",
                "gc_q1": f"{quantile(gcs, 0.25):.1f}",
                "gc_median": f"{statistics.median(gcs):.1f}",
                "gc_q3": f"{quantile(gcs, 0.75):.1f}",
                "gc_max": f"{max(gcs):.1f}",
                "gc_mean": f"{statistics.mean(gcs):.1f}",
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "length_gc_per_sequence.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = list(per_sequence[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(per_sequence)
    with (args.output_dir / "length_gc_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = list(summary[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    for row in summary:
        print(
            row["category"],
            "n=", row["n"],
            "length=", row["length_min"], row["length_median"], row["length_max"],
            "GC=", row["gc_min"], row["gc_median"], row["gc_max"],
        )


if __name__ == "__main__":
    main()
