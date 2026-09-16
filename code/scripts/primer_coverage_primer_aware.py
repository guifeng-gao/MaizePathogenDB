#!/usr/bin/env python3
"""Compute primer coverage on cleaned records containing primer-flanking anchors."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("validation_module", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_fasta(path: Path):
    records = []
    header = None
    sequence = []
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


def calculate(records, validation_module):
    sequences = defaultdict(list)
    for header, sequence in records:
        sequences[header.split("|")[3]].append(sequence)
    results = {}
    for pair_name, (forward, reverse) in validation_module.PRIMER_PAIRS.items():
        if pair_name.startswith("16S"):
            target_sequences = sequences["bacteria"]
            category = "bacteria"
        else:
            target_sequences = sequences["fungi"] + sequences["oomycetes"]
            category = "fungi+oomycetes"
        reverse_primers = {reverse} if reverse.endswith("R") or reverse in ("ITS2", "ITS4", "ITS4ngs") else set()
        covered = 0
        for sequence in target_sequences:
            forward_ok = validation_module.primer_site(validation_module.PRIMERS[forward], sequence)
            reverse_sequence = validation_module.revcomp(sequence) if reverse in reverse_primers else sequence
            reverse_ok = validation_module.primer_site(validation_module.PRIMERS[reverse], reverse_sequence)
            covered += int(forward_ok and reverse_ok)
        results[pair_name] = {
            "category": category,
            "n": len(target_sequences),
            "covered": covered,
            "coverage_pct": round(covered / len(target_sequences) * 100, 1) if target_sequences else None,
        }
    return {"validation": "primer_coverage", "definition": "primer-aware cleaned records", "results": results}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation-module", type=Path, required=True)
    parser.add_argument("--cleaned-records", type=Path, required=True)
    parser.add_argument("--marker-records", type=Path, required=True)
    parser.add_argument("--marker-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    validation_module = load_module(args.validation_module)
    with args.marker_manifest.open(encoding="utf-8", newline="") as handle:
        retained_ids = {row["seqid"] for row in csv.DictReader(handle, delimiter="\t")}

    cleaned = [
        (header, sequence)
        for header, sequence in read_fasta(args.cleaned_records)
        if header.split("|", 1)[0] in retained_ids
    ]
    marker = read_fasta(args.marker_records)
    aware_result = calculate(cleaned, validation_module)
    marker_result = calculate(marker, validation_module)
    marker_result["definition"] = "marker-only records; diagnostic only for ITS primers"

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, data in (
        ("primer_coverage.json", aware_result),
        ("primer_coverage_marker_only.json", marker_result),
    ):
        with (args.output_dir / filename).open("w", encoding="utf-8") as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2)

    for label, data in (("primer-aware", aware_result), ("marker-only", marker_result)):
        print(label)
        for pair, row in data["results"].items():
            print(pair, row["covered"], row["n"], row["coverage_pct"])


if __name__ == "__main__":
    main()
