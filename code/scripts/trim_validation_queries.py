#!/usr/bin/env python3
"""Trim a validation query set to the same marker space as the corrected MPDB."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


def load_extractor(path: Path):
    spec = importlib.util.spec_from_file_location("mpdb_extractor", path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extractor", type=Path, required=True)
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--query-meta", type=Path, required=True)
    parser.add_argument("--mpdb-marker", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-meta", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--itsx", default="ITSx")
    parser.add_argument("--barrnap", default="barrnap")
    parser.add_argument("--blastn", default="blastn")
    args = parser.parse_args()

    ex = load_extractor(args.extractor)
    records = ex.read_fasta(args.query_fasta)
    with args.query_meta.open(encoding="utf-8", newline="") as handle:
        meta = {row["qid"]: row for row in csv.DictReader(handle, delimiter="\t")}

    by_category: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for header, sequence in records:
        by_category[ex.category(header)].append((header, sequence))

    args.work_dir.mkdir(parents=True, exist_ok=True)
    category_dir = args.work_dir / "category_fasta"
    category_dir.mkdir(parents=True, exist_ok=True)
    for category_name, category_records in by_category.items():
        ex.write_fasta(category_records, category_dir / f"{category_name}.fasta")

    itsx_data = {}
    for category_name in ("fungi", "oomycetes"):
        if category_name not in by_category:
            continue
        directory = args.work_dir / f"itsx_{category_name}"
        directory.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                args.itsx,
                "-i",
                str(category_dir / f"{category_name}.fasta"),
                "-o",
                str(directory / "itsx"),
                "-t",
                "all",
                "--save_regions",
                "all",
                "--preserve",
                "T",
                "--cpu",
                "8",
                "--silent",
                "T",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        itsx_data[category_name] = {
            "regions": ex.load_its_regions(directory),
            "positions": ex.load_its_positions(directory),
        }

    barrnap_16s = {}
    if "bacteria" in by_category:
        barrnap_dir = args.work_dir / "barrnap"
        barrnap_dir.mkdir(parents=True, exist_ok=True)
        barrnap_fasta = barrnap_dir / "bacteria.16S.fasta"
        with barrnap_fasta.open("w", encoding="utf-8") as handle:
            subprocess.run(
                [
                    args.barrnap,
                    "--kingdom",
                    "bac",
                    "--threads",
                    "8",
                    "--quiet",
                    "--outseq",
                    str(barrnap_fasta),
                    str(category_dir / "bacteria.fasta"),
                ],
                check=True,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
            )
        barrnap_16s = ex.parse_barrnap_16s(barrnap_fasta)

    mpdb_records = ex.read_fasta(args.mpdb_marker)
    references_by_category: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for header, sequence in mpdb_records:
        references_by_category[ex.category(header)].append((header, sequence))

    selected: list[tuple[str, str]] = []
    extracted_meta: list[dict[str, str]] = []
    dropped: list[dict[str, str]] = []
    unresolved_its: dict[str, list[tuple[str, str]]] = defaultdict(list)
    method_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()

    for header, sequence in records:
        qid = ex.seqid(header)
        category_name = ex.category(header)
        row = meta.get(qid, {})
        title = row.get("title", header)
        marker_sequence = ""
        method = ""
        region = ""

        if category_name in ("fungi", "oomycetes"):
            selected_its = ex.select_its_sequence(itsx_data[category_name]["regions"].get(qid, {}))
            if selected_its:
                marker_sequence, region = selected_its
                method = "ITSx-region-extraction"
            else:
                unresolved_its[category_name].append((header, sequence))
                continue
        elif category_name == "bacteria":
            if qid in barrnap_16s:
                marker_sequence = barrnap_16s[qid]
                region = "16S rRNA"
                method = "barrnap-coordinate-extraction"
            elif re.search(r"16S", title, re.I) and not re.search(r"methylase|methyltransferase|plasmid", title, re.I):
                marker_sequence = sequence
                region = "16S rRNA partial"
                method = "marker-only-record-retained"
            else:
                dropped.append({"qid": qid, "category": category_name, "reason": "No clean 16S rRNA marker"})
                continue
        elif category_name == "viruses":
            marker_sequence = sequence
            region = "complete viral genome/segment"
            method = "viral-genome-retained"
        else:
            raise SystemExit(f"Unknown category: {category_name}")

        selected.append((header, marker_sequence))
        extracted_meta.append(
            {
                "qid": qid,
                "category": category_name,
                "source_length": str(len(sequence)),
                "marker_length": str(len(marker_sequence)),
                "marker_region": region,
                "extraction_method": method,
            }
        )
        method_counts[method] += 1
        category_counts[category_name] += 1

    for category_name in ("fungi", "oomycetes"):
        remapped = ex.blast_remap_unresolved(
            unresolved_its[category_name],
            references_by_category[category_name],
            args.blastn,
            args.work_dir / f"blast_remap_{category_name}",
        )
        for header, sequence in unresolved_its[category_name]:
            qid = ex.seqid(header)
            title = meta.get(qid, {}).get("title", header)
            coordinate_fallback = ex.extract_longest_contiguous_its(
                sequence, itsx_data[category_name]["positions"].get(qid, {})
            )
            if qid in remapped:
                marker_sequence, region = remapped[qid]
                method = "BLAST-remapped-ITS-fragment"
            elif coordinate_fallback is not None:
                marker_sequence, region = coordinate_fallback
                method = "ITSx-coordinate-extraction"
            elif ex.marker_only_its_description(title):
                marker_sequence = sequence
                region = "marker-only partial ITS"
                method = "marker-only-record-retained"
            else:
                dropped.append({"qid": qid, "category": category_name, "reason": "No clean ITS marker"})
                continue
            selected.append((header, marker_sequence))
            extracted_meta.append(
                {
                    "qid": qid,
                    "category": category_name,
                    "source_length": str(len(sequence)),
                    "marker_length": str(len(marker_sequence)),
                    "marker_region": region,
                    "extraction_method": method,
                }
            )
            method_counts[method] += 1
            category_counts[category_name] += 1

    ex.write_fasta(selected, args.output_fasta)
    with args.output_meta.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "qid",
            "category",
            "source_length",
            "marker_length",
            "marker_region",
            "extraction_method",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(extracted_meta)
    dropped_path = args.output_meta.with_name(args.output_meta.stem + "_dropped.tsv")
    with dropped_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["qid", "category", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(dropped)

    print(args.query_fasta.name, "input", len(records), "kept", len(selected), "dropped", len(dropped))
    print("category counts:", dict(sorted(category_counts.items())))
    print("methods:", dict(sorted(method_counts.items())))


if __name__ == "__main__":
    main()
