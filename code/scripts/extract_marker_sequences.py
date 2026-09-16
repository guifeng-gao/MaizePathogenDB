#!/usr/bin/env python3
"""Build marker-only MPDB sequences from cleaned full-length records."""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


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


def seqid(header: str) -> str:
    return header.split("|", 1)[0]


def taxid(header: str) -> str:
    return header.split("|")[1]


def category(header: str) -> str:
    return header.split("|")[3]


def load_its_regions(directory: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = defaultdict(dict)
    for region, filename in (
        ("ITS1", "itsx.ITS1.fasta"),
        ("5.8S", "itsx.5_8S.fasta"),
        ("ITS2", "itsx.ITS2.fasta"),
    ):
        path = directory / filename
        for header, sequence in read_fasta(path):
            result[seqid(header)][region] = sequence
    return result


def load_its_positions(directory: Path) -> dict[str, dict[str, tuple[int, int]]]:
    result: dict[str, dict[str, tuple[int, int]]] = defaultdict(dict)
    path = directory / "itsx.positions.txt"
    interval_re = re.compile(r"(ITS1|5\.8S|ITS2):\s*(\d+)-(\d+)")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            qid = line.split("|", 1)[0].strip()
            for region, start, end in interval_re.findall(line):
                result[qid][region.replace("5.8S", "5.8S")] = (int(start), int(end))
    return result


def extract_longest_contiguous_its(
    sequence: str, intervals: dict[str, tuple[int, int]], minimum_length: int = 50
) -> tuple[str, str] | None:
    ordered = [
        (region, intervals[region])
        for region in ("ITS1", "5.8S", "ITS2")
        if region in intervals
    ]
    if not ordered:
        return None
    blocks: list[tuple[str, int, int]] = []
    current_regions: list[str] = []
    current_start = current_end = 0
    for region, (start, end) in ordered:
        if not current_regions:
            current_regions = [region]
            current_start, current_end = start, end
        elif start <= current_end + 1:
            current_regions.append(region)
            current_end = max(current_end, end)
        else:
            blocks.append(("-".join(current_regions), current_start, current_end))
            current_regions = [region]
            current_start, current_end = start, end
    if current_regions:
        blocks.append(("-".join(current_regions), current_start, current_end))
    region_label, start, end = max(blocks, key=lambda item: item[2] - item[1] + 1)
    extracted = sequence[start - 1 : end]
    if len(extracted) < minimum_length:
        return None
    return extracted, region_label


def select_its_sequence(
    regions: dict[str, str], minimum_region_length: int = 50
) -> tuple[str, str] | None:
    presence = [region for region in ("ITS1", "5.8S", "ITS2") if region in regions]
    if not presence:
        return None

    # A short HMM hit is not sufficient evidence for a clean marker segment.
    reliable = {
        region: sequence
        for region, sequence in regions.items()
        if len(sequence) >= minimum_region_length
    }
    if not reliable:
        return None

    if all(region in reliable for region in ("ITS1", "5.8S", "ITS2")):
        return (
            reliable["ITS1"] + reliable["5.8S"] + reliable["ITS2"],
            "ITS1-5.8S-ITS2",
        )
    if all(region in reliable for region in ("ITS1", "5.8S")):
        return reliable["ITS1"] + reliable["5.8S"], "ITS1-5.8S"
    if all(region in reliable for region in ("5.8S", "ITS2")):
        return reliable["5.8S"] + reliable["ITS2"], "5.8S-ITS2"
    if len(reliable) == 1:
        region, sequence = next(iter(reliable.items()))
        return sequence, region

    # When the internal 5.8S anchor is absent, do not join ITS1 and ITS2.
    if "ITS1" in reliable and "ITS2" in reliable:
        region = max(("ITS1", "ITS2"), key=lambda item: len(reliable[item]))
        return reliable[region], f"{region}-only"
    return None


def marker_only_its_description(description: str) -> bool:
    has_its = bool(
        re.search(r"internal transcribed spacer|(?:^|[^A-Za-z])ITS[12]?(?:[^A-Za-z]|$)|5\.8S", description)
    )
    has_flank = bool(
        re.search(r"18S|28S|small subunit|large subunit|SSU|LSU", description, re.I)
    )
    return has_its and not has_flank


def parse_barrnap_16s(path: Path) -> dict[str, str]:
    best: dict[str, tuple[int, str, str]] = {}
    for header, sequence in read_fasta(path):
        parts = header.split("::", 1)
        if len(parts) != 2 or not parts[0].startswith("16S_rRNA"):
            continue
        source_id = parts[1].split("|", 1)[0]
        current = best.get(source_id)
        if current is None or len(sequence) > current[0]:
            best[source_id] = (len(sequence), sequence, header)
    return {key: value[1] for key, value in best.items()}


def blast_remap_unresolved(
    unresolved: list[tuple[str, str]],
    references: list[tuple[str, str]],
    blastn: str,
    work_dir: Path,
) -> dict[str, tuple[str, str]]:
    if not unresolved or not references:
        return {}
    work_dir.mkdir(parents=True, exist_ok=True)
    query_path = work_dir / "unresolved.fasta"
    subject_path = work_dir / "marker_references.fasta"
    db_path = work_dir / "marker_references"
    write_fasta(unresolved, query_path)
    write_fasta(references, subject_path)
    subprocess.run(
        ["makeblastdb", "-in", str(subject_path), "-dbtype", "nucl", "-out", str(db_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run(
        [
            blastn,
            "-query",
            str(query_path),
            "-db",
            str(db_path),
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend bitscore",
            "-max_target_seqs",
            "20",
            "-evalue",
            "1e-10",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    query_lookup = {seqid(header): sequence for header, sequence in unresolved}
    query_taxid = {seqid(header): taxid(header) for header, _ in unresolved}
    best: dict[str, tuple[float, int, int, str]] = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 7:
            continue
        qid, subject, identity, length, qstart, qend, bitscore = fields
        qid = qid.split("|", 1)[0]
        if taxid(subject) != query_taxid[qid]:
            continue
        if float(identity) < 90.0 or int(length) < 50:
            continue
        candidate = (float(bitscore), int(qstart), int(qend), subject)
        if qid not in best or candidate[0] > best[qid][0]:
            best[qid] = candidate

    remapped: dict[str, tuple[str, str]] = {}
    for qid, (_, qstart, qend, subject) in best.items():
        sequence = query_lookup[qid]
        start = min(qstart, qend) - 1
        end = max(qstart, qend)
        extracted = sequence[start:end]
        if len(extracted) >= 50:
            remapped[qid] = (extracted, f"BLAST-remapped from {subject}")
    return remapped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cleaned-dir", type=Path, required=True)
    parser.add_argument("--itsx-root", type=Path, required=True)
    parser.add_argument("--barrnap-fasta", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--blastn", default="blastn")
    args = parser.parse_args()

    source_records = read_fasta(args.cleaned_dir / "maize_pathogens_all.fasta")
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        manifest = {row["seqid"]: row for row in csv.DictReader(handle, delimiter="\t")}

    selected: list[tuple[str, str]] = []
    manifest_rows: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    method_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter()

    fungal_region_data = {
        category_name: load_its_regions(args.itsx_root / category_name)
        for category_name in ("fungi", "oomycetes")
    }
    fungal_position_data = {
        category_name: load_its_positions(args.itsx_root / category_name)
        for category_name in ("fungi", "oomycetes")
    }
    barrnap_16s = parse_barrnap_16s(args.barrnap_fasta)

    selected_by_category: dict[str, list[tuple[str, str]]] = defaultdict(list)
    unresolved_its: dict[str, list[tuple[str, str]]] = defaultdict(list)
    pending: dict[str, tuple[str, str, str]] = {}

    for header, sequence in source_records:
        sid = seqid(header)
        cat = category(header)
        row = manifest[sid]
        method = ""
        marker_region = ""
        marker_sequence = ""

        if cat in ("fungi", "oomycetes"):
            selected_its = select_its_sequence(fungal_region_data[cat].get(sid, {}))
            if selected_its:
                marker_sequence, marker_region = selected_its
                method = "ITSx-region-extraction"
            else:
                pending[sid] = (header, sequence, cat)
                unresolved_its[cat].append((header, sequence))
                continue
        elif cat == "bacteria":
            if sid in barrnap_16s:
                marker_sequence = barrnap_16s[sid]
                marker_region = "16S rRNA"
                method = "barrnap-coordinate-extraction"
            else:
                description = row["original_header"].lower()
                if ("16s" in description or "16s ribosomal" in description) and "methyl" not in description:
                    marker_sequence = sequence
                    marker_region = "16S rRNA partial"
                    method = "marker-only-record-retained"
                else:
                    excluded.append(
                        {
                            "seqid": sid,
                            "accession": row["accession"],
                            "species": row["species"],
                            "category": cat,
                            "source_length": str(len(sequence)),
                            "reason": "No 16S rRNA gene could be located",
                        }
                    )
                    continue
        elif cat == "viruses":
            marker_sequence = sequence
            marker_region = "complete viral genome/segment"
            method = "viral-genome-retained"
        else:
            raise SystemExit(f"Unknown category: {cat}")

        selected.append((header, marker_sequence))
        selected_by_category[cat].append((header, marker_sequence))
        method_counts[method] += 1
        category_counts[cat] += 1
        manifest_rows.append(
            {
                "seqid": sid,
                "taxid": row["taxid"],
                "species": row["species"],
                "category": cat,
                "accession": row["accession"],
                "source_length": str(len(sequence)),
                "marker_length": str(len(marker_sequence)),
                "marker_region": marker_region,
                "extraction_method": method,
            }
        )

    # Recover marker fragments that ITSx could not anchor but that align to a
    # reliable same-TaxID ITS sequence already extracted above.
    for cat in ("fungi", "oomycetes"):
        references = [
            (header, sequence)
            for header, sequence in selected_by_category[cat]
            if header.split("|")[3] == cat
        ]
        remapped = blast_remap_unresolved(
            unresolved_its[cat],
            references,
            args.blastn,
            args.work_dir / f"blast_remap_{cat}",
        )
        for header, sequence in unresolved_its[cat]:
            sid = seqid(header)
            row = manifest[sid]
            coordinate_fallback = extract_longest_contiguous_its(
                sequence, fungal_position_data[cat].get(sid, {})
            )
            if sid not in remapped and coordinate_fallback is not None:
                marker_sequence, marker_region = coordinate_fallback
                method = "ITSx-coordinate-extraction"
            elif sid not in remapped and marker_only_its_description(row["original_header"]):
                marker_sequence = sequence
                marker_region = "marker-only partial ITS"
                method = "marker-only-record-retained"
            elif sid not in remapped:
                excluded.append(
                    {
                        "seqid": sid,
                        "accession": row["accession"],
                        "species": row["species"],
                        "category": cat,
                        "source_length": str(len(sequence)),
                        "reason": "No clean ITS region could be extracted",
                    }
                )
                continue
            else:
                marker_sequence, method_detail = remapped[sid]
                marker_region = method_detail
                method = "BLAST-remapped-ITS-fragment"
            selected.append((header, marker_sequence))
            selected_by_category[cat].append((header, marker_sequence))
            method_counts[method] += 1
            category_counts[cat] += 1
            manifest_rows.append(
                {
                    "seqid": sid,
                    "taxid": row["taxid"],
                    "species": row["species"],
                    "category": cat,
                    "accession": row["accession"],
                    "source_length": str(len(sequence)),
                    "marker_length": str(len(marker_sequence)),
                    "marker_region": marker_region,
                    "extraction_method": method,
                }
            )

    output_dir = args.output_dir
    write_fasta(selected, output_dir / "maize_pathogens_all.fasta")
    for cat, records in sorted(selected_by_category.items()):
        write_fasta(records, output_dir / f"maize_pathogens_{cat}.fasta")

    with (output_dir / "sequence_manifest_marker.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = [
            "seqid",
            "taxid",
            "species",
            "category",
            "accession",
            "source_length",
            "marker_length",
            "marker_region",
            "extraction_method",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    with (output_dir / "marker_extraction_exclusions.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = ["seqid", "accession", "species", "category", "source_length", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded)

    print("retained records:", len(selected))
    print("category counts:", dict(sorted(category_counts.items())))
    print("extraction methods:", dict(sorted(method_counts.items())))
    print("excluded records:", len(excluded))


if __name__ == "__main__":
    main()
