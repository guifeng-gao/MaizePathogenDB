#!/usr/bin/env python3
"""Reclassify the 260-sample ASV set against the corrected MPDB database."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from collections import Counter
from pathlib import Path


def genus_of(species: str) -> str:
    tokens = species.replace("/", " ").replace("(", " ").replace(")", " ").split()
    return tokens[0].strip("|,.;") if tokens else ""


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["seqid"]: row for row in csv.DictReader(handle, delimiter="\t")}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blastn", required=True)
    parser.add_argument("--asv-fasta", type=Path, required=True)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--marker-manifest", type=Path, required=True)
    parser.add_argument("--old-joined", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    manifest = read_manifest(args.marker_manifest)
    fields, old_rows = read_tsv(args.old_joined)
    cmd = [
        args.blastn,
        "-query",
        str(args.asv_fasta),
        "-db",
        str(args.db),
        "-outfmt",
        "6 qseqid sseqid pident qcovs",
        "-max_target_seqs",
        "1",
        "-evalue",
        "1e-5",
        "-num_threads",
        "8",
    ]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=7200)
    hits: dict[str, tuple[str, float, float]] = {}
    for line in proc.stdout.splitlines():
        fields_hit = line.split("\t")
        if len(fields_hit) < 4:
            continue
        qid, sseqid, pident, qcovs = fields_hit
        if qid not in hits:
            hits[qid] = (sseqid, float(pident), float(qcovs))

    updated: list[dict[str, object]] = []
    for row in old_rows:
        qid = row["qseqid"]
        hit = hits.get(qid)
        if hit:
            sseqid, pident, qcovs = hit
            sid = sseqid.split("|", 1)[0]
            marker = manifest[sid]
            species = marker["species"]
            genus = genus_of(species)
            hit_flag = 1
        else:
            sseqid = ""
            pident = 0.0
            qcovs = 0.0
            species = ""
            genus = ""
            hit_flag = 0
        species_ok = bool(hit_flag and pident >= 99.0 and qcovs >= 90.0)
        genus_ok = bool(hit_flag and pident >= 95.0 and qcovs >= 70.0)
        unite_genus = row.get("unite_genus", "")
        row.update(
            {
                "mpdb_subject": species,
                "mpdb_pident": pident,
                "mpdb_qcovs": qcovs,
                "mpdb_species": species,
                "mpdb_genus": genus,
                "mpdb_species_hit_ok": species_ok,
                "mpdb_genus_hit_ok": genus_ok,
                "mpdb_pathogen_species": species_ok,
                "mpdb_pathogen_genus": genus_ok,
                "mpdb_unite_genus_agree": bool(genus_ok and genus and unite_genus and genus.lower() == unite_genus.lower()),
            }
        )
        updated.append(row)

    write_tsv(args.output, fields, updated)
    both = sum(1 for r in updated if str(r["mpdb_pathogen_genus"]) == "True" and str(r["unite_pathogen_genus"]) == "True")
    mpdb_only = sum(1 for r in updated if str(r["mpdb_pathogen_genus"]) == "True" and str(r["unite_pathogen_genus"]) == "False")
    unite_only = sum(1 for r in updated if str(r["mpdb_pathogen_genus"]) == "False" and str(r["unite_pathogen_genus"]) == "True")
    rescued = [
        r
        for r in updated
        if (not r.get("unite_genus") or "Incertae_sedis" in str(r.get("unite_genus")))
        and str(r.get("mpdb_pathogen_genus")) == "True"
    ]

    def totals(predicate):
        selected = [r for r in updated if predicate(r)]
        return len(selected), sum(int(r["reads"]) for r in selected)

    result = {
        "n_asvs": len(updated),
        "mpdb_species_asvs": totals(lambda r: str(r["mpdb_pathogen_species"]) == "True"),
        "mpdb_genus_asvs": totals(lambda r: str(r["mpdb_pathogen_genus"]) == "True"),
        "unite_species_asvs": totals(lambda r: str(r["unite_pathogen_species"]) == "True"),
        "unite_genus_asvs": totals(lambda r: str(r["unite_pathogen_genus"]) == "True"),
        "agreement": {"both": both, "mpdb_only": mpdb_only, "unite_only": unite_only},
        "rescued_asvs": len(rescued),
        "rescued_reads": sum(int(r["reads"]) for r in rescued),
        "top_rescued": [],
    }
    # Rebuild top rescued from ASV/read counts because Counter over tuples is not a genus aggregate.
    genus_asv = Counter(str(r.get("mpdb_genus") or "Unknown") for r in rescued)
    genus_reads = Counter()
    for r in rescued:
        genus_reads[str(r.get("mpdb_genus") or "Unknown")] += int(r["reads"])
    result["top_rescued"] = [
        {"genus": genus, "asv": genus_asv[genus], "reads": genus_reads[genus]}
        for genus, _ in genus_reads.most_common(10)
    ]
    with args.summary.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
