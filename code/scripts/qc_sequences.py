#!/usr/bin/env python3
"""Run first-pass sequence QC on the standardized release FASTA."""

import csv
import os
import re
import subprocess
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEQ_DIR = os.environ.get(
    "SEQ_DIR", os.path.join(ROOT, "release", "sequences")
)
MANIFEST = os.environ.get(
    "MANIFEST", os.path.join(ROOT, "data", "sequence_manifest.tsv")
)
OUT_TSV = os.environ.get(
    "OUT_TSV", os.path.join(ROOT, "data", "sequence_qc_report.tsv")
)
QC_BLAST = os.path.join(ROOT, "data", "qc_blast")

BLASTN = os.environ.get("BLASTN", "blastn")
MAKEBLASTDB = os.environ.get("MAKEBLASTDB", "makeblastdb")

HEADER_RE = re.compile(
    r"^MPDB\d+\|\d+\|[^|]+\|[^|]+\|[A-Z]{1,2}_?\d+(?:\.\d+)?$"
)


def load_records():
    with open(MANIFEST, encoding="utf-8") as fh:
        manifest = {r["seqid"]: r for r in csv.DictReader(fh, delimiter="\t")}
    records = {}
    for path in os.listdir(SEQ_DIR):
        if not path.endswith(".fasta"):
            continue
        full = os.path.join(SEQ_DIR, path)
        header = None
        seq = []
        with open(full, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if header:
                        records[header.split("|", 1)[0]] = {
                            "header": header,
                            "sequence": "".join(seq),
                            "file": path,
                            **manifest[header.split("|", 1)[0]],
                        }
                    header = line[1:].strip()
                    seq = []
                else:
                    seq.append(line.strip())
        if header:
            records[header.split("|", 1)[0]] = {
                "header": header,
                "sequence": "".join(seq),
                "file": path,
                **manifest[header.split("|", 1)[0]],
            }
    return records


def marker_mismatch(record):
    description = record.get("original_header", "")
    category = record["category"]
    if category == "bacteria":
        return not re.search(r"16S", description, re.I)
    if category == "viruses":
        return not re.search(r"genome|segment|polyprotein|complete", description, re.I)
    if category in ("fungi", "oomycetes"):
        return not re.search(r"internal transcribed spacer|ITS|5\.8S|18S|28S", description, re.I)
    return True


def build_db(name, fasta_path):
    out = os.path.join(QC_BLAST, name)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = [MAKEBLASTDB, "-in", fasta_path, "-dbtype", "nucl", "-out", out]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    return out


def cross_category_hits():
    groups = {
        "bacteria": "maize_pathogens_bacteria.fasta",
        "viruses": "maize_pathogens_viruses.fasta",
        "fungi": "maize_pathogens_fungi.fasta",
    }
    dbs = {}
    for name, fname in groups.items():
        dbs[name] = build_db(name, os.path.join(SEQ_DIR, fname))

    suspicious = []
    for query_name, query_file in groups.items():
        for db_name in groups:
            if query_name == db_name:
                continue
            cmd = [
                BLASTN, "-query", os.path.join(SEQ_DIR, query_file),
                "-db", dbs[db_name], "-outfmt", "6 qseqid sseqid pident qcovs",
                "-max_target_seqs", "1", "-evalue", "1e-5",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            for line in result.stdout.strip().splitlines():
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                qseqid, sseqid, pident, qcovs = parts[0], parts[1], float(parts[2]), float(parts[3])
                if pident >= 70 and qcovs >= 50:
                    suspicious.append({
                        "query": qseqid, "subject": sseqid,
                        "query_db": query_name, "subject_db": db_name,
                        "pident": pident, "qcovs": qcovs,
                    })
    return suspicious


def main():
    records = load_records()
    seen_seqs = defaultdict(list)
    rows = []
    for seqid, record in records.items():
        seq = record["sequence"].upper()
        seen_seqs[seq].append(seqid)
        ambiguous = re.findall(r"[^ACGTN]", seq)
        rows.append({
            "seqid": seqid,
            "taxid": record["taxid"],
            "species": record["species"],
            "category": record["category"],
            "accession": record["accession"],
            "length": len(seq),
            "n_count": seq.count("N"),
            "ambiguous_count": len(ambiguous),
            "ambiguous_bases": "".join(sorted(set(ambiguous))),
            "header_format_ok": bool(HEADER_RE.match(record["header"])),
            "marker_consistent": not marker_mismatch(record),
            "notes": "",
        })

    duplicate_groups = {seq: ids for seq, ids in seen_seqs.items() if len(ids) > 1}
    for row in rows:
        if row["seqid"] in [i for ids in duplicate_groups.values() for i in ids]:
            group = next(ids for seq, ids in duplicate_groups.items() if row["seqid"] in ids)
            taxids = {records[i]["taxid"] for i in group}
            if len(taxids) == 1:
                row["notes"] = (
                    f"DUPLICATE_SEQUENCE_WITHIN_TAXID (n={len(group)}); "
                    "identical sequences from independent GenBank records, retained"
                )
            else:
                row["notes"] = (
                    f"DUPLICATE_SEQUENCE_CROSS_TAXID with "
                    f"{', '.join(i for i in group if i != row['seqid'])}"
                )

    suspicious = cross_category_hits()
    suspicious_by_query = {s["query"]: s for s in suspicious}
    for row in rows:
        if row["seqid"] in suspicious_by_query:
            hit = suspicious_by_query[row["seqid"]]
            row["notes"] += f"; CROSS_CATEGORY_HIT {hit['subject']} ({hit['query_db']}->{hit['subject_db']}, {hit['pident']:.1f}%, {hit['qcovs']:.0f}% qcov)".lstrip("; ")

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    columns = list(rows[0].keys())
    with open(OUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in columns) + "\n")

    print(f"wrote {OUT_TSV} ({len(rows)} records)")
    print("bad header format:", sum(1 for r in rows if not r["header_format_ok"]))
    print("marker inconsistent:", sum(1 for r in rows if not r["marker_consistent"]))
    print("duplicate groups:", len(duplicate_groups))
    within = sum(1 for ids in duplicate_groups.values()
                 if len({records[i]["taxid"] for i in ids}) == 1)
    print("duplicate groups within taxid:", within)
    print("duplicate groups cross taxid:", len(duplicate_groups) - within)
    print("cross-category hits:", len(suspicious))
    lengths = [r["length"] for r in rows]
    print("length min/median/max:", min(lengths), sorted(lengths)[len(lengths)//2], max(lengths))


if __name__ == "__main__":
    main()
