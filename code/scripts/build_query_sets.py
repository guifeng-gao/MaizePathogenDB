#!/usr/bin/env python3
"""Build fixed query sets Q2/Q3/Q4 for the MaizePathogenDB release."""

import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECIES_LIST = os.path.join(ROOT, "data", "species_list.tsv")
REFERENCE_FASTA = os.path.join(
    ROOT, "release", "sequences", "maize_pathogens_all.fasta"
)
OUT_DIR = os.path.join(ROOT, "docs", "validation", "query_sets")

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"
TOOL = "maize_pathogen_db_query_builder"
DELAY = 0.6
Q2_MAX_PER_SPECIES = int(os.environ.get("Q2_MAX_PER_SPECIES", "10"))
Q4_MAX_PER_SPECIES = int(os.environ.get("Q4_MAX_PER_SPECIES", "5"))
Q4_FETCH_MAX = int(os.environ.get("Q4_FETCH_MAX", "100"))

MARKER = {
    "bacteria": '("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title])',
    "viruses": '("complete genome"[Title] OR "polyprotein"[Gene])',
    "fungi": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])',
    "oomycetes": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])',
}


def eutils_get(endpoint, params):
    url = f"{ENTREZ_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": f"MaizePathogenDB/{TOOL} ({EMAIL})"}
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.read().decode("utf-8", errors="replace")
        except Exception as exc:
            last = exc
            time.sleep(4 + attempt * 4)
    raise RuntimeError(f"eutils failed: {url}: {last}")


def normalize(seq):
    return re.sub(r"[^ACGTNacgtn]", "", seq).upper()


def load_reference_seqs():
    refs = set()
    header = None
    seq = []
    with open(REFERENCE_FASTA, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    refs.add(normalize("".join(seq)))
                header = line
                seq = []
            else:
                seq.append(line.strip())
        if header:
            refs.add(normalize("".join(seq)))
    return refs


def marker_ok(description, category):
    if category == "bacteria":
        return bool(re.search(r"16S", description, re.I))
    if category == "viruses":
        return bool(re.search(r"genome|segment|polyprotein|complete", description, re.I))
    return bool(re.search(r"internal transcribed spacer|ITS|5\.8S|18S", description, re.I))


def fetch_records(
    taxid,
    category,
    max_records=50,
    date_expr='"2020/01/01"[PDAT] : "2026/08/25"[PDAT]',
):
    query = (
        f'txid{taxid}[Organism] AND {MARKER[category]} '
        f"AND ({date_expr})"
    )
    params = {
        "db": "nucleotide", "term": query, "retmode": "json", "retmax": str(max_records),
        "email": EMAIL, "tool": TOOL,
    }
    text = eutils_get("esearch.fcgi", params)
    ids = json.loads(text)["esearchresult"]["idlist"]
    time.sleep(DELAY)
    if not ids:
        return []
    params2 = {
        "db": "nucleotide", "id": ",".join(ids), "rettype": "fasta",
        "retmode": "text", "email": EMAIL, "tool": TOOL,
    }
    fasta = eutils_get("efetch.fcgi", params2)
    time.sleep(DELAY)
    records = []
    header = None
    seq = []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq)))
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if header:
        records.append((header, "".join(seq)))
    return records


def build_q2(species_rows, refs):
    records = []
    meta = []
    seen_seqs = set()
    for row in species_rows:
        taxid = row["taxid"]
        category = row["category"]
        fetched = fetch_records(taxid, category)
        kept = 0
        for header, seq in fetched:
            if kept >= Q2_MAX_PER_SPECIES or len(seq) < 100:
                continue
            norm = normalize(seq)
            if norm in refs or norm in seen_seqs:
                continue
            accession = re.match(r"(\S+)", header)
            accession = accession.group(1) if accession else ""
            if not marker_ok(header, category):
                continue
            seen_seqs.add(norm)
            qid = f"e{len(records) + 1:05d}"
            records.append({
                "qid": qid, "taxid": taxid, "species": row["species"],
                "category": category, "accession": accession,
                "header": f"{qid}|{taxid}|{row['species']}|{category}|{accession}",
                "seq": seq,
            })
            meta.append({
                "qid": qid, "taxid": taxid, "species": row["species"],
                "category": category, "accession": accession,
                "title": header, "seq_len": len(seq),
            })
            kept += 1
    return records, meta


def write_fasta(records, path):
    with open(path, "w", encoding="utf-8") as fh:
        for record in records:
            fh.write(">" + record["header"] + "\n")
            seq = record["seq"]
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_meta(meta, path, fields):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(fields) + "\n")
        for item in meta:
            fh.write("\t".join(str(item.get(f, "")) for f in fields) + "\n")


def build_q3(refs):
    src_fasta = os.path.join(
        ROOT, "Figshare", "docs", "validation", "classification_benchmark",
        "negative_queries.fasta",
    )
    src_meta = os.path.join(
        ROOT, "Figshare", "docs", "validation", "classification_benchmark",
        "negative_queries_meta.json",
    )
    records = []
    header = None
    seq = []
    with open(src_fasta, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            records.append((header, "".join(seq)))
    meta_map = {}
    if os.path.exists(src_meta):
        for item in json.load(open(src_meta, encoding="utf-8")):
            meta_map[item["qid"]] = item
    kept = []
    overlaps = 0
    for header, seq in records:
        id_key = header.split(" ", 1)[0]
        qid = id_key.split("|", 1)[0]
        if normalize(seq) in refs:
            overlaps += 1
            continue
        meta = meta_map.get(id_key, {})
        taxid = meta.get("taxid", id_key.split("|")[1] if "|" in id_key else "")
        species = meta.get("species", "")
        category = meta.get("category", "")
        accession = meta.get("accession", "")
        new_header = f"{qid}|{taxid}|{species}|{category}|{accession}"
        kept.append({"qid": qid, "header": new_header, "seq": seq, "meta": meta})
    return kept, overlaps


def parse_archived_fasta(path, source):
    records = []
    header = None
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            records.append((header, "".join(seq)))
    out = []
    for header, seq in records:
        parts = header.split("|", 3)
        taxid = parts[0].strip()
        species = parts[1].strip() if len(parts) > 1 else ""
        rest = parts[2].strip() if len(parts) > 2 else ""
        category = rest.split()[0] if rest else ""
        category = "oomycetes" if "Oomycota" in header else category
        accession = re.search(r"[A-Z]{1,2}_?\d+(?:\.\d+)?", header)
        out.append({
            "taxid": taxid, "species": species, "category": category,
            "accession": accession.group(0) if accession else "",
            "header": header, "seq": seq, "source": source,
        })
    return out


def build_q4(species_rows, refs, q2_records):
    allowed = {row["taxid"] for row in species_rows}
    name_to_taxid = {}
    for row in species_rows:
        name_to_taxid.setdefault(
            re.sub(r"[^a-z0-9]+", "", row["species"].lower()), row["taxid"]
        )
        name_to_taxid.setdefault(
            re.sub(r"[^a-z0-9]+", "", row.get("full_name", "").lower()),
            row["taxid"],
        )
    candidates = []
    candidates += parse_archived_fasta(
        os.path.join(ROOT, "Figshare", "docs", "validation", "external",
                     "silva_unite_cross_queries.fasta"),
        "silva_unite",
    )
    candidates += parse_archived_fasta(
        os.path.join(ROOT, "Figshare", "docs", "validation", "external",
                     "refseq_queries.fasta"),
        "refseq",
    )
    kept = []
    seen_seqs = {
        normalize(record["seq"]) for record in q2_records
    }
    per_taxid = {}
    for record in candidates:
        taxid = record["taxid"]
        if taxid not in allowed:
            taxid = name_to_taxid.get(
                re.sub(r"[^a-z0-9]+", "", record["species"].lower())
            )
        if not taxid:
            continue
        if per_taxid.get(taxid, 0) >= Q4_MAX_PER_SPECIES:
            continue
        norm = normalize(record["seq"])
        if norm in refs or norm in seen_seqs:
            continue
        seen_seqs.add(norm)
        per_taxid[taxid] = per_taxid.get(taxid, 0) + 1
        qid = f"x{len(kept) + 1:05d}"
        kept.append({
            "qid": qid,
            "taxid": taxid,
            "species": record["species"],
            "category": record["category"],
            "accession": record["accession"],
            "header": f"{qid}|{taxid}|{record['species']}|{record['category']}|{record['accession']}",
            "seq": record["seq"],
            "source": record["source"],
            "title": record["header"],
        })

    for row in species_rows:
        taxid = row["taxid"]
        category = row["category"]
        if per_taxid.get(taxid, 0) >= Q4_MAX_PER_SPECIES:
            continue
        fetched = fetch_records(
            taxid, category, max_records=Q4_FETCH_MAX,
            date_expr='"1900/01/01"[PDAT] : "2026/08/25"[PDAT]',
        )
        for header, seq in fetched:
            if per_taxid.get(taxid, 0) >= Q4_MAX_PER_SPECIES:
                break
            norm = normalize(seq)
            if norm in refs or norm in seen_seqs or len(seq) < 100:
                continue
            accession_match = re.match(r"(\S+)", header)
            accession = accession_match.group(1) if accession_match else ""
            if not marker_ok(header, category):
                continue
            seen_seqs.add(norm)
            per_taxid[taxid] = per_taxid.get(taxid, 0) + 1
            qid = f"x{len(kept) + 1:05d}"
            kept.append({
                "qid": qid,
                "taxid": taxid,
                "species": row["species"],
                "category": category,
                "accession": accession,
                "header": f"{qid}|{taxid}|{row['species']}|{category}|{accession}",
                "seq": seq,
                "source": "ncbi",
                "title": header,
            })
    return kept


def main():
    with open(SPECIES_LIST, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    species_rows = [
        r for r in rows
        if r["sequence_status"] == "present" and r["taxid"] != "NOT_VERIFIED"
    ]
    refs = load_reference_seqs()
    os.makedirs(OUT_DIR, exist_ok=True)

    # Q2
    q2_fasta_path = os.path.join(OUT_DIR, "external_positives.fasta")
    q2_meta_path = os.path.join(OUT_DIR, "external_positives_meta.tsv")
    if os.environ.get("REUSE_QUERY_SETS") == "1" and os.path.exists(q2_fasta_path):
        q2_records = []
        header = None
        seq = []
        with open(q2_fasta_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if header:
                        parts = header.split("|")
                        q2_records.append({
                            "qid": parts[0], "taxid": parts[1],
                            "species": parts[2], "category": parts[3],
                            "accession": parts[4], "header": header,
                            "seq": "".join(seq),
                        })
                    header = line[1:].strip()
                    seq = []
                else:
                    seq.append(line.strip())
            if header:
                parts = header.split("|")
                q2_records.append({
                    "qid": parts[0], "taxid": parts[1],
                    "species": parts[2], "category": parts[3],
                    "accession": parts[4], "header": header,
                    "seq": "".join(seq),
                })
        q2_meta = []
        with open(q2_meta_path, encoding="utf-8") as fh:
            for item in csv.DictReader(fh, delimiter="\t"):
                q2_meta.append(item)
    else:
        q2_records, q2_meta = build_q2(species_rows, refs)
    write_fasta(q2_records, q2_fasta_path)
    write_meta(q2_meta, q2_meta_path,
               ["qid", "taxid", "species", "category", "accession", "title", "seq_len"])
    print("Q2 external positives:", len(q2_records),
          dict(Counter(r["category"] for r in q2_records)))

    # Q3
    q3_records, q3_overlaps = build_q3(refs)
    write_fasta(q3_records, os.path.join(OUT_DIR, "negatives.fasta"))
    write_meta(
        [{"qid": r["qid"], "header": r["header"], "seq_len": len(r["seq"]),
          **{k: v for k, v in r["meta"].items() if k != "seq"}}
         for r in q3_records],
        os.path.join(OUT_DIR, "negatives_meta.tsv"),
        ["qid", "header", "seq_len", "taxid", "species", "category", "accession", "title"],
    )
    print("Q3 negatives:", len(q3_records), "overlaps removed:", q3_overlaps,
          dict(Counter(r["meta"].get("category", "?") for r in q3_records)))

    # Q4
    q4_records = build_q4(species_rows, refs, q2_records)
    write_fasta(q4_records, os.path.join(OUT_DIR, "cross_db_queries.fasta"))
    write_meta(
        [{"qid": r["qid"], "taxid": r["taxid"], "species": r["species"],
          "category": r["category"], "accession": r["accession"],
          "source": r["source"], "title": r["header"], "seq_len": len(r["seq"])}
         for r in q4_records],
        os.path.join(OUT_DIR, "cross_db_queries_meta.tsv"),
        ["qid", "taxid", "species", "category", "accession", "source", "title", "seq_len"],
    )
    print("Q4 cross-db queries:", len(q4_records),
          dict(Counter(r["category"] for r in q4_records)))
    print("query sets written to", OUT_DIR)


if __name__ == "__main__":
    main()
