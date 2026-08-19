#!/usr/bin/env python3
"""External validation: BLAST 2020-2026 NCBI submissions against MaizePathogenDB."""

import json
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict

import requests
from Bio import SeqIO

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
BLAST_BIN = "/Users/gfgao/Desktop/blacksoil_metaG/tools/ncbi-blast-2.17.0+/bin/blastn"
DB_PATH = os.path.join(BASE, "blast_db", "maize_pathogens_all")
OUT = os.path.join(BASE, "docs", "validation", "external")
OUT_FASTA = os.path.join(OUT, "external_queries.fasta")
OUT_META = os.path.join(OUT, "external_queries_meta.json")
OUT_RAW = os.path.join(OUT, "blast_results_raw.json")
OUT_FINAL = os.path.join(OUT, "external_validation_final.json")
OUT_SUMMARY = os.path.join(OUT, "external_validation_2020_2026_summary.md")

CATALOG = json.load(open("/Users/gfgao/Desktop/blacksoil_metaG/Figshare/taxonomy.json"))
REFERENCE_FASTA = os.path.join(BASE, "sequences", "maize_pathogens_all.fasta")

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"
TOOL = "maize_pathogen_db"
NCBI_DELAY = 0.4
MAX_PER_TAXON = 3
MIN_LEN = 100
START_DATE = "2020/01/01"

GENE_QUERIES = {
    "bacteria": '("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title])',
    "viruses": '("complete genome"[Title] OR "polyprotein"[Gene])',
    "fungi": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])',
    "oomycetes": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])',
}

CAT_ORDER = ["bacteria", "viruses", "fungi", "oomycetes"]


def normalized_seq(seq):
    return re.sub(r"[^ACGTNacgtn]", "", seq).upper()


def load_reference_seqs():
    refs = set()
    for rec in SeqIO.parse(REFERENCE_FASTA, "fasta"):
        refs.add(normalized_seq(str(rec.seq)))
    return refs


def fetch_new_sequences(records):
    external = []
    ref_seqs = load_reference_seqs()
    idx = 0

    for rec in records:
        taxid = str(rec["taxid"])
        species = rec["species"]
        cat = {"virus": "viruses"}.get(rec["category"].lower(), rec["category"].lower())
        if cat not in GENE_QUERIES:
            continue
        if taxid == "NOT_VERIFIED":
            continue

        query = (
            f'txid{taxid}[Organism] AND {GENE_QUERIES[cat]} '
            f'AND ("{START_DATE}"[PDAT] : "3000/12/31"[PDAT])'
        )
        params = {
            "db": "nucleotide", "term": query, "retmax": MAX_PER_TAXON,
            "retmode": "json", "email": EMAIL, "tool": TOOL,
            "sort": "pub_date",
        }
        try:
            r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
            if r.status_code == 429:
                time.sleep(3)
                r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
        except Exception as e:
            print(f"  esearch fail {taxid} {species[:40]}: {e}", flush=True)
            ids = []

        if ids:
            params2 = {
                "db": "nucleotide", "id": ",".join(ids),
                "rettype": "fasta", "retmode": "text",
                "email": EMAIL, "tool": TOOL,
            }
            try:
                r2 = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params2, timeout=60)
                if r2.status_code == 429:
                    time.sleep(3)
                    r2 = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params2, timeout=60)
                r2.raise_for_status()
            except Exception as e:
                print(f"  efetch fail {taxid} {species[:40]}: {e}", flush=True)
                r2 = None

            if r2 is not None:
                for block in r2.text.strip().split("\n\n"):
                    lines = [ln.strip() for ln in block.split("\n") if ln.strip()]
                    if len(lines) < 2 or not lines[0].startswith(">"):
                        continue
                    seq = "".join(lines[1:])
                    if len(seq) < MIN_LEN:
                        continue
                    norm = normalized_seq(seq)
                    if norm in ref_seqs:
                        continue
                    idx += 1
                    qid = f"q{idx:04d}|{taxid}|{cat}"
                    original = lines[0][1:]
                    accession = original.split()[0] if original else ""
                    external.append({
                        "qid": qid, "taxid": taxid, "species": species,
                        "category": cat, "accession": accession,
                        "title": original, "seq": seq,
                    })
        print(f"  [{len(external)}] {taxid}: {species[:50]} -> {len(ids)} ids",
              flush=True)
        time.sleep(NCBI_DELAY)

    return external


def save_fasta(external):
    with open(OUT_FASTA, "w") as f:
        for rec in external:
            f.write(f">{rec['qid']} {rec['title']}\n")
            for i in range(0, len(rec["seq"]), 80):
                f.write(rec["seq"][i:i + 80] + "\n")
    meta = [
        {"qid": r["qid"], "taxid": r["taxid"], "species": r["species"],
         "category": r["category"], "accession": r["accession"],
         "title": r["title"], "seq_len": len(r["seq"])}
        for r in external
    ]
    with open(OUT_META, "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)


def run_blast():
    cmd = (
        f'"{BLAST_BIN}" -query "{OUT_FASTA}" -db "{DB_PATH}" '
        f'-outfmt "6 qseqid sseqid pident length qcovs bitscore evalue" '
        f'-max_target_seqs 3 -num_threads 4 -evalue 1e-5'
    )
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
    hits = []
    for line in proc.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 7:
            hits.append({
                "qseqid": parts[0], "sseqid": parts[1],
                "pident": float(parts[2]), "length": int(parts[3]),
                "qcovs": float(parts[4]), "bitscore": float(parts[5]),
                "evalue": parts[6],
            })
    with open(OUT_RAW, "w") as f:
        json.dump(hits, f, indent=2)
    return hits


def analyze(hits, external):
    meta = {r["qid"]: r for r in external}
    by_query = defaultdict(list)
    for h in hits:
        by_query[h["qseqid"]].append(h)

    results = []
    for rec in external:
        qid = rec["qid"]
        top = by_query[qid][0] if by_query.get(qid) else None
        if top is None:
            results.append({**rec, "hit": False, "top_taxid": None,
                            "top_sseqid": "", "pident": 0.0, "qcovs": 0.0,
                            "correct": False})
            continue
        top_taxid = top["sseqid"].split("|")[0]
        results.append({**rec, "hit": True, "top_taxid": top_taxid,
                        "top_sseqid": top["sseqid"], "pident": top["pident"],
                        "qcovs": top["qcovs"], "correct": top_taxid == rec["taxid"]})

    per_cat = {}
    for cat in CAT_ORDER:
        rows = [r for r in results if r["category"] == cat]
        n = len(rows)
        correct = sum(1 for r in rows if r.get("correct"))
        per_cat[cat] = {
            "n_queries": n,
            "n_correct": correct,
            "accuracy": round(correct / n * 100, 1) if n else 0.0,
        }

    total = len(results)
    total_correct = sum(1 for r in results if r.get("correct"))
    final = {
        "validation_type": "external",
        "description": f"Sequences submitted to NCBI since {START_DATE}, not present in MaizePathogenDB",
        "run_date": "2026-08-19",
        "total_external_queries": total,
        "total_taxa_with_new_seqs": len({r["taxid"] for r in results}),
        "results": per_cat,
        "overall": {
            "n_queries": total,
            "n_correct": total_correct,
            "accuracy": round(total_correct / total * 100, 1) if total else 0.0,
        },
        "details": [
            {k: r[k] for k in ("qid", "taxid", "species", "category", "accession",
                               "hit", "top_taxid", "top_sseqid", "pident", "qcovs", "correct")}
            for r in results
        ],
    }
    with open(OUT_FINAL, "w") as f:
        json.dump(final, f, indent=2, ensure_ascii=False)
    return final


def write_summary(final):
    ref_species = defaultdict(set)
    for rec in SeqIO.parse(REFERENCE_FASTA, "fasta"):
        tid = rec.id.split("|")[0]
        sp = rec.description.split("|", 2)[1] if rec.description.count("|") >= 2 else rec.description
        ref_species[tid].add(sp.strip())

    def genus(sp):
        m = re.match(r"([a-zA-Z]+)", sp.strip())
        return m.group(1).lower() if m else ""

    genus_ok = Counter()
    genus_n = Counter()
    for d in final["details"]:
        cat = d["category"]
        genus_n[cat] += 1
        qgenus = genus(d["species"])
        top_genera = {genus(x) for x in ref_species.get(d.get("top_taxid"), set())}
        if qgenus and qgenus in top_genera:
            genus_ok[cat] += 1

    per_taxon = {}
    for d in final["details"]:
        per_taxon.setdefault(d["taxid"], []).append(d)
    one_strict = {cat: [0, 0] for cat in CAT_ORDER}
    for rows in per_taxon.values():
        first = rows[0]
        one_strict[first["category"]][1] += 1
        if first.get("correct"):
            one_strict[first["category"]][0] += 1

    lines = [
        f"# External Validation Summary: NCBI New Sequences since {START_DATE}",
        "",
        "Run date: 2026-08-19",
        "",
        "Strict criterion: top-1 BLAST hit has the same NCBI TaxID as the query species.",
        "",
        "| Category | n | Correct | Accuracy |",
        "|---|---:|---:|---:|",
    ]
    for cat in CAT_ORDER:
        r = final["results"][cat]
        lines.append(f"| {cat} | {r['n_queries']} | {r['n_correct']} | {r['accuracy']:.1f}% |")
    o = final["overall"]
    lines.append(f"| **Overall** | **{o['n_queries']}** | **{o['n_correct']}** | **{o['accuracy']:.1f}%** |")
    lines.append("")
    lines.append("Genus-level criterion (top-1 hit shares query genus):")
    lines.append("")
    lines.append("| Category | n | Genus correct | Genus accuracy |")
    lines.append("|---|---:|---:|---:|")
    for cat in CAT_ORDER:
        n = genus_n.get(cat, 0)
        g = genus_ok.get(cat, 0)
        lines.append(f"| {cat} | {n} | {g} | {g / n * 100:.1f}% |" if n else f"| {cat} | 0 | 0 | - |")
    lines.append("")
    lines.append("One sequence per taxon (strict):")
    lines.append("")
    lines.append("| Category | n | Correct | Accuracy |")
    lines.append("|---|---:|---:|---:|")
    for cat in CAT_ORDER:
        n = one_strict[cat][1]
        c = one_strict[cat][0]
        lines.append(f"| {cat} | {n} | {c} | {c / n * 100:.1f}% |" if n else f"| {cat} | 0 | 0 | - |")
    lines.append("")
    lines.append("Note: this run uses the current 225-species catalog and only sequences not already in "
                 "MaizePathogenDB. The previous 94.1% external result used a different catalog/test set and a "
                 "more lenient correctness judgment, so the two numbers are not directly comparable.")
    lines.append("")
    errors = [d for d in final["details"] if not d.get("correct")]
    lines.append(f"Errors ({len(errors)}):")
    for d in errors:
        lines.append(f"- {d['species']} ({d['category']}, pident={d['pident']:.1f}%, top={d['top_sseqid']})")
    with open(OUT_SUMMARY, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    records = CATALOG
    print(f"Catalog: {len(records)} species")
    reuse = "--reuse" in sys.argv
    if reuse and os.path.exists(OUT_META) and os.path.exists(OUT_FASTA):
        meta = json.load(open(OUT_META))
        seqs = {rec.id.split()[0]: str(rec.seq)
                for rec in SeqIO.parse(OUT_FASTA, "fasta")}
        external = [{**m, "seq": seqs[m["qid"]]} for m in meta if m["qid"] in seqs]
        print(f"Reusing {len(external)} downloaded external sequences")
    else:
        print(f"Downloading NCBI sequences since {START_DATE}...")
        external = fetch_new_sequences(records)
        print(f"\nDownloaded {len(external)} external sequences")
    if not external:
        print("No new sequences found.")
        return

    save_fasta(external)
    print(f"Saved queries: {OUT_FASTA}")
    print("Running BLAST against MaizePathogenDB...")
    hits = run_blast()
    print(f"Raw BLAST hits: {len(hits)}")

    final = analyze(hits, external)
    write_summary(final)
    print("\nFinal:")
    for cat in CAT_ORDER:
        r = final["results"][cat]
        print(f"  {cat}: {r['n_correct']}/{r['n_queries']} ({r['accuracy']:.1f}%)")
    o = final["overall"]
    print(f"  overall: {o['n_correct']}/{o['n_queries']} ({o['accuracy']:.1f}%)")


if __name__ == "__main__":
    main()
