#!/usr/bin/env python3
"""External validation: BLAST newly submitted sequences against our database."""

import os, json, sys, time, subprocess
import requests

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
BLST = os.path.join(BASE, "blast_db")
OUT = os.path.join(BASE, "docs", "validation", "external")
OUT_FASTA = os.path.join(OUT, "external_queries")
os.makedirs(OUT, exist_ok=True)
os.makedirs(OUT_FASTA, exist_ok=True)

BLAST = "/tmp/ncbi-blast-2.17.0+/bin/blastn"
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"

# Load TaxIDs
with open(os.path.join(BASE, "taxonomy", "all_taxids.json")) as f:
    records = json.load(f)

# Download new sequences (2025-2026) per TaxID
print("=" * 60)
print("PHASE 1: Downloading external validation sequences (2025-2026)")
print("=" * 60)

external_seqs = {}  # clean_header -> sequence
downloaded = 0
skipped = 0

for rec in records:
    taxid = rec["tax_id"]
    cat = rec["category"]
    sp = rec["species"]
    
    # Query: sequences from THIS taxid, published 2025-2026, matching marker genes
    if cat == "bacteria":
        gene_query = '("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title])'
    elif cat == "fungi":
        gene_query = '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])'
    elif cat == "viruses":
        gene_query = '("complete genome"[Title] OR "polyprotein"[Gene])'
    
    query = f'txid{taxid}[Organism] AND {gene_query} AND ("2025/01/01"[PDAT] : "3000/12/31"[PDAT])'
    
    try:
        params = {"db": "nucleotide", "term": query, "retmax": 3, "retmode": "json",
                  "email": EMAIL, "tool": "maize_pathogen_db", "sort": "pub_date"}
        r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
        if r.status_code == 429:
            time.sleep(3)
            r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        
        if not ids:
            skipped += 1
            continue
        
        # Fetch sequences
        params2 = {"db": "nucleotide", "id": ",".join(ids), "rettype": "fasta", "retmode": "text",
                   "email": EMAIL, "tool": "maize_pathogen_db"}
        r2 = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params2, timeout=60)
        if r2.status_code == 429:
            time.sleep(3)
            r2 = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params2, timeout=60)
        r2.raise_for_status()
        
        for block in r2.text.strip().split("\n\n"):
            if not block.strip():
                continue
            lines = block.strip().split("\n")
            if len(lines) < 2:
                continue
            header = lines[0].strip()
            seq = "".join(l.strip() for l in lines[1:])
            
            if len(seq) >= 100:
                clean_hdr = f">{taxid}|{sp[:60]}|{cat}"
                full_hdr = f"{clean_hdr} {header[1:]}"
                external_seqs[full_hdr] = seq
                downloaded += 1
        
        print(f"  [{downloaded}] TaxID {taxid}: {sp[:50]} -> {len(ids)} new seqs")
        
    except Exception as e:
        print(f"  ✗ TaxID {taxid}: {str(e)[:60]}")
    
    time.sleep(0.4)

# Save external queries
ext_fasta = os.path.join(OUT_FASTA, "external_queries.fasta")
with open(ext_fasta, "w") as f:
    for hdr, seq in external_seqs.items():
        f.write(hdr + "\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

print(f"\nDownloaded {downloaded} new sequences from {len(records)-skipped}/{len(records)} TaxIDs")
print(f"Saved to: {ext_fasta}")

if downloaded == 0:
    print("\n⚠ No new sequences found. Maybe broaden date range?")
    print("Try manually: change date filter to 2020-2026")
    sys.exit(0)

# Phase 2: BLAST external against our database
print(f"\n{'=' * 60}")
print("PHASE 2: BLAST external sequences against MaizePathogenDB")
print("=" * 60)

blast_results = []
for db_name in ["maize_pathogens_all", "maize_pathogens_bacteria", "maize_pathogens_viruses", "maize_pathogens_fungi"]:
    db_path = os.path.join(BLST, db_name)
    if not os.path.exists(f"{db_path}.nsq"):
        continue
    
    print(f"  Running BLAST against {db_name}...")
    
    # Only BLAST relevant category
    cmd = f'"{BLAST}" -query "{ext_fasta}" -db "{db_path}" -outfmt "6 qseqid sseqid pident length qcovs" -max_target_seqs 3 -num_threads 4'
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    
    hits = []
    for line in proc.stdout.strip().split("\n"):
        f = line.strip().split("\t")
        if len(f) >= 5:
            hits.append({
                "qseqid": f[0], "sseqid": f[1],
                "pident": float(f[2]), "length": int(f[3]),
                "qcovs": float(f[4]), "db": db_name,
            })
    
    blast_results.extend(hits)
    print(f"    {len(hits)} hits ({db_name})")

# Save raw BLAST results
with open(os.path.join(OUT, "blast_results_raw.json"), "w") as f:
    json.dump(blast_results, f, indent=2)
print(f"  Saved {len(blast_results)} total BLAST hits")

# Phase 3: Analyze accuracy
print(f"\n{'=' * 60}")
print("PHASE 3: External validation accuracy analysis")
print("=" * 60)

from collections import defaultdict

# Group by query
queries = defaultdict(list)
for h in blast_results:
    queries[h["qseqid"]].append(h)

# Accuracy per category
results_by_cat = {"bacteria": [], "viruses": [], "fungi": []}
for qid, hits in queries.items():
    if not hits or "|" not in qid:
        continue
    
    # Parse query TaxID from header
    parts = qid.split("|")
    if len(parts) < 3:
        continue
    q_tax = parts[0].replace(">", "")
    cat = parts[2] if len(parts) > 2 else "unknown"
    
    if cat not in results_by_cat:
        continue
    
    # Top-1 hit TaxID
    top1 = hits[0]
    top1_tax = top1["sseqid"].split(".")[0]
    is_correct = (q_tax == top1_tax)
    
    results_by_cat[cat].append({
        "q_tax": q_tax, "top1_tax": top1_tax,
        "correct": is_correct, "pident": top1["pident"],
        "qcovs": top1["qcovs"], "db": top1.get("db", ""),
    })

# Calculate accuracy
print("\nExternal Validation Results (2025-2026 new sequences):")
print("-" * 50)
for cat in ["bacteria", "viruses", "fungi"]:
    data = results_by_cat[cat]
    n = len(data)
    if n == 0:
        print(f"  {cat}: No external queries available")
        continue
    correct = sum(1 for d in data if d["correct"])
    pident_values = [d["pident"] for d in data]
    acc = correct / n * 100
    
    print(f"  {cat}: {correct}/{n} correct ({acc:.1f}%)")
    print(f"          Mean identity: {sum(pident_values)/len(pident_values):.1f}%")

# Save validation summary
summary = {
    "validation_type": "external",
    "description": "Sequences submitted to NCBI in 2025-2026, not present in MaizePathogenDB v1.0",
    "total_external_queries": downloaded,
    "total_taxa_with_new_seqs": len(records) - skipped,
    "results": {
        cat: {
            "n_queries": len(data),
            "n_correct": sum(1 for d in data if d["correct"]),
            "accuracy": sum(1 for d in data if d["correct"]) / len(data) * 100 if data else 0,
            "mean_pident": sum(d["pident"] for d in data) / len(data) if data else 0,
        } for cat, data in results_by_cat.items()
    }
}

with open(os.path.join(OUT, "external_validation_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

print(f"\n✓ External validation complete")
print(f"  Results saved to: {OUT}/")
