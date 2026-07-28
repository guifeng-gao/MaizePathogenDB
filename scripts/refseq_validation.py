#!/usr/bin/env python3
"""Plan C: Validate against NCBI RefSeq reference sequences."""

import os, json, time, subprocess, requests
from collections import defaultdict

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
OUT = os.path.join(BASE, "docs", "validation", "external")
os.makedirs(OUT, exist_ok=True)

BLAST = "/tmp/ncbi-blast-2.17.0+/bin/blastn"
BLST = os.path.join(BASE, "blast_db")
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"

with open(os.path.join(BASE, "taxonomy", "all_taxids.json")) as f:
    records = json.load(f)

print("=" * 60)
print("PLAN C: NCBI RefSeq Reference Sequence Validation")
print("=" * 60)

refseq_seqs = {}
downloaded = 0

for i, rec in enumerate(records):
    taxid = rec["tax_id"]
    cat = rec["category"]
    sp = rec["species"]
    
    if cat == "bacteria":
        query = f'txid{taxid}[Organism] AND ("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title]) AND refseq[filter] AND 1000:2000[SLEN]'
    elif cat == "fungi":
        query = f'txid{taxid}[Organism] AND ("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS"[Title]) AND refseq[filter] AND 200:2000[SLEN]'
    elif cat == "viruses":
        query = f'txid{taxid}[Organism] AND ("complete genome"[Title]) AND refseq[filter]'
    else:
        continue
    
    try:
        params = {"db": "nucleotide", "term": query, "retmax": 3, "retmode": "json",
                  "email": EMAIL, "tool": "maize_pathogen_db", "sort": "sequence_length"}
        r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
        if r.status_code == 429:
            time.sleep(3)
            r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        
        if not ids:
            continue
        
        params2 = {"db": "nucleotide", "id": ",".join(ids[:2]), "rettype": "fasta", "retmode": "text",
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
                clean_hdr = f">{taxid}|{sp[:50]}|{cat} {header[1:]}"
                refseq_seqs[clean_hdr] = {"seq": seq, "taxid": taxid, "cat": cat}
                downloaded += 1
        
        print(f"  [{downloaded}] {cat}: {sp[:50]} -> RefSeq records found")
    
    except:
        pass
    
    if (i + 1) % 30 == 0:
        print(f"  Progress: {i+1}/{len(records)} TaxIDs, {downloaded} RefSeq seqs")
    time.sleep(0.35)

# Save RefSeq queries
refseq_fasta = os.path.join(OUT, "refseq_queries.fasta")
with open(refseq_fasta, "w") as f:
    for hdr, info in refseq_seqs.items():
        f.write(hdr + "\n")
        s = info["seq"]
        for j in range(0, len(s), 80):
            f.write(s[j:j+80] + "\n")

print(f"\nDownloaded {downloaded} RefSeq sequences from {len(set(info['taxid'] for info in refseq_seqs.values()))} species")

# BLAST against database
print(f"\nBLASTing RefSeq queries against MaizePathogenDB...")
db_path = os.path.join(BLST, "maize_pathogens_all")
cmd = f'"{BLAST}" -query "{refseq_fasta}" -db "{db_path}" -outfmt "6 qseqid sseqid pident length qcovs" -max_target_seqs 1 -num_threads 4'
proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)

blast_hits = []
for line in proc.stdout.strip().split("\n"):
    f = line.strip().split("\t")
    if len(f) >= 5:
        blast_hits.append({"qseqid": f[0], "sseqid": f[1], "pident": float(f[2]), "length": int(f[3]), "qcovs": float(f[4])})

print(f"  {len(blast_hits)} BLAST hits")

# Analyze using TaxID map
taxid_to_cat = {r["tax_id"]: r["category"] for r in records}
cats = {"bacteria": {"correct": 0, "total": 0}, "viruses": {"correct": 0, "total": 0}, "fungi": {"correct": 0, "total": 0}}

for h in blast_hits:
    q_tax = h["qseqid"].split("|")[0].replace(">", "").strip()
    if not q_tax.isdigit():
        continue
    cat = taxid_to_cat.get(q_tax, "unknown")
    if cat not in cats:
        continue
    
    s_tax = h["sseqid"].split(".")[0]
    cats[cat]["total"] += 1
    if q_tax == s_tax:
        cats[cat]["correct"] += 1

print(f"\n{'=' * 60}")
print("REFSEQ VALIDATION RESULTS")
print(f"{'=' * 60}")
for cat in ["bacteria", "viruses", "fungi"]:
    d = cats[cat]
    if d["total"] > 0:
        acc = d["correct"] / d["total"] * 100
        print(f"  {cat}: {d['correct']}/{d['total']} ({acc:.1f}%)")

total_c = sum(d["correct"] for d in cats.values())
total_t = sum(d["total"] for d in cats.values())
print(f"  OVERALL: {total_c}/{total_t} ({total_c/total_t*100:.1f}%)")

# Save
refseq_results = {"validation_type": "refseq", "results": {}}
for cat in cats:
    d = cats[cat]
    if d["total"] > 0:
        refseq_results["results"][cat] = {"n_queries": d["total"], "n_correct": d["correct"], "accuracy_pct": round(d["correct"]/d["total"]*100, 1)}

with open(os.path.join(OUT, "refseq_validation.json"), "w") as f:
    json.dump(refseq_results, f, indent=2)

print(f"\n✓ Results saved to {OUT}/refseq_validation.json")
