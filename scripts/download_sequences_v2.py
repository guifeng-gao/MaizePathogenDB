#!/usr/bin/env python3
"""Download marker genes using direct HTTP requests (no BioPython)."""
import requests, json, time, os, sys

BASE_DIR = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
SEQ_DIR = os.path.join(BASE_DIR, "sequences")
TAX_DIR = os.path.join(BASE_DIR, "taxonomy")
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"

def esearch(db, term, retmax=5):
    """Search NCBI Entrez."""
    params = {"db": db, "term": term, "retmax": retmax, "retmode": "json",
              "email": EMAIL, "tool": "maize_pathogen_db", "sort": "sequence_length"}
    r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
    if r.status_code == 429:
        time.sleep(5)
        r = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", []), int(r.json().get("esearchresult", {}).get("count", 0))

def efetch(db, ids, rettype="fasta"):
    """Fetch sequences."""
    params = {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": "text",
              "email": EMAIL, "tool": "maize_pathogen_db"}
    r = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=60)
    if r.status_code == 429:
        time.sleep(5)
        r = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=60)
    r.raise_for_status()
    return r.text

def efetch_taxonomy_xml(taxid):
    """Fetch taxonomy lineage."""
    params = {"db": "taxonomy", "id": taxid, "retmode": "xml",
              "email": EMAIL, "tool": "maize_pathogen_db"}
    r = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=30)
    if r.status_code == 429:
        time.sleep(5)
        r = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=30)
    r.raise_for_status()
    import xml.etree.ElementTree as ET
    root = ET.fromstring(r.text)
    for taxon in root.findall(".//Taxon"):
        lineage = taxon.findtext("Lineage", "")
        sci = taxon.findtext("ScientificName", "")
        return f"{lineage}; {sci}"
    return ""

# Load TaxIDs
with open(os.path.join(TAX_DIR, "all_taxids.json")) as f:
    records = json.load(f)

print(f"Downloading marker genes for {len(records)} maize pathogens")
print(f"Strategy: Bacteria=16S, Fungi=ITS, Viruses=genome")
print("=" * 60)

all_seqs = {}
stats = {"bacteria": 0, "viruses": 0, "fungi": 0}
tax_table = []
downloaded_taxids = set()

for i, rec in enumerate(records):
    taxid = rec["tax_id"]
    sp = rec["species"]
    cat = rec["category"]
    disease = rec["disease_cn"]
    
    print(f"[{i+1}/{len(records)}] {cat}: {sp[:55]} (txid{taxid})")
    sys.stdout.flush()
    
    # Build query based on category
    if cat == "bacteria":
        query = f'txid{taxid}[Organism] AND ("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title]) AND 1000:2000[SLEN]'
    elif cat == "fungi":
        query = f'txid{taxid}[Organism] AND ("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title] OR "18S ribosomal RNA"[Title]) AND 400:2000[SLEN]'
    elif cat == "viruses":
        query = f'txid{taxid}[Organism] AND ("complete genome"[Title] OR "segment"[Title] OR "polyprotein"[Gene] OR "coat protein"[Gene]) NOT ("clone"[Title] OR "vector"[Title] OR "synthetic"[Title])'
    
    try:
        ids, count = esearch("nucleotide", query, retmax=3)
        print(f"  Found {count} records, downloading top {len(ids)}...")
        
        if ids:
            fasta_text = efetch("nucleotide", ids)
            
            n_seqs = 0
            for record_block in fasta_text.strip().split("\n\n"):
                if not record_block.strip():
                    continue
                lines = record_block.strip().split("\n")
                if len(lines) < 2:
                    continue
                header = lines[0].strip()
                seq = "".join(l.strip() for l in lines[1:])
                
                if len(seq) >= 100:
                    clean_hdr = f">{taxid}|{sp[:80]}|{cat}|{disease[:80]}"
                    full_hdr = f"{clean_hdr} {header[1:]}"
                    all_seqs[full_hdr] = seq
                    n_seqs += 1
                    downloaded_taxids.add(taxid)
            
            if n_seqs > 0:
                lineage = efetch_taxonomy_xml(taxid)
                tax_table.append({
                    "tax_id": taxid, "species": sp, "category": cat,
                    "disease_cn": disease, "ncbi_lineage": lineage,
                    "seq_count": n_seqs,
                })
                stats[cat] += 1
                print(f"  ✓ Downloaded {n_seqs} sequences, total seqs: {len(all_seqs)}")
            else:
                print(f"  ⚠ No valid sequences after filtering")
        else:
            print(f"  ⚠ No records found")
    
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
    
    # Respect NCBI rate limits: 3/sec without key, 10/sec with key
    time.sleep(0.5)

# Save everything
print(f"\n{'=' * 60}")
print(f"DOWNLOAD COMPLETE")
print(f"  Total sequences: {len(all_seqs)}")
print(f"  Species covered: {len(downloaded_taxids)}/{len(records)}")
for cat in ["bacteria", "viruses", "fungi"]:
    n = len([r for r in records if r["category"] == cat])
    print(f"  {cat}: {stats[cat]}/{n}")

# Write FASTA
fasta_path = os.path.join(SEQ_DIR, "maize_pathogens_all.fasta")
with open(fasta_path, "w") as f:
    for hdr, seq in all_seqs.items():
        f.write(hdr + "\n")
        for j in range(0, len(seq), 80):
            f.write(seq[j:j+80] + "\n")
print(f"\n  FASTA: {fasta_path}")

# Per-category FASTA files
for cat in ["bacteria", "viruses", "fungi"]:
    cat_seqs = {h: s for h, s in all_seqs.items() if f"|{cat}|" in h}
    cat_path = os.path.join(SEQ_DIR, cat, f"maize_pathogens_{cat}.fasta")
    with open(cat_path, "w") as f:
        for hdr, seq in cat_seqs.items():
            f.write(hdr + "\n")
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + "\n")
    print(f"  {cat}: {len(cat_seqs)} seqs -> {cat_path}")

# Save taxonomy table
with open(os.path.join(TAX_DIR, "taxonomy_table.json"), "w") as f:
    json.dump(tax_table, f, indent=2, ensure_ascii=False)

with open(os.path.join(BASE_DIR, "docs", "download_stats.json"), "w") as f:
    json.dump({"categories": stats, "total_species": len(downloaded_taxids), "total_sequences": len(all_seqs)}, f, indent=2)

print("\n✓ All output files saved.")
