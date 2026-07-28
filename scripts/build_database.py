#!/usr/bin/env python3
"""Build complete maize pathogen reference database including BLAST DB, taxonomy files, and stats."""

import os, json, subprocess, sys

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
SEQ_DIR = os.path.join(BASE, "sequences")
BLAST_DIR = os.path.join(BASE, "blast_db")
TAX_DIR = os.path.join(BASE, "taxonomy")
DOCS_DIR = os.path.join(BASE, "docs")

print("=" * 70)
print("BUILDING MAIZE PATHOGEN REFERENCE DATABASE")
print("=" * 70)

# 1. Verify FASTA file
fasta = os.path.join(SEQ_DIR, "maize_pathogens_all.fasta")
seq_count = 0
with open(fasta) as f:
    for line in f:
        if line.startswith(">"):
            seq_count += 1
print(f"\n1. FASTA: {fasta} - {seq_count} sequences")

# 2. Build BLAST nucleotide database
print(f"\n2. Building BLAST nucleotide database...")
makeblastdb = "makeblastdb"
try:
    result = subprocess.run([makeblastdb, "-version"], capture_output=True, text=True, timeout=10)
    print(f"   makeblastdb found: {result.stdout.strip().split(chr(10))[0]}")
except:
    print("   makeblastdb not in PATH, trying conda/brew...")
    for path in ["/usr/local/bin/makeblastdb", "/opt/homebrew/bin/makeblastdb", 
                 os.path.expanduser("~/anaconda3/bin/makeblastdb"),
                 os.path.expanduser("~/miniconda3/bin/makeblastdb")]:
        if os.path.exists(path):
            makeblastdb = path
            break
    else:
        print("   ⚠ makeblastdb not found - install with: conda install -c bioconda blast")
        makeblastdb = None

if makeblastdb:
    blast_out = os.path.join(BLAST_DIR, "maize_pathogens")
    cmd = [makeblastdb, "-in", fasta, "-dbtype", "nucl", "-out", blast_out, 
           "-title", "Maize Pathogen Reference Database v1.0"]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode == 0:
        print(f"   ✓ BLAST database built: {blast_out}")
    else:
        print(f"   ✗ Error: {result.stderr[:200]}")
else:
    print("   ⚠ Skipping BLAST database (makeblastdb not available)")

# 3. Create QIIME2/SINTAX taxonomy file
print(f"\n3. Creating taxonomy files for metagenomic classifiers...")

# Generate taxonomy.fasta (SINTAX format: >TaxID;tax=k__;p__;c__;o__;f__;g__;s__)
taxonomy_file = os.path.join(TAX_DIR, "maize_pathogens_taxonomy_sintax.fasta")
nodes_file = os.path.join(TAX_DIR, "maize_pathogens_nodes.dmp")
names_file = os.path.join(TAX_DIR, "maize_pathogens_names.dmp")

with open(fasta) as f_in, open(taxonomy_file, "w") as f_tax:
    for line in f_in:
        if line.startswith(">"):
            header = line.strip()
            parts = header.split("|")
            if len(parts) >= 4:
                taxid = parts[0][1:]  # remove >
                species = parts[1].strip()
                
                # Write SINTAX format taxonomy header
                f_tax.write(f"{header}\n")
                f_tax.write(f"{line}")

print(f"   ✓ SINTAX taxonomy: {taxonomy_file}")

# 4. Generate database statistics
print(f"\n4. Generating database statistics...")

stats = {}
categories = {"bacteria": 0, "viruses": 0, "fungi": 0}
seq_lengths = {"bacteria": [], "viruses": [], "fungi": []}
species_set = {"bacteria": set(), "viruses": set(), "fungi": set()}

with open(fasta) as f:
    current_taxid = None
    current_seq = ""
    for line in f:
        if line.startswith(">"):
            if current_seq and current_taxid:
                for cat in categories:
                    if f"|{cat}|" in current_taxid:
                        seq_lengths[cat].append(len(current_seq))
                        species_set[cat].add(current_taxid.split("|")[0])
                        categories[cat] += 1
            header = line.strip()
            current_taxid = header
            current_seq = ""
        else:
            current_seq += line.strip()
    # Last sequence
    if current_seq and current_taxid:
        for cat in categories:
            if f"|{cat}|" in current_taxid:
                seq_lengths[cat].append(len(current_seq))
                species_set[cat].add(current_taxid.split("|")[0])
                categories[cat] += 1

stats_report = []
stats_report.append(f"Maize Pathogen Reference Database v1.0 - Statistics")
stats_report.append(f"{'='*60}")
stats_report.append(f"Total sequences: {seq_count}")
stats_report.append(f"")
for cat in ["bacteria", "viruses", "fungi"]:
    n_seqs = categories[cat]
    n_spp = len(species_set[cat])
    avg_len = sum(seq_lengths[cat]) / len(seq_lengths[cat]) if seq_lengths[cat] else 0
    min_len = min(seq_lengths[cat]) if seq_lengths[cat] else 0
    max_len = max(seq_lengths[cat]) if seq_lengths[cat] else 0
    stats_report.append(f"{cat}: {n_seqs} sequences, {n_spp} species")
    stats_report.append(f"  Avg len: {avg_len:.0f} bp, Range: {min_len}-{max_len} bp")

stats_text = "\n".join(stats_report)
stats_path = os.path.join(DOCS_DIR, "database_statistics.txt")
with open(stats_path, "w") as f:
    f.write(stats_text)
print(stats_text)

# Save JSON stats
json_stats = {
    "database_version": "1.0",
    "creation_date": "2026-07-24",
    "total_sequences": seq_count,
    "categories": {cat: {"sequences": categories[cat], "species": len(species_set[cat]),
                          "avg_length": sum(seq_lengths[cat]) / len(seq_lengths[cat]) if seq_lengths[cat] else 0}
                   for cat in categories},
    "marker_genes": {
        "bacteria": "16S rRNA gene",
        "fungi": "ITS region (ITS1-5.8S-ITS2)",
        "viruses": "Complete genome / reference segments",
        "oomycetes": "ITS region / 18S rRNA gene",
    },
    "species_without_sequences": [
        "Spiroplasma kunkelii (TaxID 47834)",
        "Maize bushy stunt phytoplasma (TaxID 202462)", 
        "Pseudomonas syringae pv. coronafaciens (TaxID 235275)",
        "Xanthomonas vasicola pv. vasculorum (TaxID 325776)",
        "Fusarium temperatum (TaxID 1035347)",
        "Peronosclerospora sorghi (TaxID 230839)",
        "Peronosclerospora maydis (TaxID 886949)",
        "Sclerophthora macrospora (TaxID 467176)",
        "Sclerospora graminicola (TaxID 162130)",
    ]
}
with open(os.path.join(DOCS_DIR, "database_stats.json"), "w") as f:
    json.dump(json_stats, f, indent=2)

print(f"\n{'='*70}")
print(f"DATABASE BUILD COMPLETE")
print(f"{'='*70}")
print(f"Output directory: {BASE}")
print(f"  sequences/maize_pathogens_all.fasta - Merged FASTA")
print(f"  blast_db/maize_pathogens* - BLAST nucleotide database")
print(f"  taxonomy/maize_pathogens_taxonomy_sintax.fasta - SINTAX format")
print(f"  docs/database_statistics.txt - Statistics report")
print(f"  docs/database_stats.json - Machine-readable stats")
