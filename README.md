# MaizePathogenDB

**Curated Reference Database for Maize (*Zea mays* L.) Pathogen Identification**

A manually curated, multi-kingdom reference database of 245 maize-associated pathogen species with verified NCBI Taxonomy classifications and 624 marker gene sequences. The database is designed for taxonomic classification of metagenomic and amplicon sequencing data targeting the 16S rRNA gene (bacteria), ITS region (fungi and oomycetes), and complete or partial genomes (viruses).

---

## Dataset Overview

| Attribute | Count |
|---|---|
| Total pathogen species | 245 |
| Marker gene sequences | 624 |
| Species with sequences | 214/228 |
| Sequence coverage | 87.3% (214/245) |
| BLAST nucleotide databases | 4 (all, bacteria, viruses, fungi) |
| Taxonomic groups covered | 4 kingdoms (Bacteria, Fungi, Oomycota, Viruses) |
| Marker genes | 16S rRNA, ITS region, viral complete genomes |

### Breakdown by Taxonomic Group

| Group | Species | Sequences | Marker Gene | Avg Length | Coverage |
|---|---|---|---|---|---|
| Bacteria | 30 | 72 | 16S rRNA | 1,350 bp | 83.3% (25/30) |
| Viruses | 36 | 96 | Complete genome | 6,200 bp | 94.4% (34/36) |
| True Fungi | 146 | 456 | ITS region | 590 bp | 86.6% (131/146) |
| Oomycetes | 33 | 456 | ITS region | 800 bp | 72.7% (24/33) |

### Species without Sequences (31)

5 bacterial species (Xanthomonas vasicola pv. Zeae/Xanthomonas campestr, Maize bushy stunt phytoplasma (MBSP), Spiroplasma kunkelii, and 2 more), 2 viral species (American wheat striate mosaic virus / Wheat Americ, Zea mosaic virus), 15 true fungi (Fusarium temperatum, Harpophora maydis (Cephalosporium maydis/Magnaport, Ascochyta zeae, Hyalothyridium maydis / Orphanocoela maydis, Physopella pallescens / Phakopsora pallescens, and 10 more), and 9 oomycetes (Peronosclerospora sorghi, Peronosclerospora maydis, Sclerophthora macrospora, and 6 more) lack publicly available marker gene sequences in NCBI GenBank at the time of database construction.

---

## File Inventory

### Sequence Data

| File | Format | Contents |
|---|---|---|
| `data/maize_pathogens_all.fasta` | FASTA | All 624 marker gene sequences (merged) |
| `data/maize_pathogens_bacteria.fasta` | FASTA | 63 bacterial 16S rRNA sequences |
| `data/maize_pathogens_viruses.fasta` | FASTA | 70 viral genome sequences |
| `data/maize_pathogens_fungi.fasta` | FASTA | 191 fungal and oomycete ITS sequences |
| `data/all_taxids.json` | JSON | Complete species catalog with NCBI TaxIDs |
| `data/taxonomy.json` | JSON | Full NCBI taxonomic lineage per species |

### BLAST Databases

| File | Description |
|---|---|
| `blast_db/maize_pathogens_all.*` | Combined BLAST nucleotide database (324 sequences) |
| `blast_db/maize_pathogens_bacteria.*` | Bacteria-only BLAST database (63 sequences) |
| `blast_db/maize_pathogens_viruses.*` | Viruses-only BLAST database (70 sequences) |
| `blast_db/maize_pathogens_fungi.*` | Fungi+Oomycetes BLAST database (191 sequences) |

BLAST databases were built with NCBI BLAST+ 2.17.0 using `makeblastdb -dbtype nucl`.

### Taxonomy Files

| File | Format | Description |
|---|---|---|
| `taxonomy/maize_pathogens_taxonomy_sintax.fasta` | SINTAX | QIIME2/VSEARCH-compatible taxonomy file with full lineage |
| `taxonomy/taxids_bacteria.txt` | TXT | List of bacterial TaxIDs |
| `taxonomy/taxids_viruses.txt` | TXT | List of viral TaxIDs |
| `taxonomy/taxids_fungi.txt` | TXT | List of fungal TaxIDs |

### Source Data (Literature Mining)

| File | Description |
|---|---|
| `source_data/玉米细菌性病原菌数据库.xlsx` | Bacterial pathogens: taxonomy, diseases, literature references |
| `source_data/玉米病毒性病原菌数据库.xlsx` | Viral pathogens: taxonomy, diseases, literature references |
| `source_data/玉米真菌性病原菌数据库.xlsx` | Fungal and oomycete pathogens: taxonomy, diseases, literature references |

### Validation Results

| File | Description |
|---|---|
| `validation/validation_summary.json` | Multi-method validation results |
| `validation/ncbi_nt_comparison_v2.json` | MaizePathogenDB vs NCBI-nt benchmark (200 queries) |
| `validation/external/external_validation_final.json` | External validation on 2025-2026 NCBI sequences |

### Web Platform

| File | Description |
|---|---|
| `web_platform/index.html` | Single-file web application for species and sequence search |

---

## Usage

### BLAST Search

```bash
# Nucleotide BLAST against the combined database
blastn -query your_sequences.fasta \
  -db blast_db/maize_pathogens_all \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -max_target_seqs 5 \
  -num_threads 4 \
  -out results.txt

# BLAST against a specific taxonomic group
blastn -query your_sequences.fasta \
  -db blast_db/maize_pathogens_bacteria \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -max_target_seqs 5 \
  -out bacteria_results.txt
```

### QIIME2 Classification

```bash
# Import the SINTAX taxonomy file
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path taxonomy/maize_pathogens_taxonomy_sintax.fasta \
  --output-path maize_pathogens_taxonomy.qza

# Train a Naive Bayes classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads reference_seqs.qza \
  --i-reference-taxonomy maize_pathogens_taxonomy.qza \
  --o-classifier maize_pathogens_classifier.qza

# Classify your sequences
qiime feature-classifier classify-sklearn \
  --i-classifier maize_pathogens_classifier.qza \
  --i-reads your_sequences.qza \
  --o-classification classification.qza
```

### VSEARCH / USEARCH SINTAX

```bash
vsearch --sintax your_sequences.fasta \
  --db taxonomy/maize_pathogens_taxonomy_sintax.fasta \
  --sintax_cutoff 0.8 \
  --tabbedout sintax_results.txt
```

### Web Platform

A web-based search interface is available for species lookup by name or taxonomy, and for sequence similarity search via 25-mer fingerprint matching.

**Local usage:**
```bash
cd web_platform
python3 -m http.server 8766
# Open http://localhost:8766 in a browser
```

**Online deployment:** The web platform is hosted at [https://USERNAME.github.io/MaizePathogenDB](https://USERNAME.github.io/MaizePathogenDB).

---

## FASTA Header Format

```
>TaxID|Species|Category|Disease_CN GenBank_Accession Organism_Description
```

**Example:**
```
>66269|Pantoea stewartii subsp. stewartii|bacteria|斯图尔特细菌性枯萎病/叶枯病 LC928073.1 Pantoea stewartii 16S rRNA
```

Fields:
1. **TaxID**: NCBI Taxonomy identifier
2. **Species**: Full scientific name (with subspecies or synonym if applicable)
3. **Category**: One of `bacteria`, `viruses`, or `fungi`
4. **Disease_CN**: Chinese common name of the associated disease
5. **Accession**: NCBI GenBank accession number
6. **Description**: Original sequence description from NCBI

---

## Taxonomy Reference

All entries were manually verified against the [NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy). Bacterial nomenclature follows the International Code of Nomenclature of Prokaryotes (ICNP). Viral taxonomy follows the International Committee on Taxonomy of Viruses (ICTV). Oomycetes are taxonomically placed within Chromista (Stramenopiles) and annotated separately from true Fungi.

The complete NCBI lineage (kingdom → phylum → class → order → family → genus → species) is provided in `taxonomy.json` for every entry.

---

## Validation Summary

| Method | Bacteria | Viruses | Fungi/Oomycetes | Overall |
|---|---|---|---|---|
| Internal self-hit | 100% (63/63) | 100% (70/70) | 99.5% (190/191) | 99.8% (323/324) |
| External (2025-2026) | 100% (10/10) | 100% (5/5) | 93.4% (57/61) | 94.7% (72/76) |
| SILVA/UNITE cross-val | 95.0% (19/20) | N/A | 90.6% (58/64) | 91.7% (77/84) |
| vs NCBI-nt (n=200) | 100% (23/23) | 100% (32/32) | 98.6% (143/145) | 99.0% (198/200) |
| Genus-level (external) | 100% | 100% | 95.1% | 96.1% |
| Genus-level (SILVA/UNITE) | 95.0% | N/A | 96.9% | — |

---

## Data Sources

Pathogen species were compiled from six independent literature sources:
1. NCBI PubMed (MeSH term + keyword search)
2. Web of Science Core Collection (2000-2026 maize disease reviews)
3. Google Scholar (grey literature and regional disease reports)
4. Chinese academic databases (CNKI, Wanfang, Baidu Scholar)
5. University research group publications
6. CABI Crop Protection Compendium and APS Compendium of Corn Diseases

---


## Code Availability

All scripts used for database construction and validation are available in the `scripts/` directory.

| Script | Function |
|---|---|
| `download_sequences_v2.py` | NCBI marker gene sequence download via Entrez API |
| `build_database.py` | BLAST database construction and SINTAX formatting |
| `ncbi_nt_comparison_v2.py` | Stratified benchmark against NCBI-nt database (200 queries) |
| `external_validation.py` | External validation using 2025–2026 NCBI sequences |
| `refseq_validation.py` | RefSeq cross-validation |

**Dependencies:** Python 3.8+ · Biopython · openpyxl · requests · matplotlib · numpy · BLAST+ 2.17.0


---

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).

---

## Links

- **GitHub repository**: [https://github.com/guifeng-gao/MaizePathogenDB](https://github.com/guifeng-gao/MaizePathogenDB)
- **Web platform**: [https://guifeng-gao.github.io/MaizePathogenDB/web/](https://guifeng-gao.github.io/MaizePathogenDB/web/)

---

## Contact

For questions, bug reports, or suggestions, please open an issue on the [GitHub repository](https://github.com/guifeng-gao/MaizePathogenDB) or contact the corresponding author at gfgao@issas.ac.cn.
