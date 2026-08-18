# MaizePathogenDB

**Curated Reference Database for Maize (*Zea mays* L.) Pathogen Identification**

A manually curated, multi-kingdom reference database of 225 maize-associated pathogen species with verified NCBI Taxonomy classifications and 573 marker gene sequences. The database is designed for taxonomic classification of metagenomic and amplicon sequencing data targeting the 16S rRNA gene (bacteria), ITS region (fungi and oomycetes), and complete or partial genomes (viruses).

---

## Dataset Overview

| Attribute | Count |
|---|---|
| Total pathogen species | 225 |
| Marker gene sequences | 573 |
| Species with sequences | 198/225 |
| Sequence coverage | 88.0% (198/225) |
| BLAST nucleotide databases | 4 (all, bacteria, viruses, fungi) |
| Taxonomic groups covered | 4 kingdoms (Bacteria, Fungi, Oomycota, Viruses) |
| Marker genes | 16S rRNA, ITS region, viral complete genomes |

### Breakdown by Taxonomic Group

| Group | Species | Sequences | Marker Gene | Coverage |
|---|---|---|---|---|
| Bacteria | 21 | 51 | 16S rRNA | 90.5% (19/21) |
| Viruses | 36 | 96 | Complete genome | 94.4% (34/36) |
| True Fungi | 136 | 426 | ITS region | — |
| Oomycetes | 32 | — | ITS region | — |

Note: True Fungi and Oomycetes are stored together in the fungi FASTA and BLAST database; 145 of 168 fungi/oomycete taxa have reference sequences.

### Species without Sequences (27)

Fourteen taxa lack NCBI Taxonomy entries (kept in the catalog with TaxID `NOT_VERIFIED`) and an additional thirteen taxa lack publicly available marker gene sequences in NCBI GenBank at the time of database construction.

---

## File Inventory

### Sequence Data

| File | Format | Contents |
|---|---|---|
| `maize_pathogens_all.fasta` | FASTA | All 573 marker gene sequences (merged) |
| `maize_pathogens_bacteria.fasta` | FASTA | 51 bacterial 16S rRNA sequences |
| `maize_pathogens_viruses.fasta` | FASTA | 96 viral genome sequences |
| `maize_pathogens_fungi.fasta` | FASTA | 426 fungal and oomycete ITS sequences |
| `all_taxids.json` | JSON | Complete species catalog with NCBI TaxIDs |
| `taxonomy.json` | JSON | Full NCBI taxonomic lineage per species |

### BLAST Databases

| File | Description |
|---|---|
| `blast_db/maize_pathogens_all.*` | Combined BLAST nucleotide database (573 sequences) |
| `blast_db/maize_pathogens_bacteria.*` | Bacteria-only BLAST database (51 sequences) |
| `blast_db/maize_pathogens_viruses.*` | Viruses-only BLAST database (96 sequences) |
| `blast_db/maize_pathogens_fungi.*` | Fungi+Oomycetes BLAST database (426 sequences) |

BLAST databases were built with NCBI BLAST+ 2.17.0 using `makeblastdb -dbtype nucl`.

### Taxonomy Files

| File | Format | Description |
|---|---|---|
| `maize_pathogens_taxonomy_sintax.fasta` | SINTAX | QIIME2/VSEARCH-compatible taxonomy file with full lineage |
| `all_taxids.json` | JSON | Species catalog with NCBI TaxIDs and disease annotations |

### Source Data

| File | Description |
|---|---|
| `Maize Pathogen.xlsx` | Curated pathogen catalog with three sheets (bacteria, virus, fungi and Oomycota), taxonomy, and literature references |

### Web Platform

| File | Description |
|---|---|
| `web/index.html` | Single-file web application for species search, sequence search, pathogen annotation, and new-pathogen submission |

---

## Usage

### BLAST Search

```bash
blastn -query your_sequences.fasta \
  -db blast_db/maize_pathogens_all \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -max_target_seqs 5 \
  -num_threads 4 \
  -out results.txt
```

### QIIME2 / VSEARCH Integration

```bash
vsearch --sintax rep_seqs.fasta \
  --db maize_pathogens_taxonomy_sintax.fasta \
  --sintax_cutoff 0.8 \
  --output taxonomy.txt
```

### Web Platform

A standalone web application is available at [https://guifeng-gao.github.io/MaizePathogenDB/web/](https://guifeng-gao.github.io/MaizePathogenDB/web/). It provides:

- Species Search: real-time keyword matching in English and Chinese
- Sequence Search: 25-mer fingerprint matching against 573 reference sequences
- Pathogen Annotation: upload a species CSV/TSV and receive an annotated file marking matched maize pathogens
- Submit Pathogen: report new species for database review

All computation runs client-side and requires no installation.

---

## FASTA Header Format

```
>TaxID|Species|Category|Disease GenBank_Accession Original_Description
```

Fields: TaxID (NCBI Taxonomy), Species (scientific name with synonym if applicable), Category (`bacteria`, `viruses`, or `fungi`), Disease, Accession, and original description.

---

## Taxonomy Reference

All entries were manually verified against the [NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy). Bacterial nomenclature follows ICNP. Viral taxonomy follows ICTV. Oomycetes are placed within Chromista (Stramenopiles) and annotated separately from true Fungi.

---

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
