# Maize Pathogen Database (MPDB)

MaizePathogenDB (MPDB) is a manually curated, multi-kingdom reference database
for molecular identification of maize-associated pathogens. It integrates
bacteria, viruses, fungi, and oomycetes in one resource.

## Release

- Catalog entries: 225
- Entries with reference sequences: 201
- Reference sequences: 6,114
- Bacteria: 392
- Fungi: 4,500
- Oomycetes: 792
- Viruses: 430

## Repository structure

```text
sequences/   Marker-cleaned FASTA files and SINTAX taxonomy
blast_db/    Prebuilt BLAST+ databases
taxonomy/    NCBI Taxonomy-verified taxonomy JSON
web/         Standalone web search platform
data/        Species catalog, marker manifest, QC, validation and traceability data
validation/  Validation protocol, query sets, and corrected results
curation/    Literature audit, PRISMA flow, expert review record
code/        Database-construction and validation scripts
```

## Quick start

Build a BLAST database:

```bash
makeblastdb -in sequences/maize_pathogens_all.fasta \
  -dbtype nucl -out blast_db/maize_pathogens_all -title "Maize Pathogen Database (MPDB)"
```

Classify a query:

```bash
blastn -query query.fasta -db blast_db/maize_pathogens_all \
  -outfmt "6 qseqid sseqid pident qcovs staxids" -out hits.tsv
```

Use SINTAX classification with VSEARCH:

```bash
vsearch --sintax rep_seqs.fasta \
  --db sequences/maize_pathogens_taxonomy_sintax.fasta \
  --sintax_cutoff 0.8 --output taxonomy.txt
```

QIIME2 integration is described in `validation/QIIME2_INSTALL.md`.

## Validation summary

- Independent positive queries: 675; retained negatives: 487;
  cross-database queries: 560.
- External retrieval: species 76.1% (514/675); genus 96.7% (653/675).
- Species classification: sensitivity 61.6%; specificity 93.8%; precision
  93.3%; F1 74.2; balanced accuracy 77.7%.
- Cross-database consistency: species 77.7% (435/560); genus 95.9%
  (537/560).
- Fixed-threshold validation half: species sensitivity 57.9%, specificity
  94.6%, F1 71.5; genus sensitivity 92.6%, specificity 91.9%, F1 93.3.
- NCBI-nt was not rerun because the fixed local 4.3-Tbp snapshot is not
  available in the release environment.

The official recommended thresholds are species-level `pident >= 99,
qcovs >= 90` and genus-level `pident >= 95, qcovs >= 70`. The database is
intended as a focused confirmation and screening resource, not as a universal
substitute for general-purpose databases.

## Literature and evidence audit

The catalog was compiled from a documented multi-source search strategy.
The PubMed arm was reproduced on 2026-08-25. Web of Science, CNKI, Google
Scholar, Wanfang, and Baidu Scholar searches were archived on 2026-08-26.
`data/source/Maize Pathogen.xlsx` is the curated source workbook, and the
expert review record is in `curation/EXPERT_REVIEW.md`.

## License

This work is licensed under CC BY 4.0. Full legal text:
https://creativecommons.org/licenses/by/4.0/legalcode

## Citation

Gao, G.-F. MPDB: a multi-kingdom reference database for the identification of
maize (*Zea mays* L.) associated pathogens. figshare
https://doi.org/10.6084/m9.figshare.33415123 (2026).
