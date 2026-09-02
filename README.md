# Maize Pathogen Database (MPDB)

MaizePathogenDB (MPDB) is a manually curated, multi-kingdom reference database
for molecular identification of maize-associated pathogens. It integrates
bacteria, viruses, fungi, and oomycetes in one resource.

## Release

- Release date: 2026-08-25
- Frozen release files use version-neutral names.
- Catalog species: 225
- Species with reference sequences: 201
- Reference sequences: 6,133 (402 bacterial 16S rRNA, 5,301 fungal/oomycete
  ITS, 430 complete virus genomes)
- Candidate evidence citations: 260 (207 Correct; 40 Incorrect; 13 Not
  found; all 260 retained for audit)

## Repository structure

```text
sequences/   Multi-FASTA reference files and SINTAX taxonomy
blast_db/    Prebuilt BLAST+ databases
taxonomy/    NCBI Taxonomy-verified taxonomy JSON
web/         Standalone web search platform
data/        Species catalog, sequence manifest, QC report, evidence audit
validation/  Validation protocol and results
curation/    Literature audit, PRISMA flow, expert review record
code/        Reproducible scripts, configuration, and environment files
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

- Independent external positives: 675 (fungi 537, oomycetes 76, bacteria 32,
  viruses 30); negatives: 500; cross-database queries: 565.
- Species-level retrieval: 78.8%; genus-level retrieval: 97.0%.
- Species-level classification (pident >= 99, qcovs >= 90): 60.9%;
  genus-level (pident >= 95, qcovs >= 70): 92.4%; specificity 93.8%.
- Cross-database consistency (SILVA/UNITE/RefSeq + NCBI): 75.8% species /
  95.0% genus overall.
- NCBI ITS_eukaryote and ITS_RefSeq comparisons and a UNITE QIIME2 comparison
  are included in `validation/results/`.

The official recommended thresholds are species-level `pident >= 99,
qcovs >= 90` and genus-level `pident >= 95, qcovs >= 70`. The database is
intended as a focused confirmation and screening resource, not as a universal
substitute for general-purpose databases.

## Literature and evidence audit

The catalog was compiled from a documented multi-source search strategy.
The PubMed arm was reproduced on 2026-08-25. The full 18,014-record Web of
Science export and CNKI ENW exports are kept in the source working repository;
the Figshare package contains the analysis summaries, query records, and
review workbook instead of the raw exports. Google Scholar, Wanfang, and
Baidu Scholar hit counts were reproduced on 2026-08-26 as approximate
identification counts. `Maize Pathogen.xlsx` is the raw expert-curated source
workbook; the expert review record is in `curation/EXPERT_REVIEW.md`
(bacteria Xiaolong Shao, fungi and oomycetes Lingmin Meng, viruses Zihao
Xia).

## License

This work is licensed under CC BY 4.0. Full legal text:
https://creativecommons.org/licenses/by/4.0/legalcode

## Citation

Please cite this release when using MPDB. A DOI will be assigned on upload.
