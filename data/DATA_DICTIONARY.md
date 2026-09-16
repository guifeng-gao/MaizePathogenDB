# MaizePathogenDB Data Dictionary

## species_list.tsv

| Field | Type | Description |
|---|---|---|
| species | string | NCBI scientific name (canonical for most entries; pending taxa kept as-is) |
| full_name | string | Name as recorded in the curated Excel |
| synonyms | string | Parenthesised synonyms from the species field, `;` separated |
| category | string | `bacteria` / `viruses` / `fungi` / `oomycetes` |
| taxid | string | NCBI Taxonomy ID or `NOT_VERIFIED` |
| disease_en / disease_cn | string | Disease names in English / Chinese |
| evidence_level | string | `confirmed` / `opportunistic` / `secondary` (provisional, from disease annotations) |
| evidence_refs_raw | string | Raw citation block from the curated Excel |
| sequence_status | string | `present` or `missing` |
| sequence_count | integer | Number of marker sequences in the release |
| accessions | string | GenBank accessions, `;` separated |
| missing_reason | string | Reason if `sequence_status == missing` |
| notes | string | Data-quality notes (TAXID_PENDING, duplicate names, version fixes) |

## sequence_manifest.tsv

| Field | Type | Description |
|---|---|---|
| seqid | string | Stable internal id (`MPDB######`) |
| taxid | string | NCBI Taxonomy ID |
| species | string | Catalog species name |
| category | string | Taxonomic category |
| accession | string | GenBank accession |
| length | integer | Sequence length (bp) |
| source | string | `reference` (all sequences are part of the final release) |
| original_header | string | Original GenBank FASTA header or prior-version header |
| source_length | integer | Length of the original GenBank record before marker extraction |
| marker_region | string | Extracted marker interval or retained marker-only region |
| extraction_method | string | ITSx, Barrnap, BLAST remapping, marker-only retention, or viral retention |

## sequence_qc_report.tsv

| Field | Type | Description |
|---|---|---|
| seqid | string | Internal id |
| taxid / species / category / accession | string | Record identity |
| length | integer | Sequence length |
| n_count | integer | Number of `N` bases |
| ambiguous_count | integer | Number of non-ACGTN bases |
| ambiguous_bases | string | Distinct ambiguous bases |
| header_format_ok | boolean | Header matches `MPDB{id}\|taxid\|species\|category\|accession` |
| marker_consistent | boolean | Original description matches expected marker |
| notes | string | Duplicate / cross-category flags |

## marker_extraction_exclusions.tsv

Records excluded because a reliable pure ITS interval could not be extracted.
All excluded records are fungal; no catalog entry lost all reference sequences.

## removed_false_positives.tsv

The ten bacterial records removed from the original release because their
GenBank descriptions indicated 16S rRNA methylase or methyltransferase genes,
not 16S rRNA genes.

## length_gc_summary.tsv

Marker-level length and GC summaries by category.

## length_gc_per_sequence.tsv

Per-sequence marker length, GC content, extraction method, and source length.

## validation_summary.tsv

Long-form summary of corrected internal-completeness, primer-coverage,
retrieval, classification, cross-database, fixed-split, and method-comparison
results.

## realworld_260_corrected.tsv

Per-ASV UNITE calls and corrected MPDB calls for the fixed 260-sample fungal
ITS ASV set.
