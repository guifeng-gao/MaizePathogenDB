# MPDB marker-cleaned reanalysis (2026-09-16)

This release is a marker-cleaned rebuild derived from the 6,133-sequence MPDB
release.

## Data changes

- Removed 10 bacterial false positives whose titles described 16S rRNA
  methylase or methyltransferase genes, not 16S rRNA genes.
- Extracted ITS1-5.8S-ITS2 for fungal and oomycete records using ITSx 1.1.3.
- Extracted 16S rRNA regions for bacterial records using Barrnap 0.9, with
  marker-only partial records retained when they were already trimmed.
- Excluded 9 fungal records for which no reliable pure ITS interval could be
  obtained. No species lost all reference sequences.

Final marker database: 6,114 sequences.

| Category | Sequences |
|---|---:|
| Bacteria | 392 |
| Fungi | 4,500 |
| Oomycetes | 792 |
| Viruses | 430 |

## Marker length and GC content

| Category | Length min | Length median | Length max | GC min | GC median | GC max |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 120 | 1,022 | 1,533 | 45.0% | 54.8% | 58.2% |
| Fungi | 56 | 482 | 992 | 18.8% | 49.6% | 72.8% |
| Oomycetes | 89 | 781 | 1,574 | 30.0% | 45.6% | 65.4% |
| Viruses | 720 | 5,676.5 | 13,800 | 30.4% | 44.9% | 61.7% |

Viral values are unchanged because complete genomes or genome segments were
retained rather than trimmed.

## Primer coverage

Primer coverage was calculated on the cleaned source records after restricting
them to the marker-clean sequence IDs. This is necessary because several ITS
primers bind outside the ITS1-5.8S-ITS2 interval itself.

| Primer pair | Covered | n | Coverage |
|---|---:|---:|---:|
| 16S V3-V4 (338F/806R) | 235 | 392 | 59.9% |
| 16S V4 (515F/806R) | 264 | 392 | 67.3% |
| 16S V3-V8 (338F/1392R) | 133 | 392 | 33.9% |
| 16S full length (27F/1492R) | 4 | 392 | 1.0% |
| 16S V4-V5 (515F/926R) | 244 | 392 | 62.2% |
| ITS1F/ITS4 | 167 | 5,292 | 3.2% |
| ITS5/ITS4 | 198 | 5,292 | 3.7% |
| ITS1/ITS4 | 381 | 5,292 | 7.2% |
| ITS86F/ITS4 | 851 | 5,292 | 16.1% |
| fITS7/ITS4 | 878 | 5,292 | 16.6% |
| ITS1F/ITS2 | 295 | 5,292 | 5.6% |
| ITS3/ITS4 | 835 | 5,292 | 15.8% |
| ITS9mun/ITS4ngs | 56 | 5,292 | 1.1% |

## Internal completeness

Top-1 self-hit rate against the same-category marker database:

| Category | Correct | n | Rate |
|---|---:|---:|---:|
| Bacteria | 392 | 392 | 100.0% |
| Viruses | 430 | 430 | 100.0% |
| Fungi | 4,023 | 4,500 | 89.4% |
| Oomycetes | 779 | 792 | 98.4% |

The reduction in fungi is partly caused by identical ITS sequences being
shared across taxa, which creates BLAST ties at species level.

## Validation

The validation queries were trimmed to the same marker space as the database.
All 675 positive queries were retained. The negative set retained 487 of 500,
and the cross-database set retained 560 of 565.

| Validation | Overall result |
|---|---:|
| Species retrieval, no threshold | 76.1% (514/675) |
| Genus retrieval, no threshold | 96.7% (653/675) |
| Species classification, 99/90 | 61.6% |
| Genus classification, 95/70 | 93.8% |
| Species sensitivity/specificity/F1/balanced | 61.6% / 93.8% / 74.2% / 77.7% |
| Fixed split species sensitivity/specificity/F1 | 57.9% / 94.6% / 71.5% |
| Fixed split genus sensitivity/specificity/F1 | 92.6% / 91.9% / 93.3% |
| Cross-database species retrieval | 77.7% (435/560) |
| Cross-database genus retrieval | 95.9% (537/560) |

Fungal-only sensitivity in the fixed-database comparison was 60.1% at species
level for MPDB, 17.9% for NCBI ITS_eukaryote, 16.6% for ITS_RefSeq_Fungi, and
37.8% for UNITE at confidence >= 0.7.

## Performance

On 19,276 ASV queries, MPDB completed in 75.447 s with 369.1 MB peak RSS
(3.914 ms/query). NCBI ITS_eukaryote completed in 528.216 s with 154.0 MB
(27.403 ms/query). This is an informational benchmark only.

## NCBI-nt comparison

The fixed NCBI-nt snapshot is not present locally and is approximately
4,315,573,502,030 bp. It was therefore not searched again. The previous
NCBI-nt comparison should not be combined with the corrected MPDB numbers
unless that snapshots is restored and both sides are rerun.
