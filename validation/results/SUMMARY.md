# MaizePathogenDB marker-clean-2026-09-16 Validation Results

- Protocol: `validation/PROTOCOL.md`
- Run: 2026-08-26
- Independent positive queries: 675; negative queries: 500; cross-database queries: 565

## Internal completeness

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 392 | 392 | 100.0 |
| viruses | 430 | 430 | 100.0 |
| fungi | 4500 | 4022 | 89.4 |
| oomycetes | 792 | 779 | 98.4 |

## Primer / region coverage

| Primer pair | Category | n | Covered | Coverage % |
|---|---:|---:|---:|---:|
| 16S_V3V4 | bacteria | 392 | 236 | 60.2 |
| 16S_V4 | bacteria | 392 | 267 | 68.1 |
| 16S_V3V8 | bacteria | 392 | 132 | 33.7 |
| 16S_full_27F_1492R | bacteria | 392 | 4 | 1.0 |
| 16S_V4V5 | bacteria | 392 | 246 | 62.8 |
| ITS_ITS1F_ITS4 | fungi+oomycetes | 5292 | 0 | 0.0 |
| ITS_ITS5_ITS4 | fungi+oomycetes | 5292 | 1 | 0.0 |
| ITS_ITS1_ITS4 | fungi+oomycetes | 5292 | 2 | 0.0 |
| ITS_ITS86F_ITS4 | fungi+oomycetes | 5292 | 4 | 0.1 |
| ITS_fITS7_ITS4 | fungi+oomycetes | 5292 | 4 | 0.1 |
| ITS_ITS1F_ITS2 | fungi+oomycetes | 5292 | 0 | 0.0 |
| ITS_ITS3_ITS4 | fungi+oomycetes | 5292 | 4 | 0.1 |
| ITS_ITS9mun_ITS4ngs | fungi+oomycetes | 5292 | 0 | 0.0 |

## External retrieval against MPDB

### Species-level retrieval (top-1, no threshold)

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 30 | 93.8 |
| viruses | 30 | 30 | 100.0 |
| fungi | 537 | 384 | 71.5 |
| oomycetes | 76 | 70 | 92.1 |

### Genus-level retrieval

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 31 | 96.9 |
| viruses | 30 | 30 | 100.0 |
| fungi | 537 | 516 | 96.1 |
| oomycetes | 76 | 76 | 100.0 |

### Species-level classification (pident>=99, qcovs>=90)

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 18 | 56.2 |
| viruses | 30 | 18 | 60.0 |
| fungi | 537 | 323 | 60.1 |
| oomycetes | 76 | 57 | 75.0 |

### Genus-level classification (pident>=95, qcovs>=70)

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 30 | 93.8 |
| viruses | 30 | 27 | 90.0 |
| fungi | 537 | 501 | 93.3 |
| oomycetes | 76 | 75 | 98.7 |

## Classification benchmark (species 99/90)

| Category | n_pos | n_neg | TP | FP | FN | TN | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bacteria | 32 | 93 | 18 | 0 | 14 | 93 | 56.2 | 100.0 | 100.0 | 72.0 | 78.1 |
| viruses | 30 | 100 | 18 | 0 | 12 | 100 | 60.0 | 100.0 | 100.0 | 75.0 | 80.0 |
| fungi | 537 | 246 | 323 | 30 | 214 | 216 | 60.1 | 87.8 | 91.5 | 72.6 | 74.0 |
| oomycetes | 76 | 48 | 57 | 0 | 19 | 48 | 75.0 | 100.0 | 100.0 | 85.7 | 87.5 |
| Overall | 675 | 487 | 416 | 30 | 259 | 457 | 61.6 | 93.8 | 93.3 | 74.2 | 77.7 |

## Cross-database consistency

### Species-level retrieval

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 26 | 23 | 88.5 |
| viruses | 20 | 20 | 100.0 |
| fungi | 434 | 318 | 73.3 |
| oomycetes | 80 | 74 | 92.5 |

### Genus-level retrieval

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 26 | 26 | 100.0 |
| viruses | 20 | 20 | 100.0 |
| fungi | 434 | 413 | 95.2 |
| oomycetes | 80 | 78 | 97.5 |

## Additional comparisons

- NCBI-nt head-to-head: `NCBI_NT_COMPARISON.md`.
- NCBI ITS_eukaryote / ITS_RefSeq comparison: `NCBI_ITS_COMPARISON.md` and `ncbi_its_comparison.json`.
- UNITE comparison: `UNITE_COMPARISON.md` and `unite_comparison.json`.
- Performance (Usage Notes only): `PERFORMANCE.md` and `performance.json`.
