# MaizePathogenDB release Validation Results

- Protocol: `docs/validation/PROTOCOL.md`
- Run: 2026-08-26
- Independent positive queries: 675; negative queries: 500; cross-database queries: 565

## Internal completeness

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 402 | 399 | 99.3 |
| viruses | 430 | 430 | 100.0 |
| fungi | 4509 | 4263 | 94.5 |
| oomycetes | 792 | 788 | 99.5 |

## Primer / region coverage

| Primer pair | Category | n | Covered | Coverage % |
|---|---:|---:|---:|---:|
| 16S_V3V4 | bacteria | 402 | 235 | 58.5 |
| 16S_V4 | bacteria | 402 | 264 | 65.7 |
| 16S_V3V8 | bacteria | 402 | 133 | 33.1 |
| 16S_full_27F_1492R | bacteria | 402 | 4 | 1.0 |
| 16S_V4V5 | bacteria | 402 | 244 | 60.7 |
| ITS_ITS1F_ITS4 | fungi+oomycetes | 5301 | 167 | 3.2 |
| ITS_ITS5_ITS4 | fungi+oomycetes | 5301 | 198 | 3.7 |
| ITS_ITS1_ITS4 | fungi+oomycetes | 5301 | 381 | 7.2 |
| ITS_ITS86F_ITS4 | fungi+oomycetes | 5301 | 851 | 16.1 |
| ITS_fITS7_ITS4 | fungi+oomycetes | 5301 | 878 | 16.6 |
| ITS_ITS1F_ITS2 | fungi+oomycetes | 5301 | 295 | 5.6 |
| ITS_ITS3_ITS4 | fungi+oomycetes | 5301 | 835 | 15.8 |
| ITS_ITS9mun_ITS4ngs | fungi+oomycetes | 5301 | 56 | 1.1 |

## External retrieval against MPDB

### Species-level retrieval (top-1, no threshold)

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 31 | 96.9 |
| viruses | 30 | 30 | 100.0 |
| fungi | 537 | 399 | 74.3 |
| oomycetes | 76 | 72 | 94.7 |

### Genus-level retrieval

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 32 | 100.0 |
| viruses | 30 | 30 | 100.0 |
| fungi | 537 | 517 | 96.3 |
| oomycetes | 76 | 76 | 100.0 |

### Species-level classification (pident>=99, qcovs>=90)

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 17 | 53.1 |
| viruses | 30 | 18 | 60.0 |
| fungi | 537 | 318 | 59.2 |
| oomycetes | 76 | 58 | 76.3 |

### Genus-level classification (pident>=95, qcovs>=70)

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 32 | 29 | 90.6 |
| viruses | 30 | 27 | 90.0 |
| fungi | 537 | 493 | 91.8 |
| oomycetes | 76 | 75 | 98.7 |

## Classification benchmark (species 99/90)

| Category | n_pos | n_neg | TP | FP | FN | TN | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bacteria | 32 | 100 | 17 | 0 | 15 | 100 | 53.1 | 100.0 | 100.0 | 69.4 | 76.6 |
| viruses | 30 | 100 | 18 | 0 | 12 | 100 | 60.0 | 100.0 | 100.0 | 75.0 | 80.0 |
| fungi | 537 | 250 | 318 | 31 | 219 | 219 | 59.2 | 87.6 | 91.1 | 71.8 | 73.4 |
| oomycetes | 76 | 50 | 58 | 0 | 18 | 50 | 76.3 | 100.0 | 100.0 | 86.6 | 88.2 |
| Overall | 675 | 500 | 411 | 31 | 264 | 469 | 60.9 | 93.8 | 93.0 | 73.6 | 77.3 |

## Cross-database consistency

### Species-level retrieval

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 30 | 25 | 83.3 |
| viruses | 20 | 20 | 100.0 |
| fungi | 435 | 310 | 71.3 |
| oomycetes | 80 | 73 | 91.2 |

### Genus-level retrieval

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 30 | 27 | 90.0 |
| viruses | 20 | 20 | 100.0 |
| fungi | 435 | 412 | 94.7 |
| oomycetes | 80 | 78 | 97.5 |

## Head-to-head against NCBI-nt

NCBI-nt snapshot 2026-08-23; top-1; same thresholds as MPDB. The Q2
head-to-head uses the 260 non-self queries (415 of 675 queries had their own
accession as the NCBI-nt top hit and were excluded to avoid self-match
inflation).

### Non-self Q2 comparison (n=260)

| Metric | MPDB | NCBI-nt |
|---|---:|---:|
| Species-level retrieval | 71.5% | 63.1% |
| Genus-level retrieval | 96.2% | 91.9% |
| Species-level classification (99/90) | 66.2% | 63.1% |
| Genus-level classification (95/70) | 95.4% | 91.9% |

Per category (species-level retrieval):

| Category | n | MPDB % | NCBI-nt % |
|---|---:|---:|---:|
| bacteria | 7 | 100.0 | 71.4 |
| viruses | 0 | - | - |
| fungi | 214 | 66.4 | 59.3 |
| oomycetes | 39 | 94.9 | 82.1 |

### Negative specificity (Q3, n=500)

| Category | n | Catalog false positives |
|---|---:|---:|
| bacteria | 100 | 0 |
| viruses | 100 | 0 |
| fungi | 250 | 3 |
| oomycetes | 50 | 0 |
| Overall | 500 | 3 (0.6%) |

## Pending

- NCBI-nt head-to-head: `NCBI_NT_COMPARISON.md`。
- NCBI ITS_eukaryote / ITS_RefSeq comparison: `NCBI_ITS_COMPARISON.md` / `ncbi_its_comparison.json`。
- UNITE comparison: `UNITE_COMPARISON.md` / `unite_comparison.json`（QIIME2 2026.7）。
- Performance (Usage Notes only): `PERFORMANCE.md` / `performance.json`。
