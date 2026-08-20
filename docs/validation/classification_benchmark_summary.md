# Maize-Pathogen Classification Benchmark: Our DB vs NCBI-nt

Date: 2026-08-20

Positive set: 573 sequences from MaizePathogenDB (51 bacteria, 96 viruses, 361 true fungi, 65 oomycetes).
Negative set: 500 sequences downloaded from NCBI (100 bacteria, 100 viruses, 250 true fungi, 50 oomycetes), selected from non-maize-pathogen taxa with no sequence overlap with the 573 positives.

Classification rules:
- MaizePathogenDB: top-1 BLAST hit with pident >= 97% and qcovs >= 90%.
- NCBI-nt: top-1 web BLAST hit with pident >= 97% and hit TaxID belonging to the maize-pathogen catalog.

## MaizePathogenDB

| Category | n_pos | n_neg | Sensitivity | Specificity | F1 | Balanced accuracy |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 51 | 100 | 100.0% | 100.0% | 100.0% | 100.0% |
| Viruses | 96 | 100 | 100.0% | 100.0% | 100.0% | 100.0% |
| True Fungi | 361 | 250 | 100.0% | 89.6% | 96.5% | 94.8% |
| Oomycetes | 65 | 50 | 100.0% | 100.0% | 100.0% | 100.0% |
| **Overall** | **573** | **500** | **100.0%** | **94.8%** | **97.8%** | **97.4%** |

## NCBI-nt

| Category | n_pos | n_neg | Sensitivity | Specificity | F1 | Balanced accuracy |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 51 | 100 | 78.4% | 100.0% | 87.9% | 89.2% |
| Viruses | 96 | 100 | 96.9% | 100.0% | 98.4% | 98.4% |
| True Fungi | 361 | 250 | 75.3% | 99.2% | 85.7% | 87.3% |
| Oomycetes | 65 | 50 | 87.7% | 100.0% | 93.4% | 93.8% |
| **Overall** | **573** | **500** | **80.6%** | **99.6%** | **89.1%** | **90.1%** |

Notes:
- Because the 573 positives are from MaizePathogenDB itself, its sensitivity is expected to be 100%; this makes the sensitivity comparison inherently favorable to MaizePathogenDB.
- The more informative comparison is specificity and F1 on the independent negative set.
- 70 NCBI-nt negative queries timed out or returned no TaxID and were treated as negative; all were negatives, so NCBI-nt specificity may be slightly optimistic.

Files:
- `classification_benchmark/benchmark_metrics.json`
- `classification_benchmark/benchmark_ourdb_predictions.json`
- `classification_benchmark/benchmark_ncbi_predictions.json`
- `classification_benchmark/negative_queries.fasta`
- `classification_benchmark/negative_queries_meta.json`
