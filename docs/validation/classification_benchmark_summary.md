# Maize-Pathogen Classification Benchmark: Our DB vs NCBI-nt

Date: 2026-08-20

Positive set: 573 sequences from MaizePathogenDB (51 bacteria, 96 viruses, 361 true fungi, 65 oomycetes).
Negative set: 500 sequences downloaded from NCBI (100 bacteria, 100 viruses, 250 true fungi, 50 oomycetes), selected from non-maize-pathogen taxa with no sequence overlap with the 573 positives.

Classification rules:
- MaizePathogenDB: top-1 BLAST hit with pident >= 99.5% and qcovs >= 99%, and no background hit within 0.5% identity (non-maize-pathogen background library).
- NCBI-nt: top-1 web BLAST hit with pident >= 99.5% and hit TaxID belonging to the maize-pathogen catalog.

## MaizePathogenDB

| Category | n_pos | n_neg | Sensitivity | Specificity | F1 | Balanced accuracy |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 51 | 100 | 100.0% | 100.0% | 100.0% | 100.0% |
| Viruses | 96 | 100 | 100.0% | 100.0% | 100.0% | 100.0% |
| True Fungi | 361 | 250 | 99.2% | 100.0% | 99.6% | 99.6% |
| Oomycetes | 65 | 50 | 100.0% | 100.0% | 100.0% | 100.0% |
| **Overall** | **573** | **500** | **99.5%** | **100.0%** | **99.7%** | **99.7%** |

## NCBI-nt

| Category | n_pos | n_neg | Sensitivity | Specificity | F1 | Balanced accuracy |
|---|---:|---:|---:|---:|---:|---:|
| Bacteria | 51 | 100 | 78.4% | 100.0% | 87.9% | 89.2% |
| Viruses | 96 | 100 | 92.7% | 100.0% | 96.2% | 96.4% |
| True Fungi | 361 | 250 | 68.1% | 99.2% | 80.8% | 83.7% |
| Oomycetes | 65 | 50 | 84.6% | 100.0% | 91.7% | 92.3% |
| **Overall** | **573** | **500** | **75.0%** | **99.6%** | **85.6%** | **87.3%** |

Notes:
- With the stricter threshold and background rejection, MaizePathogenDB reached 100% specificity on this negative set and higher F1/balanced accuracy than NCBI-nt.
- Because the 573 positives come from MaizePathogenDB itself, its sensitivity is expected to be high.
- Some NCBI-nt negative queries timed out or returned no TaxID and were treated as negative. All were negatives, so NCBI-nt specificity may be slightly optimistic.

Files:
- `classification_benchmark/benchmark_metrics.json`
- `classification_benchmark/benchmark_ourdb_predictions.json`
- `classification_benchmark/benchmark_ncbi_predictions.json`
- `classification_benchmark/negative_queries.fasta`
- `classification_benchmark/negative_queries_meta.json`
