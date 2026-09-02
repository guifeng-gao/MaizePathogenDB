# Classification vs UNITE (QIIME2)

**Data**: release; positive/negative fungal queries only (UNITE dynamic fungi has no oomycete classes).
**Classifier**: UNITE v10.0 dynamic fungi, release 2025-02-19, `unite_ver2025-02-19_dynamic_fungi-Q2-2026.4.qza`.
**QIIME2**: QIIME2 2026.7 (rachis 2026.7.0).
**Confidence threshold**: >= 0.7; calls below threshold are treated as no call.
**Queries**: 787 (positive 537, negative 250).

| Level | TP | FP | FN | TN | Sensitivity | Specificity | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| species | 208 | 1 | 329 | 249 | 38.7% | 99.6% | 55.8 | 69.2% |
| genus | 449 | 45 | 88 | 205 | 83.6% | 82.0% | 87.1 | 82.8% |

## Notes

- Oomycetes are excluded because the UNITE dynamic fungi classifier contains no oomycete reference classes; forcing oomycete queries through it would report them as Fungi and produce misleading metrics.
- Species matching uses NCBI taxonomy names and synonyms for each query TaxID; UNITE species names that are not present in NCBI taxonomy are treated as mismatches and should be reviewed manually.
- Full per-query calls are not included in the public package.

## Comparison on fungal queries only (same query set)

UNITE is a naive-Bayes classifier with confidence >= 0.7; BLAST databases use species 99/90 and genus 95/70. Metrics are not threshold-equivalent.

| Method | Level | Sensitivity | Specificity | F1 | Balanced |
|---|---|---:|---:|---:|---:|
| MaizePathogenDB (BLAST 99/90, 95/70) | species | 59.2% | 87.6% | 71.8 | 73.4% |
| MaizePathogenDB (BLAST 99/90, 95/70) | genus | 91.8% | 82.8% | 91.9 | 87.3% |
| NCBI ITS_eukaryote (BLAST 99/90, 95/70) | species | 14.2% | 99.2% | 24.7 | 56.7% |
| NCBI ITS_eukaryote (BLAST 99/90, 95/70) | genus | 82.1% | 84.0% | 86.6 | 83.1% |

Source for MPDB/NCBI rows: `ncbi_its_comparison.json`.
