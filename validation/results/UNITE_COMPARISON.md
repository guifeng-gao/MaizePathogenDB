# Classification vs UNITE (QIIME2)

**Data**: release; positive/negative fungal queries only (UNITE dynamic fungi has no oomycete classes).
**Classifier**: UNITE v10.0 dynamic fungi, release 2025-02-19, `unite_ver2025-02-19_dynamic_fungi-Q2-2026.4.qza`.
**QIIME2**: QIIME2 2026.7 (rachis 2026.7.0).
**Confidence threshold**: >= 0.7; calls below threshold are treated as no call.
**Queries**: 783 (positive 537, negative 246).

| Level | TP | FP | FN | TN | Sensitivity | Specificity | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| species | 203 | 1 | 334 | 245 | 37.8% | 99.6% | 54.8 | 68.7% |
| genus | 452 | 46 | 85 | 200 | 84.2% | 81.3% | 87.3 | 82.7% |

## Notes

- Oomycetes are excluded because the UNITE dynamic fungi classifier contains no oomycete reference classes; forcing oomycete queries through it would report them as Fungi and produce misleading metrics.
- Species matching uses NCBI taxonomy names and synonyms for each query TaxID; UNITE species names that are not present in NCBI taxonomy are treated as mismatches and should be reviewed manually.
- Full per-query calls are not included in the public package.

## Comparison on fungal queries only (same query set)

UNITE is a naive-Bayes classifier with confidence >= 0.7; BLAST databases use species 99/90 and genus 95/70. Metrics are not threshold-equivalent.

| Method | Level | Sensitivity | Specificity | F1 | Balanced |
|---|---|---:|---:|---:|---:|
| MaizePathogenDB (BLAST 99/90, 95/70) | species | 60.1% | 87.8% | 72.6 | 74.0% |
| MaizePathogenDB (BLAST 99/90, 95/70) | genus | 93.3% | 82.5% | 92.7 | 87.9% |
| NCBI ITS_eukaryote (BLAST 99/90, 95/70) | species | 17.9% | 99.2% | 30.2 | 58.5% |
| NCBI ITS_eukaryote (BLAST 99/90, 95/70) | genus | 83.6% | 83.7% | 87.5 | 83.7% |

Source for MPDB/NCBI rows: `ncbi_its_comparison.json`.