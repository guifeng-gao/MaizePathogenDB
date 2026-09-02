# Fixed-Threshold Validation Split

**Data**: MaizePathogenDB release (6,133 sequences), Q2 positives n=675,
Q3 negatives n=500.  
**Method**: seed 42; split by TaxID, stratified by category; fixed recommended
thresholds evaluated on both halves. No threshold selection. No species
appears in both halves.

Official and recommended thresholds:

- Species: `pident>=99, qcovs>=90`
- Genus: `pident>=95, qcovs>=70`

## Species-level (99/90), validation half

| Category | n_pos | n_neg | TP | FP | FN | TN | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bacteria | 9 | 40 | 3 | 0 | 6 | 40 | 33.3% | 100.0% | 100.0% | 50.0 | 66.7% |
| viruses | 16 | 44 | 8 | 0 | 8 | 44 | 50.0% | 100.0% | 100.0% | 66.7 | 75.0% |
| fungi | 243 | 114 | 141 | 12 | 102 | 102 | 58.0% | 89.5% | 92.2% | 71.2 | 73.7% |
| oomycetes | 31 | 23 | 22 | 0 | 9 | 23 | 71.0% | 100.0% | 100.0% | 83.0 | 85.5% |
| **Overall** | 299 | 221 | 174 | 12 | 125 | 209 | **58.2%** | **94.6%** | **93.5%** | **71.8** | **76.4%** |

## Genus-level (95/70), validation half

| Category | n_pos | n_neg | TP | FP | FN | TN | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| bacteria | 9 | 40 | 8 | 0 | 1 | 40 | 88.9% | 100.0% | 100.0% | 94.1 | 94.4% |
| viruses | 16 | 44 | 14 | 0 | 2 | 44 | 87.5% | 100.0% | 100.0% | 93.3 | 93.8% |
| fungi | 243 | 114 | 221 | 18 | 22 | 96 | 90.9% | 84.2% | 92.5% | 91.7 | 87.6% |
| oomycetes | 31 | 23 | 31 | 0 | 0 | 23 | 100.0% | 100.0% | 100.0% | 100.0 | 100.0% |
| **Overall** | 299 | 221 | 274 | 18 | 25 | 203 | **91.6%** | **91.9%** | **93.8%** | **92.7** | **91.7%** |

## Notes

- Bacteria and virus validation halves are small (9-16 positives); their
  numbers are indicative only.
- Fungal genus specificity at 95/70 is 84.2%; the strict species tier
  (99/90) is the recommended default for species-level calls.
- Full per-half metrics: `fixed_threshold_validation_split.json`.
