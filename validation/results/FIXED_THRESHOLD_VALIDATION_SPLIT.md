# Fixed-threshold validation split

TaxID-stratified split, seed 42. The fixed thresholds were evaluated once on
the held-out validation half without further optimization.

## Held-out validation half

| Level | Sensitivity | Specificity | Precision | F1 | Balanced accuracy |
|---|---:|---:|---:|---:|---:|
| Species | 57.9% | 94.6% | 93.5% | 71.5 | 76.2% |
| Genus | 92.6% | 91.9% | 93.9% | 93.3 | 92.2% |

The validation half contains 299 positive and 221 negative queries.
