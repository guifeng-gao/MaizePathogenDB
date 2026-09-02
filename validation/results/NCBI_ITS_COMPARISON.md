# Classification vs Fixed NCBI ITS Databases

**Data**: MaizePathogenDB release (6,133 sequences); independent positives and
negatives for fungi + oomycetes.  
**Thresholds**: species `pident>=99, qcovs>=90`; genus `pident>=95, qcovs>=70`.  
**Fixed NCBI databases**: `ITS_eukaryote_sequences` and `ITS_RefSeq_Fungi`,
both downloaded 2026-08-21.  
**Queries**: positives 613 (fungi 537, oomycetes 76); negatives 300
(fungi 250, oomycetes 50).

## Overall (fungi + oomycetes)

| Database | Level | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|
| MaizePathogenDB | species | 61.3% | 89.7% | 92.4% | 73.7 | 75.5% |
| MaizePathogenDB | genus | 92.7% | 85.3% | 92.8% | 92.7 | 89.0% |
| NCBI ITS_eukaryote | species | 18.6% | 99.3% | 98.3% | 31.3 | 59.0% |
| NCBI ITS_eukaryote | genus | 82.5% | 86.3% | 92.5% | 87.2 | 84.4% |
| NCBI ITS_RefSeq_Fungi | species | 12.2% | 99.3% | 97.4% | 21.7 | 55.8% |
| NCBI ITS_RefSeq_Fungi | genus | 64.8% | 84.0% | 89.2% | 75.0 | 74.4% |

## Per category

### MaizePathogenDB

| Category | Level | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|
| fungi | species | 59.2% | 87.6% | 91.1% | 71.8 | 73.4% |
| fungi | genus | 91.8% | 82.8% | 92.0% | 91.9 | 87.3% |
| oomycetes | species | 76.3% | 100.0% | 100.0% | 86.6 | 88.2% |
| oomycetes | genus | 98.7% | 98.0% | 98.7% | 98.7 | 98.3% |

### NCBI ITS_eukaryote

| Category | Level | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|
| fungi | species | 14.2% | 99.2% | 97.4% | 24.7 | 56.7% |
| fungi | genus | 82.1% | 84.0% | 91.7% | 86.6 | 83.1% |
| oomycetes | species | 50.0% | 100.0% | 100.0% | 66.7 | 75.0% |
| oomycetes | genus | 85.5% | 98.0% | 98.5% | 91.5 | 91.8% |

### NCBI ITS_RefSeq_Fungi

| Category | Level | Sensitivity | Specificity | Precision | F1 | Balanced |
|---|---:|---:|---:|---:|---:|---:|
| fungi | species | 14.0% | 99.2% | 97.4% | 24.4 | 56.6% |
| fungi | genus | 73.9% | 81.2% | 89.4% | 80.9 | 77.6% |
| oomycetes | species | 0.0% | 100.0% | - | - | 50.0% |
| oomycetes | genus | 0.0% | 98.0% | 0.0% | - | 49.0% |

## Notes

- ITS_RefSeq_Fungi contains no oomycete records, so oomycete sensitivity is 0
  by construction.
- The comparisons use the same positive/negative query set and the same thresholds for all
  databases.
- Full metrics: `ncbi_its_comparison.json`.
