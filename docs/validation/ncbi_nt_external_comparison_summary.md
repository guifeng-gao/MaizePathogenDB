# Fair Comparison: External 2020-2026, Our DB vs NCBI-nt

Date: 2026-08-20

Method: 125 external 2020-2026 NCBI sequences not present in MaizePathogenDB were used as the same query set for both databases. Correctness is species-level: top-1 hit TaxID must equal the query TaxID.

| Category | n | MaizePathogenDB | NCBI-nt |
|---|---:|---:|---:|
| Bacteria | 10 | 7/10 (70.0%) | 10/10 (100.0%) |
| Viruses | 3 | 3/3 (100.0%) | 3/3 (100.0%) |
| True Fungi | 106 | 82/106 (77.4%) | 80/106 (75.5%) |
| Oomycetes | 6 | 1/6 (16.7%) | 5/6 (83.3%) |
| **Overall** | **125** | **93/125 (74.4%)** | **98/125 (78.4%)** |

Notes:
- Both databases use the same 125 external query sequences and the same species-level TaxID criterion.
- NCBI-nt performs better on this external test set overall, mainly driven by bacteria and oomycetes.

Files:
- `ncbi_nt_external_comparison.json` (full per-query results)
