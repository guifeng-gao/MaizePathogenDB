# NCBI-nt 150-Query Comparison Summary

Date: 2026-08-20

Method: stratified random sampling of 150 marker sequences (seed 42, one sequence per species) from MaizePathogenDB, compared against NCBI-nt via web BLAST. Correctness is species-level: the top-1 NCBI hit TaxID must equal the query TaxID.

| Category | n | MaizePathogenDB | NCBI-nt |
|---|---:|---:|---:|
| Bacteria (16S) | 14 | 14/14 (100.0%) | 11/14 (78.6%) |
| Viruses (Genome) | 26 | 26/26 (100.0%) | 25/26 (96.2%) |
| True Fungi (ITS) | 94 | 94/94 (100.0%) | 65/94 (69.1%) |
| Oomycetes (ITS) | 16 | 16/16 (100.0%) | 14/16 (87.5%) |
| **Overall** | **150** | **150/150 (100.0%)** | **115/150 (76.7%)** |

NCBI-nt errors: 35 (bacteria 3, viruses 1, true fungi 29, oomycetes 2). See `ncbi_nt_comparison.json` for per-query details.

Notes:
- All validation methods now use species-level TaxID comparison for comparability.
- Many NCBI-nt species-level differences reflect ITS species-level resolution and NCBI nomenclature updates.

Files:
- `ncbi_nt_comparison.json` (full per-query results)
