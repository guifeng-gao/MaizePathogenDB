# Performance (Usage Notes, informational)

**Purpose**: resource characteristic only. No pathogen detection or ecological analysis.  
**Query set**: fixed ASVs from `260samples_fungi/rep-seqs.fasta` (19,276), deterministic subsets with seed 42.  
**Environment**: BLAST+ 2.17.0, 4 threads, macOS; wall time and maximum resident set size.

| n_queries | MPDB release (s) | MPDB max RSS (MB) | NCBI ITS_eukaryote (s) | NCBI max RSS (MB) |
|---:|---:|---:|---:|---:|
| 100 | 1.9 | 43.9 | 2.9 | 57.4 |
| 500 | 2.6 | 114.8 | 11.0 | 60.0 |
| 1,000 | 4.2 | 115.2 | 21.0 | 60.8 |
| 2,000 | 7.7 | 116.3 | 41.6 | 62.0 |
| 5,000 | 18.1 | 116.2 | 95.6 | 63.0 |
| 10,000 | 34.8 | 117.8 | 188.4 | 63.6 |
| 19,276 | 78.0 | 336.6 | 499.0 | 183.2 |

Per-query: MPDB ~4 ms; NCBI ITS_eukaryote ~26 ms at full scale.

Full rows: `performance.json`.
