# Head-to-head against NCBI-nt

- Protocol: `docs/validation/PROTOCOL.md`
- NCBI-nt snapshot: 2026-08-23 (130,097,750 sequences; 4,315,573,502,030 bp)
- Independent positive queries: 675; negative queries: 500
- Method: local blastn (BLAST+ 2.17.0), top-1 hit, full NCBI-nt, same
  thresholds as MPDB (species pident>=99/qcovs>=90; genus pident>=95/qcovs>=70)
- Raw BLAST hits and per-query details are not included in the public
  package.

## Fair-comparison subset: independent positives without self-matches (n=260)

Positive queries were downloaded from NCBI and 415 of the 675 top-1 NCBI-nt hits
were the query sequence itself (61.5%), so full-Q2 retrieval is inflated by
self-matching. For the head-to-head we report the 260 queries whose top-1
NCBI-nt hit was not the query's own accession.

| Metric (n=260 non-self queries) | MPDB | NCBI-nt |
|---|---:|---:|
| Species-level retrieval | 71.5% | 63.1% |
| Genus-level retrieval | 96.2% | 91.9% |
| Species-level classification (99/90) | 66.2% | 63.1% |
| Genus-level classification (95/70) | 95.4% | 91.9% |

Per category (non-self queries, species-level retrieval):

| Category | n | MPDB % | NCBI-nt % |
|---|---:|---:|---:|
| bacteria | 7 | 100.0 | 71.4 |
| viruses | 0 | - | - |
| fungi | 214 | 66.4 | 59.3 |
| oomycetes | 39 | 94.9 | 82.1 |

## Full positive-query context (n=675, includes self-matches)

| Category | n | Species retrieval % | Genus retrieval % | Species classification % | Genus classification % |
|---|---:|---:|---:|---:|---:|
| bacteria | 32 | 93.8 | 100.0 | 93.8 | 100.0 |
| viruses | 30 | 100.0 | 100.0 | 100.0 | 100.0 |
| fungi | 537 | 83.6 | 96.1 | 83.4 | 96.1 |
| oomycetes | 76 | 90.8 | 98.7 | 90.8 | 98.7 |
| Overall | 675 | 85.6 | 96.7 | 85.5 | 96.7 |

## Negative specificity (n=500)

- False positives (top-1 hit within MPDB catalog lineages): 3 / 500 (0.6%).
- Per category: bacteria 0/100, viruses 0/100, fungi 3/250, oomycetes 0/50.
