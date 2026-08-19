# NCBI-nt 150-Query Comparison Summary

Date: 2026-08-18

Method: stratified random sampling of 150 marker sequences (seed 42, one sequence per species) from MaizePathogenDB, compared against NCBI-nt via web BLAST (top-1 hit, genus-level species check).

| Category | n | MaizePathogenDB | NCBI-nt |
|---|---:|---:|---:|
| Bacteria (16S) | 14 | 14/14 (100.0%) | 14/14 (100.0%) |
| Viruses (Genome) | 26 | 26/26 (100.0%) | 23/26 (88.5%) |
| Fungi (ITS) | 110 | 110/110 (100.0%) | 100/110 (90.9%) |
| **Overall** | **150** | **150/150 (100.0%)** | **137/150 (91.3%)** |

NCBI-nt errors (13): High Plains wheat mosaic virus (HPWMoV), Maize stripe virus (MSpV), Maize-associated totivirus, Globisporangium irregulare/Pythium irregulare, Helminthosporium hawaiiensis/Curvularia hawaiiensis, Xepiculopsis graminea/Myrothecium gramineum, Epicoccum nigrum, Rhizoctonia cerealis/Ceratobasidium cereale, Curvularia lunata (Cochliobolus lunatus), Helminthosporium pedicellatum/Exserohilum pedicellatum, Fusarium sacchari, Trichoderma roseum/Trichothecium roseum, Neocosmospora solani/Fusarium solani.

Notes:
- 6 queries initially timed out at NCBI and were successfully retried with a longer wait; no timeout entries remain.
- 4 cached results were recovered during the restart to avoid counting them as failures.
- Most NCBI-nt errors are synonyms or nomenclature differences at 100% identity (e.g. Tenuivirus zeae vs Maize stripe virus, Myrothecium gramineum vs Xepiculopsis graminea), rather than true mismatches.

Files:
- `ncbi_nt_comparison_v2.json` (full per-query results)
- `Fig_NCBI_nt_Comparison_Final_v2.pdf` (summary figure)
- `Fig_NCBI_nt_Comparison.pdf` / `Fig_NCBI_nt_Comparison_Final.pdf` (copies of the summary figure)
