# Literature search archive (MPDB PRISMA audit)

The PubMed arm was reproduced on 2026-08-25 with NCBI E-utilities from the
search keywords documented in the manuscript Methods. The WoS and CNKI export
arms are archived under `wos_export/` and `cnki_export/`.

## Files

| File | Contents |
|---|---|
| `pubmed_broad_query.tsv` | Broad PubMed query, date, exact hit count. |
| `pubmed_broad_pmids.txt` | PMIDs returned by the broad query (capped at 9,999 by JSON ESearch). |
| `pubmed_species_queries.tsv` | 225 species-targeted queries with per-species hit counts. |
| `pubmed_species_pmids.txt` | Union of PMIDs from the 225 species-targeted queries. |
| `citation_pmid_map.tsv` | DOI -> PMID mapping for the 260 candidate citations. |
| `prisma_counts.json` | Machine-readable counts used in the PRISMA flow. |
| `literature_export_analysis.json` | WoS/CNKI export de-duplication and catalog-overlap analysis. |
| `reference_review/` | Completed 260-candidate DOI/PMID review workbook (`.xlsx` + `.tsv`; 207 Correct, 40 Incorrect, 13 Not found). |
| `wos_export/` | Full 18,014-record WoS RIS export. |
| `cnki_export/` | Provided CNKI EndNote exports (1,752 raw records; 1,199 unique). |

## Caveat

The Figshare package keeps the analysis summaries and review workbook, but
omits the raw WoS RIS batches/combined file and the CNKI ENW files to keep the
package focused. Those raw exports remain in the source working repository.

WoS, Google Scholar, CNKI, Wanfang, Baidu Scholar and compendium searches are
documented in the manuscript Methods. On 2026-08-25 we attempted Google
Scholar, Baidu Scholar and CNKI; all were blocked by security verification or
timed out, and WoS requires an institutional login. A second attempt on
2026-08-26 reproduced the WoS Core Collection search (18,014 raw; 14,145 for
Article/Review 2000-2026), received the CNKI advanced-search count from an
authorized connection (1,702 raw), and exported the full 18,014 WoS records to
`wos_export/`. The provided CNKI exports contain 1,752 raw records and 1,199
unique records after de-duplication by DOI or title+year+journal; 503 of the
reported 1,702 hits are not present in the provided export files. On
2026-08-26 Google Scholar, Wanfang and Baidu Scholar were reproduced in
Microsoft Edge and their hit counts were archived in
`data/literature_search_log.tsv` (Google Scholar and Baidu Scholar report
approximate "related results" counts, not de-duplicated exports).
