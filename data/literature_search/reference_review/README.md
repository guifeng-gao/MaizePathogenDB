# 260 DOI/PMID candidate human review

## Purpose

Every row in `data/reference_resolution.tsv` is an automatically generated
candidate citation for a maize-associated pathogen. All 260 candidates have
been reviewed; the decisions are recorded in this workbook and the verified
set is merged into `data/reference_resolution.tsv`. This folder contains the
review workbook that tracks the decision for each candidate.

## Files

- `260_citations_review.tsv` - main review workbook (tab-separated text).

## How to review

1. Open `260_citations_review.tsv` in a spreadsheet or text editor.
2. For each row, open `doi_url` (or `pmid_url` when present).
3. Check that the candidate DOI/PMID points to a publication about the same
   disease, host (maize/corn), and pathogen as the `raw_chunk`.
4. Set `decision` to one of: `Correct`, `Incorrect`, `Unclear`, or
   `Not found`.
5. Add a short `notes` explaining the result, plus `reviewer` and
   `review_date`.
6. Do not edit `data/reference_resolution.tsv` during review. After the
   workbook is complete, merge approved citations back into the resolution
   table.

## Current counts (2026-08-26)

- Candidate citations: 260
- Catalog species covered: 225
- Candidates with DOI: 260
- Candidates with PMID (via reference resolution or DOI->PMID map): 98
- PubMed species-search anchor: 33
- WoS anchor (DOI or title): 108
- CNKI anchor (DOI or title): 15
- Anchored by PubMed, WoS, or CNKI: 140

## Review outcome (2026-08-26)

- Correct: 207
- Incorrect: 40
- Not found: 13
- Reviewer: Gui-Feng Gao
- Review period: 2026-08-20 to 2026-08-30 (review completed 2026-08-26;
  individual confirmation dates were not recorded; dates in the workbook
  are tracking labels assigned within this period)

## Consistency with Maize Pathogen.xlsx

- Normalized comparison of `Maize Pathogen.xlsx` evidence references against
  this workbook: 260/260 identical, 0 missing, 0 extra.

`review_id` is unique even though some legacy `citation_id` values such as
`NOT_VERIFIED.1` are repeated across species.
