# External Validation Summary: 2025-2026 NCBI New Sequences

Run date: 2026-08-19

Strict criterion: top-1 BLAST hit has the same NCBI TaxID as the query species.

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 2 | 2 | 100.0% |
| viruses | 1 | 1 | 100.0% |
| fungi | 91 | 69 | 75.8% |
| oomycetes | 9 | 4 | 44.4% |
| **Overall** | **103** | **76** | **73.8%** |

Genus-level criterion (top-1 hit shares query genus):

| Category | n | Genus correct | Genus accuracy |
|---|---:|---:|---:|
| bacteria | 2 | 2 | 100.0% |
| viruses | 1 | 1 | 100.0% |
| fungi | 91 | 82 | 90.1% |
| oomycetes | 9 | 6 | 66.7% |

One sequence per taxon (strict):

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 2 | 2 | 100.0% |
| viruses | 1 | 1 | 100.0% |
| fungi | 60 | 44 | 73.3% |
| oomycetes | 5 | 3 | 60.0% |

Note: this run uses the current 225-species catalog and only sequences not already in MaizePathogenDB. The previous 94.1% external result used a different catalog/test set and a more lenient correctness judgment, so the two numbers are not directly comparable.

Errors (27):
- Fusarium graminearum (Gibberella zeae) (fungi, pident=99.3%, top=56641|Fusarium)
- Fusarium fujikuroi（Gibberella fujikuroi ） (fungi, pident=99.4%, top=948311|Fusarium)
- Bipolaris maydis (Cochliobolus heterostrophus) (fungi, pident=100.0%, top=145392|Curvularia)
- Cladosporium herbarum (fungi, pident=100.0%, top=456999|Rhizoctonia)
- Pythium graminicola (oomycetes, pident=100.0%, top=140044|Pythium)
- Pythium graminicola (oomycetes, pident=97.8%, top=140044|Pythium)
- Peronosclerospora maydis (oomycetes, pident=88.2%, top=126844|Phytopythium)
- Peronosclerospora maydis (oomycetes, pident=88.0%, top=126844|Phytopythium)
- Peronosclerospora maydis (oomycetes, pident=89.2%, top=126844|Phytopythium)
- Curvularia lunata (Cochliobolus lunatus) (fungi, pident=95.5%, top=318706|Curvularia)
- Alternaria alternata (fungi, pident=99.5%, top=1187904|Alternaria)
- Alternaria alternata (fungi, pident=99.2%, top=119927|Alternaria)
- Rhizopus stolonifer (fungi, pident=98.5%, top=97093|Trichoderma)
- Curvularia geniculata(Cochliobolus geniculatus) (fungi, pident=100.0%, top=215132|Curvularia)
- Penicillium chrysogenum (fungi, pident=97.7%, top=60171|Penicillium)
- Cladosporium cladosporioides (fungi, pident=98.6%, top=29918|Cladosporium)
- Alternaria tenuissima (fungi, pident=99.5%, top=1187904|Alternaria)
- Helminthosporium pedicellatum/Exserohilum pedicellatum (fungi, pident=98.6%, top=40559|Botrytis)
- Curvularia australiensis (fungi, pident=98.4%, top=1230527|Helminthosporium)
- Curvularia australiensis (fungi, pident=99.8%, top=145392|Curvularia)
- Aspergillus glaucus (fungi, pident=82.5%, top=69781|Penicillium)
- Fusarium chlamydosporum (fungi, pident=99.8%, top=61235|Fusarium)
- Fusarium boothii (fungi, pident=99.8%, top=56641|Fusarium)
- Penicillium aurantiogriseum/Penicillium cyclopium/Penicillium viridicatum (fungi, pident=98.7%, top=60171|Penicillium)
- Neoascochyta exitialis/ Didymella exitialis (fungi, pident=84.7%, top=40559|Botrytis)
- Neoascochyta exitialis/ Didymella exitialis (fungi, pident=0.0%, top=)
- Neoascochyta exitialis/ Didymella exitialis (fungi, pident=0.0%, top=)
