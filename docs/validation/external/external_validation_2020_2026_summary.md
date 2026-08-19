# External Validation Summary: NCBI New Sequences since 2020/01/01

Run date: 2026-08-19

Strict criterion: top-1 BLAST hit has the same NCBI TaxID as the query species.

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 10 | 7 | 70.0% |
| viruses | 3 | 3 | 100.0% |
| fungi | 106 | 82 | 77.4% |
| oomycetes | 6 | 1 | 16.7% |
| **Overall** | **125** | **93** | **74.4%** |

Genus-level criterion (top-1 hit shares query genus):

| Category | n | Genus correct | Genus accuracy |
|---|---:|---:|---:|
| bacteria | 10 | 9 | 90.0% |
| viruses | 3 | 3 | 100.0% |
| fungi | 106 | 97 | 91.5% |
| oomycetes | 6 | 3 | 50.0% |

One sequence per taxon (strict):

| Category | n | Correct | Accuracy |
|---|---:|---:|---:|
| bacteria | 7 | 5 | 71.4% |
| viruses | 2 | 2 | 100.0% |
| fungi | 68 | 51 | 75.0% |
| oomycetes | 3 | 1 | 33.3% |

Note: this run uses the current 225-species catalog and only sequences not already in MaizePathogenDB. The previous 94.1% external result used a different catalog/test set and a more lenient correctness judgment, so the two numbers are not directly comparable.

Errors (32):
- Pantoea agglomerans (bacteria, pident=99.1%, top=553|Pantoea)
- Pantoea agglomerans (bacteria, pident=96.2%, top=553|Pantoea)
- Pectobacterium carotovorum（Erwinia carotovora） (bacteria, pident=97.3%, top=615|Serratia)
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
- Curvularia australiensis (fungi, pident=99.2%, top=145392|Curvularia)
- Curvularia pallescens (fungi, pident=96.8%, top=1203437|Curvularia)
- Aspergillus glaucus (fungi, pident=82.5%, top=69781|Penicillium)
- Fusarium chlamydosporum (fungi, pident=99.8%, top=61235|Fusarium)
- Fusarium boothii (fungi, pident=99.8%, top=56641|Fusarium)
- Penicillium aurantiogriseum/Penicillium cyclopium/Penicillium viridicatum (fungi, pident=98.7%, top=60171|Penicillium)
- Neoascochyta exitialis/ Didymella exitialis (fungi, pident=84.7%, top=40559|Botrytis)
- Neoascochyta exitialis/ Didymella exitialis (fungi, pident=0.0%, top=)
- Neoascochyta exitialis/ Didymella exitialis (fungi, pident=0.0%, top=)
