#!/usr/bin/env python3
"""Assemble the public Figshare package for MaizePathogenDB.

The package contains release sequences, prebuilt BLAST databases, taxonomy,
web files, validation results, curation records, source data, and the
reproducible code. The script also verifies that no local machine paths are
present in the assembled text files.
"""

import glob
import hashlib
import os
import re
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEST = os.environ.get(
    "FIGSHARE_DEST", os.path.join(ROOT, "..", "FigureshareMPDB")
)

IGNORE = shutil.ignore_patterns(".DS_Store", "__pycache__", ".git", "*.pyc")
LOCAL_PATH_RE = re.compile(
    r"/Users/[A-Za-z0-9_.-]+|/Volumes/[A-Za-z0-9_.-]+"
)

TEXT_EXTENSIONS = {
    ".py", ".sh", ".md", ".yaml", ".yml", ".txt", ".cfg", ".ini",
    ".json", ".html", ".js",
}

README = """\
# Maize Pathogen Database (MPDB)

MaizePathogenDB (MPDB) is a manually curated, multi-kingdom reference database
for molecular identification of maize-associated pathogens. It integrates
bacteria, viruses, fungi, and oomycetes in one resource.

## Release

- Release date: 2026-08-25
- Frozen release files use version-neutral names.
- Catalog species: 225
- Species with reference sequences: 201
- Reference sequences: 6,133 (402 bacterial 16S rRNA, 5,301 fungal/oomycete
  ITS, 430 complete virus genomes)
- Candidate evidence citations: 260 (207 Correct; 40 Incorrect; 13 Not
  found; all 260 retained for audit)

## Repository structure

```text
sequences/   Multi-FASTA reference files and SINTAX taxonomy
blast_db/    Prebuilt BLAST+ databases
taxonomy/    NCBI Taxonomy-verified taxonomy JSON
web/         Standalone web search platform
data/        Species catalog, sequence manifest, QC report, evidence audit
validation/  Validation protocol and results
curation/    Literature audit, PRISMA flow, expert review record
code/        Reproducible scripts, configuration, and environment files
```

## Quick start

Build a BLAST database:

```bash
makeblastdb -in sequences/maize_pathogens_all.fasta \\
  -dbtype nucl -out blast_db/maize_pathogens_all -title "Maize Pathogen Database (MPDB)"
```

Classify a query:

```bash
blastn -query query.fasta -db blast_db/maize_pathogens_all \\
  -outfmt "6 qseqid sseqid pident qcovs staxids" -out hits.tsv
```

Use SINTAX classification with VSEARCH:

```bash
vsearch --sintax rep_seqs.fasta \\
  --db sequences/maize_pathogens_taxonomy_sintax.fasta \\
  --sintax_cutoff 0.8 --output taxonomy.txt
```

QIIME2 integration is described in `validation/QIIME2_INSTALL.md`.

## Validation summary

- Independent external positives: 675 (fungi 537, oomycetes 76, bacteria 32,
  viruses 30); negatives: 500; cross-database queries: 565.
- Species-level retrieval: 78.8%; genus-level retrieval: 97.0%.
- Species-level classification (pident >= 99, qcovs >= 90): 60.9%;
  genus-level (pident >= 95, qcovs >= 70): 92.4%; specificity 93.8%.
- Cross-database consistency (SILVA/UNITE/RefSeq + NCBI): 75.8% species /
  95.0% genus overall.
- NCBI ITS_eukaryote and ITS_RefSeq comparisons and a UNITE QIIME2 comparison
  are included in `validation/results/`.

The official recommended thresholds are species-level `pident >= 99,
qcovs >= 90` and genus-level `pident >= 95, qcovs >= 70`. The database is
intended as a focused confirmation and screening resource, not as a universal
substitute for general-purpose databases.

## Literature and evidence audit

The catalog was compiled from a documented multi-source search strategy.
The PubMed arm was reproduced on 2026-08-25. The full 18,014-record Web of
Science export and CNKI ENW exports are kept in the source working repository;
the Figshare package contains the analysis summaries, query records, and
review workbook instead of the raw exports. Google Scholar, Wanfang, and
Baidu Scholar hit counts were reproduced on 2026-08-26 as approximate
identification counts. `Maize Pathogen.xlsx` is the raw expert-curated source
workbook; the expert review record is in `curation/EXPERT_REVIEW.md`
(bacteria Xiaolong Shao, fungi and oomycetes Lingmin Meng, viruses Zihao
Xia).

## License

This work is licensed under CC BY 4.0. Full legal text:
https://creativecommons.org/licenses/by/4.0/legalcode

## Citation

Please cite this release when using MPDB. A DOI will be assigned on upload.
"""

LICENSE = """\
Creative Commons Attribution 4.0 International (CC BY 4.0)

This work is licensed under the Creative Commons Attribution 4.0
International License. To view a copy of this license, visit
https://creativecommons.org/licenses/by/4.0/legalcode
"""


def copy_file(src, rel):
    dst = os.path.join(DEST, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)


def copy_tree(src, rel):
    dst = os.path.join(DEST, rel)
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=IGNORE)


def scan_text(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    if LOCAL_PATH_RE.search(text):
        raise SystemExit(f"local path found in {path}")


def scan_tree(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__")]
        for name in filenames:
            if name == ".DS_Store":
                continue
            path = os.path.join(dirpath, name)
            if os.path.splitext(name)[1].lower() in TEXT_EXTENSIONS:
                scan_text(path)


def prune_package():
    prune_paths = [
        "data/literature_search/attempts",
        "data/literature_search/wos_export/batches",
        "data/literature_search/wos_export/wos_ris_1-1000.ris",
        "data/literature_search/wos_export/wos_ris_full.ris",
        "data/literature_search/wos_overlap.json",
        "data/literature_search/reference_review/260_citations_review_verified.tsv",
        "data/literature_search/reference_review/browser_evidence.json",
        "data/literature_search/reference_review/strict_audit.json",
        "data/literature_search/reference_review/strict_audit_openalex.json",
        "data/literature_search/reference_review/verification_audit.json",
        "validation/results/archive",
        "validation/results/unite_work",
        "validation/results/THRESHOLD_VALIDATION_SPLIT.md",
        "validation/results/threshold_tuning.json",
        "validation/results/threshold_tuning_validation_split.json",
        "validation/results/internal_completeness.json",
        "validation/results/external_retrieval_mpdb.json",
        "validation/results/cross_database_consistency.json",
        "validation/results/ncbi_nt",
        "validation/results/ncbi_nt_comparison.json",
        "validation/results/unite_per_query.tsv",
        "validation/results/unite_taxonomy.tsv",
        "validation/PROTOCOL_CHANGELOG.md",
        "data/literature_search/WOS_CNKI_FILL_TEMPLATE.md",
        "data/literature_search/reference_review/260_citations_review.xlsx",
        "data/literature_search/reference_review/maize_pathogen_consistency.json",
    ]
    for rel in prune_paths:
        path = os.path.join(DEST, rel)
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.exists(path):
            os.remove(path)
    for path in glob.glob(
        os.path.join(DEST, "data", "literature_search", "cnki_export", "CNKI-*.enw")
    ):
        os.remove(path)
    for pattern in (
        "validation/results/ncbi_nt/batch_*.fasta",
        "validation/results/ncbi_nt/hits_*.tsv",
        "validation/results/ncbi_nt/taxidlist_*.txt",
        "validation/results/ncbi_nt/nt_blast.log",
        "data/performance_query_sets/asv_query_100.fasta",
        "data/performance_query_sets/asv_query_500.fasta",
        "data/performance_query_sets/asv_query_1000.fasta",
        "data/performance_query_sets/asv_query_2000.fasta",
        "data/performance_query_sets/asv_query_5000.fasta",
        "data/performance_query_sets/asv_query_10000.fasta",
    ):
        for path in glob.glob(os.path.join(DEST, pattern)):
            os.remove(path)
    for dirpath, dirnames, filenames in os.walk(DEST):
        dirnames[:] = [d for d in dirnames if d != ".git"]
        for name in filenames:
            if name == ".DS_Store":
                os.remove(os.path.join(dirpath, name))


def write_checksums():
    lines = []
    for dirpath, dirnames, filenames in os.walk(DEST):
        dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__")]
        for name in sorted(filenames):
            if name == "CHECKSUMS.sha256" or name == ".DS_Store":
                continue
            path = os.path.join(dirpath, name)
            rel = os.path.relpath(path, DEST)
            digest = hashlib.sha256()
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                    digest.update(chunk)
            lines.append(f"{digest.hexdigest()}  {rel}")
    checksum_path = os.path.join(DEST, "CHECKSUMS.sha256")
    with open(checksum_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    os.makedirs(DEST, exist_ok=True)
    for rel in (
        "blast_db",
        "curation",
        "data",
        "sequences",
        "taxonomy",
        "validation",
    ):
        path = os.path.join(DEST, rel)
        if os.path.exists(path):
            shutil.rmtree(path)

    # Sequences and BLAST databases
    copy_tree(os.path.join(ROOT, "release", "sequences"), "sequences")
    copy_file(
        os.path.join(ROOT, "release", "maize_pathogens_taxonomy_sintax.fasta"),
        "sequences/maize_pathogens_taxonomy_sintax.fasta",
    )
    copy_tree(os.path.join(ROOT, "release", "blast_db"), "blast_db")

    # Taxonomy and web
    copy_file(
        os.path.join(ROOT, "Figshare", "web", "data", "taxonomy.json"),
        "taxonomy/taxonomy.json",
    )
    copy_tree(os.path.join(ROOT, "Figshare", "web"), "web")

    # Data
    data_files = [
        "species_list.tsv",
        "sequence_manifest.tsv",
        "sequence_qc_report.tsv",
        "evidence_audit.tsv",
        "reference_resolution.tsv",
        "taxid_resolution.tsv",
        "literature_search_log.tsv",
        "DATA_DICTIONARY.md",
    ]
    for name in data_files:
        copy_file(os.path.join(ROOT, "data", name), f"data/{name}")
    copy_tree(os.path.join(ROOT, "data", "literature_search"), "data/literature_search")
    copy_tree(
        os.path.join(ROOT, "data", "performance_query_sets"),
        "data/performance_query_sets",
    )
    copy_file(
        os.path.join(ROOT, "Maize Pathogen.xlsx"),
        "data/source/Maize Pathogen.xlsx",
    )

    # Validation
    copy_file(
        os.path.join(ROOT, "docs", "validation", "PROTOCOL.md"),
        "validation/PROTOCOL.md",
    )
    copy_file(
        os.path.join(ROOT, "docs", "validation", "PROTOCOL_CHANGELOG.md"),
        "validation/PROTOCOL_CHANGELOG.md",
    )
    copy_file(
        os.path.join(ROOT, "docs", "validation", "README.md"),
        "validation/README.md",
    )
    copy_file(
        os.path.join(ROOT, "docs", "validation", "QIIME2_INSTALL.md"),
        "validation/QIIME2_INSTALL.md",
    )
    copy_tree(
        os.path.join(ROOT, "docs", "validation", "query_sets"),
        "validation/query_sets",
    )
    copy_tree(
        os.path.join(ROOT, "docs", "validation", "results"),
        "validation/results",
    )

    # Curation
    copy_file(
        os.path.join(ROOT, "docs", "curation", "PRISMA_LITERATURE_FLOW.md"),
        "curation/PRISMA_LITERATURE_FLOW.md",
    )
    copy_file(
        os.path.join(ROOT, "docs", "curation", "PRISMA_flow.html"),
        "curation/PRISMA_flow.html",
    )
    copy_file(
        os.path.join(ROOT, "docs", "curation", "EXPERT_REVIEW.md"),
        "curation/EXPERT_REVIEW.md",
    )
    copy_file(
        os.path.join(ROOT, "docs", "curation", "expert_review_species.tsv"),
        "curation/expert_review_species.tsv",
    )
    copy_file(
        os.path.join(
            ROOT, "docs", "curation", "PRISMA_2020_flow_diagram_filled_MPDB.docx"
        ),
        "curation/PRISMA_2020_flow_diagram_filled_MPDB.docx",
    )
    # Code
    code_src = os.path.join(ROOT, "scripts")
    code_root = os.path.join(DEST, "code")
    if os.path.exists(code_root):
        shutil.rmtree(code_root)
    code_dst = os.path.join(code_root, "scripts")
    os.makedirs(code_dst, exist_ok=True)
    kept_scripts = {
        "apply_taxonomy_resolution.py",
        "build_figshare_package.py",
        "build_query_sets.py",
        "build_species_list.py",
        "build_web.py",
        "fixed_threshold_validation_split.py",
        "qc_sequences.py",
        "rebuild_sintax.py",
        "resolve_taxonomy.py",
        "run_ncbi_its.py",
        "run_ncbi_nt.py",
        "run_performance.py",
        "run_unite.py",
        "run_validation.py",
        "smoke_test.py",
        "standardize_fasta.py",
    }
    for name in sorted(os.listdir(code_src)):
        if name.endswith(".py") and name in kept_scripts:
            shutil.copy2(os.path.join(code_src, name), os.path.join(code_dst, name))
    for name in ("run_all.sh", "environment.yml", "requirements.txt"):
        copy_file(os.path.join(ROOT, name), f"code/{name}")

    # Top-level docs
    with open(os.path.join(DEST, "README.md"), "w", encoding="utf-8") as fh:
        fh.write(README)
    with open(os.path.join(DEST, "LICENSE"), "w", encoding="utf-8") as fh:
        fh.write(LICENSE)
    prune_package()
    scan_tree(DEST)
    write_checksums()
    print("package written to", DEST)


if __name__ == "__main__":
    main()
