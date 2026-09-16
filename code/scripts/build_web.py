#!/usr/bin/env python3
"""Rebuild MaizePathogenDB web page data and validation summary."""

import csv
import json
import os
import re

ROOT = os.environ.get(
    "MPDB_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
TAXONOMY_JSON = next(
    path for path in (
        os.path.join(ROOT, "Figshare", "taxonomy", "taxonomy.json"),
        os.path.join(ROOT, "taxonomy", "taxonomy.json"),
    ) if os.path.exists(path)
)
MANIFEST = os.path.join(ROOT, "data", "sequence_manifest.tsv")
FASTA = next(
    path for path in (
        os.path.join(ROOT, "release", "sequences", "maize_pathogens_all.fasta"),
        os.path.join(ROOT, "sequences", "maize_pathogens_all.fasta"),
    ) if os.path.exists(path)
)
WEB_FILES = [
    path for path in (
        os.path.join(ROOT, "maize_pathogen_web", "index.html"),
        os.path.join(ROOT, "Figshare", "web", "index.html"),
        os.path.join(ROOT, "web", "index.html"),
    ) if os.path.exists(path)
]
DATA_DIRS = [
    path for path in (
        os.path.join(ROOT, "maize_pathogen_web", "data"),
        os.path.join(ROOT, "Figshare", "web", "data"),
        os.path.join(ROOT, "web", "data"),
    ) if os.path.exists(os.path.dirname(path))
]

MAX_KMERS = 200
CAT_MAP = {"Bacteria": "Bacteria", "Virus": "Viruses", "Viruses": "Viruses",
           "Fungi": "Fungi", "Oomycetes": "Oomycetes"}


def load_sequences():
    groups = {}
    header = None
    seq = []
    with open(FASTA, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    parts = header.split("|")
                    groups.setdefault(parts[1], {"species": parts[2],
                                                 "category": parts[3],
                                                 "lengths": [], "seqs": []})
                    groups[parts[1]]["lengths"].append(len("".join(seq)))
                    groups[parts[1]]["seqs"].append("".join(seq))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            parts = header.split("|")
            groups.setdefault(parts[1], {"species": parts[2],
                                         "category": parts[3],
                                         "lengths": [], "seqs": []})
            groups[parts[1]]["lengths"].append(len("".join(seq)))
            groups[parts[1]]["seqs"].append("".join(seq))
    return groups


def build_tax():
    records = json.load(open(TAXONOMY_JSON, encoding="utf-8"))
    out = []
    for row in records:
        taxid = str(row.get("taxid"))
        category = CAT_MAP.get(str(row.get("category")), str(row.get("category")))
        keywords = " ".join([
            (row.get("species") or "").lower(),
            (row.get("disease_en") or "").lower(),
            (row.get("disease_cn") or ""),
            (row.get("kingdom") or "").lower(),
            (row.get("phylum") or "").lower(),
            (row.get("class") or "").lower(),
            (row.get("order") or "").lower(),
            (row.get("family") or "").lower(),
            (row.get("genus") or "").lower(),
        ])
        out.append({
            "species": row.get("species", ""),
            "disease_en": row.get("disease_en", ""),
            "disease_cn": row.get("disease_cn", ""),
            "kingdom": row.get("kingdom", ""),
            "phylum": row.get("phylum", ""),
            "class": row.get("class", ""),
            "order": row.get("order", ""),
            "family": row.get("family", ""),
            "genus": row.get("genus", ""),
            "taxid": taxid,
            "category": category,
            "keywords": keywords,
        })
    return out


def build_meta_fp(groups):
    meta = {}
    fp = {}
    for taxid, info in groups.items():
        meta[taxid] = {
            "s": info["species"],
            "c": info["category"],
            "l": max(info["lengths"]),
            "n": len(info["lengths"]),
        }
        kmers = set()
        for seq in info["seqs"]:
            seq = re.sub(r"[^ACGTNacgtn]", "", seq.upper())
            for i in range(0, len(seq) - 24, 5):
                kmers.add(seq[i:i + 25])
                if len(kmers) >= MAX_KMERS * 5:
                    break
            if len(kmers) >= MAX_KMERS * 5:
                break
        ordered = sorted(kmers)
        step = max(1, len(ordered) // MAX_KMERS)
        fp[taxid] = ordered[::step][:MAX_KMERS]
    return meta, fp


def validation_block():
    # Keep these results synchronized with the manuscript Technical Validation
    # section and docs/validation/results/SUMMARY.md.
    return """
<div style="margin-top:24px">
<h3 style="font-size:16px;margin-bottom:12px">Validation Summary</h3>
<table class="result-table">
<tr><th>Analysis</th><th>Result</th></tr>
<tr><td>Internal top-1 self-hit (n=6,114)</td><td>Bacteria 100.0%; viruses 100.0%; fungi 89.4%; oomycetes 98.4%</td></tr>
<tr><td>Independent external retrieval (n=675)</td><td>Species 76.1% (514/675); genus 96.7% (653/675)</td></tr>
<tr><td>External classification at fixed thresholds (675 positives)</td><td>Species 61.6%; genus 93.8%</td></tr>
<tr><td>Classification benchmark (675 positives; 487 negatives)</td><td>Sensitivity 61.6%; specificity 93.8%; precision 93.3%; F1 74.2%; balanced accuracy 77.7%</td></tr>
<tr><td>Fixed-threshold validation split (299 positives; 221 negatives)</td><td>Species: sensitivity 57.9%, specificity 94.6%, precision 93.5%, F1 71.5%, balanced accuracy 76.2%; genus: sensitivity 92.6%, specificity 91.9%, precision 93.9%, F1 93.3%, balanced accuracy 92.2%</td></tr>
<tr><td>Cross-database consistency (n=560)</td><td>Species 77.7% (435/560); genus 95.9% (537/560)</td></tr>
<tr><td>NCBI ITS comparison (613 fungal/oomycete positives)</td><td>Species sensitivity: MPDB 62.0% vs NCBI ITS_eukaryote 23.0% and ITS_RefSeq_Fungi 14.5%; genus sensitivity: MPDB 94.0% vs 84.5% and 67.0%</td></tr>
<tr><td>UNITE comparison (fungi only; confidence >= 0.7)</td><td>Species sensitivity: MPDB 60.1% vs UNITE 37.8%; genus sensitivity: MPDB 93.3% vs UNITE 84.2%</td></tr>
</table>
<p style="font-size:12px;color:#555;margin:8px 0 0">Species calls use pident&gt;=99 and query coverage&gt;=90; genus calls use pident&gt;=95 and query coverage&gt;=70. Validation queries were trimmed to the same marker space as the corrected database. The NCBI-nt head-to-head was not rerun because the fixed local snapshot is not available.</p>
</div>
"""


def main():
    tax = build_tax()
    groups = load_sequences()
    meta, fp = build_meta_fp(groups)
    tax_json = json.dumps(tax, ensure_ascii=False, separators=(",", ":"))
    meta_json = json.dumps(meta, ensure_ascii=False, separators=(",", ":"))
    fp_json = json.dumps(fp, ensure_ascii=False, separators=(",", ":"))

    for path in WEB_FILES:
        html = open(path, encoding="utf-8").read()
        html = re.sub(r"var TAX = \[.*?\];", "var TAX = " + tax_json + ";",
                      html, count=1, flags=re.S)
        html = re.sub(r"var META = \{.*?\};", "var META = " + meta_json + ";",
                      html, count=1, flags=re.S)
        html = re.sub(r"var FP = \{.*?\};", "var FP = " + fp_json + ";",
                      html, count=1, flags=re.S)
        html = html.replace(">573</div><div class=\"label\">Marker Gene Sequences",
                            ">6133</div><div class=\"label\">Marker Gene Sequences")
        html = html.replace(">198/225</div><div class=\"label\">Species with Sequences",
                            ">201/225</div><div class=\"label\">Species with Sequences")
        html = html.replace(">6133</div><div class=\"label\">Marker Gene Sequences",
                            ">6,114</div><div class=\"label\">Reference Sequences")
        html = html.replace(">6,133</div><div class=\"label\">Reference Sequences",
                            ">6,114</div><div class=\"label\">Reference Sequences")
        html = html.replace("Fingerprint database loaded (198 species with sequences)",
                            "Fingerprint database loaded (201 species with sequences)")
        html = html.replace(
            "Extracts 25-mer fingerprints and matches against 198 species reference database.",
            "Extracts 25-mer fingerprints and performs rapid screening against a "
            "species-level fingerprint index covering the 201 species with reference sequences."
        )
        html = html.replace(
            ">4</div><div class=\"label\">BLAST Databases",
            ">5</div><div class=\"label\">BLAST Databases"
        )
        html = html.replace("<title>MaizePathogenDB</title>",
                            "<title>Maize Pathogen Database (MPDB)</title>")
        html = html.replace("<h1>🌽 MaizePathogenDB</h1>",
                            "<h1>🌽 Maize Pathogen Database (MPDB)</h1>")
        html = html.replace("not yet in MaizePathogenDB",
                            "not yet in the Maize Pathogen Database (MPDB)")
        html = html.replace("[MaizePathogenDB] New Pathogen Submission",
                            "[MPDB] New Pathogen Submission")
        html = html.replace("matches each row against MaizePathogenDB",
                            "matches each row against the Maize Pathogen Database (MPDB)")
        html = html.replace("<footer><div class=\"container\">MaizePathogenDB ·",
                            "<footer><div class=\"container\">Maize Pathogen Database (MPDB) ·")
        composition = (
            '<p style="font-size:13px;color:#555;margin:16px 0 0">'
            'Sequence composition: 392 bacterial 16S rRNA · 4,500 fungal ITS · 792 oomycete ITS · '
            '430 viral genome or genome-segment sequences</p>'
        )
        html = re.sub(
            r'<p style="font-size:13px;color:#555;margin:16px 0 0">Sequence composition:.*?</p>\s*',
            '',
            html,
            flags=re.S,
        )
        html = html.replace(
            '<div id="info-species-list"',
            composition + '\n<div id="info-species-list"',
            1,
        )

        start = html.rindex('<div style="margin-top:24px">', 0, html.index("Validation Summary"))
        end = html.index("<footer>")
        html = html[:start] + validation_block().strip() + "\n" + html[end:]
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(html)
        print("updated", path, len(html), "bytes")

    for data_dir in DATA_DIRS:
        os.makedirs(data_dir, exist_ok=True)
        json.dump(tax, open(os.path.join(data_dir, "taxonomy.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=2)
        json.dump(fp, open(os.path.join(data_dir, "fingerprints.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=2)
        json.dump(meta, open(os.path.join(data_dir, "seq_meta.json"), "w", encoding="utf-8"),
                  ensure_ascii=False, indent=2)
        print("updated data in", data_dir)


if __name__ == "__main__":
    main()
