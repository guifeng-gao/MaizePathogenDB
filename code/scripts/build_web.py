#!/usr/bin/env python3
"""Rebuild MaizePathogenDB web page data and validation summary."""

import csv
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAXONOMY_JSON = os.path.join(ROOT, "Figshare", "taxonomy.json")
MANIFEST = os.path.join(ROOT, "data", "sequence_manifest.tsv")
FASTA = os.path.join(ROOT, "release", "sequences", "maize_pathogens_all.fasta")
WEB_FILES = [
    os.path.join(ROOT, "maize_pathogen_web", "index.html"),
    os.path.join(ROOT, "Figshare", "web", "index.html"),
]
DATA_DIRS = [
    os.path.join(ROOT, "maize_pathogen_web", "data"),
    os.path.join(ROOT, "Figshare", "web", "data"),
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
    return """
<div style="margin-top:24px">
<h3 style="font-size:16px;margin-bottom:12px">Validation Summary</h3>
<table class="result-table">
<tr><th>Method</th><th>Species</th><th>Genus</th></tr>
<tr><td>External retrieval (n=442)</td><td>83.3% (368/442)</td><td>97.3% (430/442)</td></tr>
<tr><td>External classification</td><td>64.7% (99/90)</td><td>92.1% (95/70)</td></tr>
<tr><td>Fixed-threshold validation split (species 99/90)</td><td>Sens 61.8%, Spec 94.6%, F1 73.8</td><td>-</td></tr>
<tr><td>Fixed-threshold validation split (genus 95/70)</td><td>-</td><td>Sens 92.0%, Spec 91.9%, F1 91.8</td></tr>
<tr><td>V4 vs NCBI ITS_eukaryote (fungi+oomycetes)</td><td>MPDB 65.4% vs 18.8%</td><td>MPDB 92.6% vs 78.9%</td></tr>
<tr><td>V4 vs NCBI ITS_RefSeq_Fungi</td><td>MPDB 65.4% vs 12.2%</td><td>MPDB 92.6% vs 60.3%</td></tr>
</table>
<p style="font-size:12px;color:#555;margin:8px 0 0">Species calls use pident>=99 and query coverage>=90; genus calls use pident>=95 and query coverage>=70.</p>
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
                            ">6,133</div><div class=\"label\">Reference Sequences")
        html = html.replace("Fingerprint database loaded (198 species with sequences)",
                            "Fingerprint database loaded (201 species with sequences)")
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
        html = html.replace(
            '<div id="info-species-list"',
            '<p style="font-size:13px;color:#555;margin:16px 0 0">'
            'Sequence composition: 402 bacterial 16S rRNA · 5,301 fungal/oomycete ITS · '
            '430 viral complete genomes</p>\n<div id="info-species-list"'
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
