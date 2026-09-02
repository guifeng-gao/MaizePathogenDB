#!/usr/bin/env python3
"""Build data/species_list.tsv from the canonical MaizePathogenDB release."""

import os
import re
import csv
from collections import defaultdict

import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXCEL = os.path.join(ROOT, "Figshare", "Maize Pathogen.xlsx")
FASTA = os.environ.get(
    "MAIZE_FASTA",
    os.path.join(ROOT, "Figshare", "maize_pathogens_all.fasta"),
)
RESOLUTION = os.path.join(ROOT, "data", "taxid_resolution.tsv")
OUT_DIR = os.path.join(ROOT, "data")
OUT_TSV = os.path.join(OUT_DIR, "species_list.tsv")

COLUMNS = [
    "species",
    "full_name",
    "synonyms",
    "category",
    "taxid",
    "disease_en",
    "disease_cn",
    "evidence_level",
    "evidence_refs_raw",
    "sequence_status",
    "sequence_count",
    "accessions",
    "missing_reason",
    "notes",
]


def clean(value):
    if value is None:
        return ""
    text = str(value).replace("（", "(").replace("）", ")").replace("\u3000", " ")
    return re.sub(r"\s+", " ", text).strip()


def synonyms_from_species(species):
    matches = re.findall(r"\(([^)]+)\)", species)
    parts = []
    for match in matches:
        for item in re.split(r"[/、，,]|\s+and\s+", match):
            item = item.strip()
            if item and item not in parts:
                parts.append(item)
    return "; ".join(parts)


def evidence_level_from_disease(disease_en, disease_cn):
    text = f"{disease_en} {disease_cn}".lower()
    if any(k in text for k in ("opportunistic", "条件致病", "机会")):
        return "opportunistic"
    if any(k in text for k in ("occasional", "secondary", "偶发", "继发")):
        return "secondary"
    return "confirmed"


def load_sequence_map():
    seq_map = defaultdict(lambda: {"count": 0, "accessions": set()})
    accession_re = re.compile(r"[A-Z]{1,2}_?\d+(?:\.\d+)?")
    current = None
    with open(FASTA, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current = line[1:]
                first = current.split("|", 1)[0].strip()
                taxid = current.split("|")[1].strip() if first.startswith("MPDB") else first
                seq_map[taxid]["count"] += 1
                rest = current.split("|", 3)[-1]
                accession = accession_re.search(rest)
                if accession:
                    seq_map[taxid]["accessions"].add(accession.group(0))
    return seq_map


def main():
    resolution = {}
    if os.path.exists(RESOLUTION):
        with open(RESOLUTION, encoding="utf-8") as fh:
            for item in csv.DictReader(fh, delimiter="\t"):
                if item["row_type"] == "not_verified":
                    resolution[("not_verified", item["original_species"])] = item
                elif item["row_type"] == "duplicate_species_name":
                    resolution[("duplicate", item["original_taxid"])] = item

    seq_map = load_sequence_map()
    wb = openpyxl.load_workbook(EXCEL, data_only=True)
    rows = []
    for ws in wb.worksheets:
        sheet_rows = list(ws.iter_rows(values_only=True))
        for record in sheet_rows[1:]:
            full_name = clean(record[1])
            disease_en = clean(record[2])
            disease_cn = clean(record[3])
            phylum = clean(record[5])
            species = clean(record[10])
            taxid = clean(record[11])
            refs = clean(record[12])
            if ws.title.lower().startswith("bacteria"):
                category = "bacteria"
            elif ws.title.lower().startswith("virus"):
                category = "viruses"
            else:
                category = "oomycetes" if "Oomycota" in phylum else "fungi"
            if not taxid:
                taxid = "NOT_VERIFIED"
            if not species:
                species = full_name

            seq_info = seq_map.get(taxid, {"count": 0, "accessions": set()})
            sequence_count = seq_info["count"]
            accessions = "; ".join(sorted(seq_info["accessions"]))
            sequence_status = "present" if sequence_count else "missing"
            if sequence_status == "missing":
                if taxid == "NOT_VERIFIED":
                    missing_reason = "TaxID NOT_VERIFIED; no marker sequence linked"
                else:
                    missing_reason = "No public marker gene sequence found in NCBI GenBank (2026-08-25)"
            else:
                missing_reason = ""

            notes = ""
            resolution_item = resolution.get(("not_verified", species))
            if taxid == "NOT_VERIFIED" and resolution_item:
                notes = "TAXID_PENDING; NCBI Taxonomy query 2026-08-25 found no exact match"
            if taxid in ("10989", "10990"):
                notes = "NCBI species node resolved; original isolate TaxID retained"
                if taxid == "10990":
                    notes += "; renamed from Fijivirus zeae to Fijivirus alporyzae"

            rows.append({
                "species": species,
                "full_name": full_name,
                "synonyms": synonyms_from_species(species),
                "category": category,
                "taxid": taxid,
                "disease_en": disease_en,
                "disease_cn": disease_cn,
                "evidence_level": evidence_level_from_disease(disease_en, disease_cn),
                "evidence_refs_raw": refs,
                "sequence_status": sequence_status,
                "sequence_count": sequence_count,
                "accessions": accessions,
                "missing_reason": missing_reason,
                "notes": notes,
            })

    from collections import Counter
    species_counts = Counter(r["species"] for r in rows)
    for row in rows:
        if species_counts[row["species"]] > 1:
            others = [r["taxid"] for r in rows
                      if r["species"] == row["species"] and r["taxid"] != row["taxid"]]
            row["notes"] = f"DUPLICATE_SPECIES_NAME; same species name also used by taxid(s): {', '.join(others)}"

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        category_order = {"bacteria": 0, "viruses": 1, "fungi": 2, "oomycetes": 3}
        for row in sorted(rows, key=lambda r: (category_order[r["category"]], r["species"])):
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")

    print(f"wrote {OUT_TSV} ({len(rows)} rows)")
    print("categories:", dict(Counter(r["category"] for r in rows)))
    print("sequences:", sum(r["sequence_count"] for r in rows))
    print("present/missing:", sum(1 for r in rows if r["sequence_status"] == "present"),
          sum(1 for r in rows if r["sequence_status"] == "missing"))
    print("NOT_VERIFIED:", sum(1 for r in rows if r["taxid"] == "NOT_VERIFIED"))
    print("unmapped taxids in fasta:", sorted(set(seq_map) - {r["taxid"] for r in rows}))


if __name__ == "__main__":
    main()
