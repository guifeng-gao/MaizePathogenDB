#!/usr/bin/env python3
"""Build standardized MaizePathogenDB FASTA, manifest and SINTAX files."""

import json
import os
import re

import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXCEL = os.path.join(ROOT, "Figshare", "Maize Pathogen.xlsx")
OLD_FASTA = os.environ.get(
    "BASE_FASTA",
    os.path.join(ROOT, "release", "sequences", "maize_pathogens_all.fasta"),
)
OLD_SINTAX = os.path.join(ROOT, "Figshare", "maize_pathogens_taxonomy_sintax.fasta")
NEW_FASTA = os.environ.get(
    "EXTRA_FASTA",
    os.path.join(ROOT, "data", "new_sequences.fasta"),
)
TAXONOMY_JSON = os.path.join(ROOT, "Figshare", "taxonomy.json")
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(ROOT, "release"))
SEQ_DIR = os.path.join(OUT_DIR, "sequences")
MANIFEST = os.environ.get("MANIFEST", os.path.join(ROOT, "data", "sequence_manifest.tsv"))
EXTRA_SOURCE = os.environ.get("EXTRA_SOURCE", "reference")

ACCESSION_RE = re.compile(r"[A-Z]{1,2}_?\d+(?:\.\d+)?")
CATEGORY_ORDER = {"bacteria": 0, "viruses": 1, "fungi": 2, "oomycetes": 3}


def clean_name(value):
    value = value.replace("\u00a0", " ").replace("\u3000", " ")
    value = value.replace("（", "(").replace("）", ")")
    return re.sub(r"\s+", " ", value).strip()


def load_catalog():
    wb = openpyxl.load_workbook(EXCEL, data_only=True)
    catalog = {}
    for ws in wb.worksheets:
        for record in ws.iter_rows(min_row=2, values_only=True):
            taxid = str(record[11]).strip() if record[11] is not None else "NOT_VERIFIED"
            species = str(record[10]).strip() if record[10] else ""
            phylum = str(record[5]).strip() if record[5] else ""
            if ws.title.lower().startswith("bacteria"):
                category = "bacteria"
            elif ws.title.lower().startswith("virus"):
                category = "viruses"
            else:
                category = "oomycetes" if "Oomycota" in phylum else "fungi"
            catalog[taxid] = {"species": clean_name(species), "category": category}
    return catalog


def load_old_sintax():
    mapping = {}
    with open(OLD_SINTAX, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                taxid = header.split(";", 1)[0]
                mapping[taxid] = header
    return mapping


def load_taxonomy_json():
    with open(TAXONOMY_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def english_tail(value):
    if not value:
        return ""
    match = re.search(r"[A-Za-z][A-Za-z0-9_ ()/.-]*$", value.strip())
    if not match:
        return ""
    text = match.group(0).strip()
    if "stramenopiles" in text.lower():
        return "Stramenopiles"
    inner = re.search(r"\(([^()]+)\)", text)
    if inner and re.search(r"[A-Za-z]", inner.group(1)):
        return inner.group(1).split("/")[0].strip()
    return text.split("/")[0].strip()


def sintax_from_taxonomy(record, taxonomy_by_taxid):
    taxon = taxonomy_by_taxid.get(record["taxid"])
    if not taxon:
        return None
    species = re.sub(r"[()（）]", "", record["species"]).replace(" ", "_")
    fields = {
        "d": english_tail(taxon.get("kingdom", "")),
        "p": english_tail(taxon.get("phylum", "")),
        "c": english_tail(taxon.get("class", "")),
        "o": english_tail(taxon.get("order", "")),
        "f": english_tail(taxon.get("family", "")),
        "g": english_tail(taxon.get("genus", "")),
        "s": species,
    }
    if not all(fields.values()):
        return None
    return ">{};tax={}".format(
        record["seqid"],
        ";".join(f"{key}:{value}" for key, value in fields.items()),
    )


def parse_records(path, source):
    records = []
    with open(path, encoding="utf-8") as fh:
        header = None
        seq = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq)))
    return records


def main():
    catalog = load_catalog()
    old_sintax = load_old_sintax()
    taxonomy = load_taxonomy_json()
    taxonomy_by_taxid = {str(r.get("taxid")): r for r in taxonomy}

    raw_records = []
    for header, seq in parse_records(OLD_FASTA, "reference"):
        raw_records.append((header, seq, "reference"))
    for header, seq in parse_records(NEW_FASTA, EXTRA_SOURCE):
        raw_records.append((header, seq, EXTRA_SOURCE))

    records = []
    for index, (header, seq, source) in enumerate(raw_records, 1):
        parts = header.split("|")
        if parts[0].startswith("MPDB"):
            taxid = parts[1].strip()
            header_species = parts[2].strip() if len(parts) > 2 else ""
            header_category = parts[3].strip() if len(parts) > 3 else ""
            rest = parts[4] if len(parts) > 4 else ""
        else:
            taxid = parts[0].strip()
            header_species = parts[1].strip() if len(parts) > 1 else ""
            header_category = parts[2].strip() if len(parts) > 2 else ""
            rest = parts[3] if len(parts) > 3 else ""
        accession_match = ACCESSION_RE.search(rest)
        accession = accession_match.group(0) if accession_match else ""
        catalog_row = catalog.get(taxid, {})
        species = catalog_row.get("species") or clean_name(
            header_species
        )
        category = catalog_row.get("category") or (
            "viruses" if header_category == "viruses" else header_category or "fungi"
        )
        seqid = f"MPDB{index:06d}"
        records.append({
            "seqid": seqid,
            "taxid": taxid,
            "species": species,
            "category": category,
            "accession": accession,
            "length": len(seq),
            "source": source,
            "original_header": header,
            "sequence": seq,
        })

    records.sort(key=lambda r: (CATEGORY_ORDER[r["category"]], r["taxid"], r["accession"]))
    for index, record in enumerate(records, 1):
        record["seqid"] = f"MPDB{index:06d}"

    os.makedirs(SEQ_DIR, exist_ok=True)
    category_files = {
        "all": os.path.join(SEQ_DIR, "maize_pathogens_all.fasta"),
        "bacteria": os.path.join(SEQ_DIR, "maize_pathogens_bacteria.fasta"),
        "viruses": os.path.join(SEQ_DIR, "maize_pathogens_viruses.fasta"),
        "fungi": os.path.join(SEQ_DIR, "maize_pathogens_fungi.fasta"),
        "oomycetes": os.path.join(SEQ_DIR, "maize_pathogens_oomycetes.fasta"),
    }
    handles = {key: open(path, "w", encoding="utf-8") for key, path in category_files.items()}
    for record in records:
        header = f">{record['seqid']}|{record['taxid']}|{record['species']}|{record['category']}|{record['accession']}"
        lines = [header]
        for i in range(0, record["length"], 80):
            lines.append(record["sequence"][i:i + 80])
        text = "\n".join(lines) + "\n"
        handles["all"].write(text)
        if record["category"] in ("bacteria", "viruses", "fungi"):
            handles[record["category"]].write(text)
        if record["category"] == "oomycetes":
            handles["oomycetes"].write(text)
            handles["fungi"].write(text)
    for handle in handles.values():
        handle.close()

    os.makedirs(os.path.dirname(MANIFEST), exist_ok=True)
    with open(MANIFEST, "w", encoding="utf-8") as fh:
        fh.write("\t".join([
            "seqid", "taxid", "species", "category", "accession",
            "length", "source", "original_header",
        ]) + "\n")
        for record in records:
            original_header = re.sub(r"[\t\r\n]+", " ", record["original_header"]).strip()
            fh.write("\t".join([
                record["seqid"], record["taxid"], record["species"],
                record["category"], record["accession"], str(record["length"]),
                record["source"], original_header,
            ]) + "\n")

    sintax_path = os.path.join(OUT_DIR, "maize_pathogens_taxonomy_sintax.fasta")
    with open(sintax_path, "w", encoding="utf-8") as fh:
        for record in records:
            generated = sintax_from_taxonomy(record, taxonomy_by_taxid)
            if generated:
                fh.write(generated + "\n")
                fh.write(record["sequence"] + "\n")
            elif record["taxid"] in old_sintax:
                old = old_sintax[record["taxid"]]
                tax_part = old.split(";", 1)[1]
                fh.write(f">{record['seqid']};{tax_part}\n")
                fh.write(record["sequence"] + "\n")
            else:
                raise RuntimeError(f"missing taxonomy for {record['seqid']} {record['taxid']}")

    from collections import Counter
    counts = Counter(r["category"] for r in records)
    species_with_seq = len({r["taxid"] for r in records})
    print(f"total sequences: {len(records)}")
    print("category counts:", dict(counts))
    print("species with sequences:", species_with_seq)
    print("manifest:", MANIFEST)
    print("sequences dir:", SEQ_DIR)
    print("sintax:", sintax_path)


if __name__ == "__main__":
    main()
