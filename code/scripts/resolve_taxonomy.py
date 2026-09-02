#!/usr/bin/env python3
"""Resolve TaxIDs and names for the MaizePathogenDB catalog audit.

Covers:
- 14 catalog rows with TaxID NOT_VERIFIED
- the Fijivirus zeae duplicate (TaxIDs 10989 / 10990)
- 13 catalog rows with TaxID but no marker sequence in the release
"""

import csv
import json
import os
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECIES_LIST = os.path.join(ROOT, "data", "species_list.tsv")
OUT_TSV = os.path.join(ROOT, "data", "taxid_resolution.tsv")

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EMAIL = "maize_pathogen_db@example.com"
TOOL = "maize_pathogen_db_resolver"
DELAY = 0.35

GENE_QUERIES = {
    "bacteria": '("16S ribosomal RNA"[Gene] OR "16S rRNA"[Title])',
    "viruses": '("complete genome"[Title] OR "polyprotein"[Gene])',
    "fungi": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])',
    "oomycetes": '("internal transcribed spacer"[Title] OR "ITS1"[Title] OR "ITS2"[Title] OR "5.8S"[Title])',
}


def norm_name(value):
    value = value.lower().replace("（", "(").replace("）", ")")
    value = re.sub(r"[^a-z0-9. ]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def eutils_get(endpoint, params):
    url = f"{ENTREZ_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
    last_error = None
    for attempt in range(5):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": f"MaizePathogenDB/{TOOL} ({EMAIL})"},
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.read().decode("utf-8", errors="replace")
        except Exception as exc:
            last_error = exc
            body = ""
            if hasattr(exc, "read"):
                try:
                    body = exc.read().decode("utf-8", errors="replace")[:200]
                except Exception:
                    pass
            print(f"eutils retry {attempt + 1}: {url} -> {exc} {body}", file=__import__("sys").stderr)
            time.sleep(4 + attempt * 4)
    raise RuntimeError(f"eutils request failed: {url}: {last_error}")


def esearch(db, term):
    text = eutils_get("esearch.fcgi", {
        "db": db, "term": term, "retmode": "json", "retmax": "10",
        "email": EMAIL, "tool": TOOL,
    })
    time.sleep(DELAY)
    return json.loads(text)["esearchresult"].get("idlist", [])


def esummary(db, ids):
    text = eutils_get("esummary.fcgi", {
        "db": db, "id": ",".join(ids), "retmode": "json",
        "email": EMAIL, "tool": TOOL,
    })
    time.sleep(DELAY)
    data = json.loads(text)["result"]
    return {uid: data[uid] for uid in data["uids"]}


def resolve_taxonomy_name(name):
    name = name.strip()
    if not name:
        return None
    ids = esearch("taxonomy", f'"{name}"[All Names]')
    if not ids:
        return None
    summaries = esummary("taxonomy", ids)
    target = norm_name(name)
    exact = [summaries[i] for i in ids if norm_name(summaries[i].get("scientificname", "")) == target]
    if exact:
        return exact[0]
    partial = [summaries[i] for i in ids
               if target in norm_name(summaries[i].get("scientificname", ""))
               or norm_name(summaries[i].get("scientificname", "")) in target]
    return partial[0] if partial else None


def lineage_species_node(taxid):
    text = eutils_get("efetch.fcgi", {
        "db": "taxonomy", "id": str(taxid), "retmode": "xml",
        "email": EMAIL, "tool": TOOL,
    })
    time.sleep(DELAY)
    root = ET.fromstring(text)
    taxa = root.findall(".//Taxon")
    species = None
    for taxon in taxa:
        rank = (taxon.findtext("Rank") or "").strip()
        if rank in ("species", "no rank"):
            species = {
                "taxid": taxon.findtext("TaxId", "").strip(),
                "name": taxon.findtext("ScientificName", "").strip(),
                "rank": rank,
            }
    return species


def count_marker_sequences(taxid, category):
    query = f'txid{taxid}[Organism] AND {GENE_QUERIES[category]}'
    ids = esearch("nucleotide", query)
    if ids or category == "viruses":
        return len(ids), ids[:3], taxid
    species_node = lineage_species_node(taxid)
    if species_node and species_node["taxid"] != str(taxid):
        species_taxid = species_node["taxid"]
        query2 = f'txid{species_taxid}[Organism] AND {GENE_QUERIES[category]}'
        ids2 = esearch("nucleotide", query2)
        return len(ids2), ids2[:3], species_taxid
    return 0, [], taxid


def name_variants(value):
    value = value.replace("（", "(").replace("）", ")")
    parts = re.split(r"\s*/\s*", value)
    variants = []
    for part in parts:
        part = re.sub(r"\s*\(.*\)\s*", " ", part).strip()
        if part and part not in variants:
            variants.append(part)
    return variants


def main():
    with open(SPECIES_LIST, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    rows_by_taxid = {r["taxid"]: r for r in rows}
    results = []

    # 1) NOT_VERIFIED rows
    for row in rows:
        if row["taxid"] != "NOT_VERIFIED":
            continue
        variants = name_variants(row["full_name"] or row["species"])
        best = None
        best_name = None
        for variant in variants:
            summary = resolve_taxonomy_name(variant)
            if summary:
                best = summary
                best_name = variant
                break
        if best:
            results.append({
                "row_type": "not_verified",
                "original_species": row["species"],
                "original_taxid": "NOT_VERIFIED",
                "category": row["category"],
                "candidate_names": "; ".join(variants),
                "matched_name": best_name,
                "resolved_taxid": best["taxid"],
                "resolved_scientific_name": best["scientificname"],
                "resolved_rank": best.get("rank", ""),
                "status": "RESOLVED",
                "sequence_hits": "",
                "notes": "",
            })
        else:
            results.append({
                "row_type": "not_verified",
                "original_species": row["species"],
                "original_taxid": "NOT_VERIFIED",
                "category": row["category"],
                "candidate_names": "; ".join(variants),
                "matched_name": "",
                "resolved_taxid": "",
                "resolved_scientific_name": "",
                "resolved_rank": "",
                "status": "PENDING",
                "sequence_hits": "",
                "notes": "no exact NCBI Taxonomy match for any candidate name",
            })

    # 2) Fijivirus zeae duplicate
    for taxid in ("10989", "10990"):
        row = rows_by_taxid[taxid]
        node = lineage_species_node(taxid)
        if node:
            results.append({
                "row_type": "duplicate_species_name",
                "original_species": row["species"],
                "original_taxid": taxid,
                "category": row["category"],
                "candidate_names": row["full_name"],
                "matched_name": "",
                "resolved_taxid": node["taxid"],
                "resolved_scientific_name": node["name"],
                "resolved_rank": node["rank"],
                "status": "RESOLVED",
                "sequence_hits": "",
                "notes": f"isolate taxid {taxid}; species-level node from NCBI lineage",
            })
        else:
            results.append({
                "row_type": "duplicate_species_name",
                "original_species": row["species"],
                "original_taxid": taxid,
                "category": row["category"],
                "candidate_names": row["full_name"],
                "matched_name": "",
                "resolved_taxid": "",
                "resolved_scientific_name": "",
                "resolved_rank": "",
                "status": "PENDING",
                "sequence_hits": "",
                "notes": "no species-level node found",
            })

    # 3) TaxID present but no sequence in the release
    for row in rows:
        if row["taxid"] == "NOT_VERIFIED" or row["sequence_status"] == "present":
            continue
        n_hits, sample_ids, query_taxid = count_marker_sequences(row["taxid"], row["category"])
        results.append({
            "row_type": "missing_sequence",
            "original_species": row["species"],
            "original_taxid": row["taxid"],
            "category": row["category"],
            "candidate_names": "",
            "matched_name": "",
            "resolved_taxid": row["taxid"],
            "resolved_scientific_name": row["species"],
            "resolved_rank": "",
            "status": "SEQUENCE_FOUND" if n_hits else "NO_SEQUENCE",
            "sequence_hits": str(n_hits),
            "notes": f"query_taxid={query_taxid}; marker hits: " + ", ".join(sample_ids),
        })

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    columns = [
        "row_type", "original_species", "original_taxid", "category",
        "candidate_names", "matched_name", "resolved_taxid",
        "resolved_scientific_name", "resolved_rank", "status",
        "sequence_hits", "notes",
    ]
    with open(OUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("\t".join(columns) + "\n")
        for result in results:
            fh.write("\t".join(str(result[c]) for c in columns) + "\n")

    from collections import Counter
    print(f"wrote {OUT_TSV} ({len(results)} rows)")
    print("status:", dict(Counter(r["status"] for r in results)))
    print("by row_type:", dict(Counter(r["row_type"] for r in results)))


if __name__ == "__main__":
    main()
