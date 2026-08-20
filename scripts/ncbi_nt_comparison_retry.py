#!/usr/bin/env python3
"""Retry timed-out NCBI-nt web BLAST queries for the 150-query comparison."""

import json
import os
import random
import re
import time
from collections import defaultdict

import requests

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
FASTA = os.path.join(BASE, "sequences", "maize_pathogens_all.fasta")
RESULTS_FILE = os.path.join(BASE, "docs", "validation", "ncbi_nt_comparison.json")

SAMPLE_N = 150
RANDOM_SEED = 42
NCBI_DELAY = 3.5
NCBI_TIMEOUT = 1800
NCBI_POLL = 5
MAX_ATTEMPTS = 5
NCBI_EMAIL = "maize_pathogen_db@example.com"
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CAT_ORDER = ["bacteria", "viruses", "fungi"]


def parse_fasta(fpath):
    records = []
    current = None
    for line in open(fpath):
        if line.startswith(">"):
            if current:
                records.append(current)
            parts = line[1:].strip().split("|", 3)
            taxid = parts[0]
            species = parts[1] if len(parts) > 1 else ""
            category = parts[2] if len(parts) > 2 else "unknown"
            current = {"header": line[1:].strip(), "seq": "",
                       "taxid": taxid, "species": species,
                       "category": category}
        elif current is not None:
            current["seq"] += line.strip()
    if current:
        records.append(current)
    return records


def reproduce_sample():
    random.seed(RANDOM_SEED)
    all_records = parse_fasta(FASTA)
    by_cat_species = defaultdict(list)
    for r in all_records:
        if r["category"] in CAT_ORDER:
            by_cat_species[(r["category"], r["species"])].append(r)

    total_species = sum(len([k for k in by_cat_species if k[0] == cat]) for cat in CAT_ORDER)
    sample = {}
    for cat in CAT_ORDER:
        cat_species = [k for k in by_cat_species if k[0] == cat]
        prop = int(round(len(cat_species) / total_species * SAMPLE_N))
        chosen = random.sample(cat_species, min(prop, len(cat_species)))
        for k in chosen:
            sample[by_cat_species[k][0]["header"]] = random.choice(by_cat_species[k])
    return sample


def submit_ncbi_blast(seq, seq_id="query"):
    params = {
        "CMD": "Put", "PROGRAM": "blastn", "DATABASE": "nt",
        "QUERY": f">{seq_id}\n{seq}",
        "HITLIST_SIZE": "3", "MEGABLAST": "on", "FILTER": "L",
    }
    try:
        r = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=params, timeout=60)
        if r.status_code == 429:
            time.sleep(5)
            r = requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data=params, timeout=60)
        r.raise_for_status()
    except Exception as e:
        return f"SUBMIT_ERROR: {e}"

    rid = None
    for line in r.text.split("\n"):
        m = re.search(r"RID\s*=\s*['\"]?([A-Z0-9-]+)['\"]?", line)
        if m:
            rid = m.group(1)
            break
    if not rid:
        return "NO_RID"

    start = time.time()
    while time.time() - start < NCBI_TIMEOUT:
        time.sleep(NCBI_POLL)
        try:
            r2 = requests.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                              params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"}, timeout=30)
            r2.raise_for_status()
        except Exception as e:
            return f"POLL_ERROR: {e}"
        if "Status=WAITING" in r2.text or "Status=UNKNOWN" in r2.text:
            continue
        if "Status=FAILED" in r2.text:
            return "BLAST_FAILED"
        if r2.text.startswith("<?xml") or "<BlastOutput" in r2.text:
            return r2.text
    return "TIMEOUT"


def extract_top_hit(xml_data):
    if isinstance(xml_data, str) and xml_data.startswith(("BLAST_ERROR", "NO_", "TIMEOUT", "SUBMIT_", "POLL_")):
        return {"description": xml_data, "pident": 0.0, "accession": ""}
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(xml_data)
        hits = root.findall(".//Hit")
        if not hits:
            return {"description": "NO_HITS", "pident": 0.0, "accession": ""}
        top = hits[0]
        desc = top.findtext("Hit_def", "UNKNOWN")
        accession = top.findtext("Hit_accession", "")
        if not accession:
            accession = desc.split()[0] if desc else ""
        hsps = top.findall(".//Hsp")
        if not hsps:
            return {"description": desc, "pident": 0.0, "accession": accession}
        identity = int(hsps[0].findtext("Hsp_identity", "0"))
        align_len = int(hsps[0].findtext("Hsp_align-len", "1"))
        pident = identity / align_len * 100 if align_len else 0.0
        return {"description": desc, "pident": pident, "accession": accession}
    except ET.ParseError as e:
        return {"description": f"XML_PARSE_ERROR: {e}", "pident": 0.0, "accession": ""}
    except Exception as e:
        return {"description": f"XML_ERROR: {e}", "pident": 0.0, "accession": ""}


def fetch_ncbi_taxid(accession):
    if not accession:
        return ""
    params = {
        "db": "nucleotide", "id": accession, "retmode": "json",
        "email": NCBI_EMAIL, "tool": "maize_pathogen_db",
    }
    try:
        r = requests.get(f"{ENTREZ_BASE}/esummary.fcgi", params=params, timeout=30)
        if r.status_code == 429:
            time.sleep(3)
            r = requests.get(f"{ENTREZ_BASE}/esummary.fcgi", params=params, timeout=30)
        r.raise_for_status()
        data = r.json().get("result", {})
        uids = data.get("uids", [])
        if uids:
            return str(data[uids[0]].get("taxid", ""))
    except Exception:
        pass
    return ""


def main():
    sample = reproduce_sample()
    sample_by_key = {}
    for rec in sample.values():
        sample_by_key[(rec["category"], rec["species"])] = rec

    with open(RESULTS_FILE) as f:
        data = json.load(f)

    retried = 0
    for result in data["results"]:
        desc = str(result.get("ncbi_top", ""))
        is_failed = (not desc) or desc.startswith(("TIMEOUT", "POLL_", "SUBMIT_", "BLAST_FAILED", "NO_RID"))
        if not is_failed:
            continue

        rec = sample_by_key.get((result["category"], result["species"]))
        if rec is None:
            print(f"WARN: no sampled sequence for {result['species']}", flush=True)
            continue

        print(f"Retrying: {result['species'][:60]} ({result['category']})", end=" ", flush=True)
        hit = None
        for attempt in range(MAX_ATTEMPTS):
            blast_result = submit_ncbi_blast(rec["seq"])
            hit = extract_top_hit(blast_result)
            if not str(hit["description"]).startswith(("TIMEOUT", "POLL_", "SUBMIT_", "BLAST_FAILED", "NO_RID", "XML_")):
                break
            if attempt < MAX_ATTEMPTS - 1:
                print(f"attempt{attempt + 1} ", end="", flush=True)
                time.sleep(10)

        result["ncbi_pident"] = hit["pident"]
        result["ncbi_top"] = hit["description"]
        result["ncbi_accession"] = hit.get("accession", "")
        result["ncbi_taxid"] = fetch_ncbi_taxid(result["ncbi_accession"])
        result["ncbi_correct"] = bool(result["ncbi_taxid"]) and result["ncbi_taxid"] == rec["taxid"]
        status = "OK" if result["ncbi_correct"] else "FAIL"
        print(f"{status} (pident={hit['pident']:.1f}%)", flush=True)
        retried += 1

        with open(RESULTS_FILE, "w") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        time.sleep(NCBI_DELAY)

    summary = {}
    for cat in CAT_ORDER:
        cat_results = [r for r in data["results"] if r["category"] == cat]
        n = len(cat_results)
        db_c = sum(1 for r in cat_results if r.get("db_correct"))
        ncbi_c = sum(1 for r in cat_results if r.get("ncbi_correct"))
        summary[cat] = {"n": n, "db_correct": db_c, "ncbi_correct": ncbi_c,
                        "db_accuracy": round(db_c / n * 100, 1) if n else 0,
                        "ncbi_accuracy": round(ncbi_c / n * 100, 1) if n else 0}
    n_total = len(data["results"])
    db_total = sum(1 for r in data["results"] if r.get("db_correct"))
    ncbi_total = sum(1 for r in data["results"] if r.get("ncbi_correct"))
    summary["overall"] = {"n": n_total, "db_correct": db_total, "ncbi_correct": ncbi_total,
                          "db_accuracy": round(db_total / n_total * 100, 1),
                          "ncbi_accuracy": round(ncbi_total / n_total * 100, 1)}
    summary["retried"] = retried
    data.update(summary)

    with open(RESULTS_FILE, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print(f"\nRetried: {retried}")
    for cat in CAT_ORDER:
        print(f"  {cat}: NCBI-nt {summary[cat]['ncbi_correct']}/{summary[cat]['n']} "
              f"({summary[cat]['ncbi_accuracy']:.1f}%)")
    print(f"  overall: {ncbi_total}/{n_total} ({summary['overall']['ncbi_accuracy']:.1f}%)")


if __name__ == "__main__":
    main()
