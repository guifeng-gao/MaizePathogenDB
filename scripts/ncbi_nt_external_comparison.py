#!/usr/bin/env python3
"""Fair head-to-head comparison: External 2020-2026 sequences vs NCBI-nt (species-level TaxID)."""

import json
import os
import re
import time
from collections import Counter

import requests
from Bio import SeqIO

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
OUT_DIR = os.path.join(BASE, "docs", "validation")
META_FILE = os.path.join(OUT_DIR, "external", "external_queries_meta.json")
FASTA_FILE = os.path.join(OUT_DIR, "external", "external_queries.fasta")
DB_FINAL = os.path.join(OUT_DIR, "external", "external_validation_final.json")
RESULTS_FILE = os.path.join(OUT_DIR, "ncbi_nt_external_comparison.json")

NCBI_DELAY = 3.5
NCBI_TIMEOUT = 600
NCBI_POLL = 5
NCBI_EMAIL = "maize_pathogen_db@example.com"
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_ATTEMPTS = 3

CAT_ORDER = ["bacteria", "viruses", "fungi", "oomycetes"]


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
    os.makedirs(OUT_DIR, exist_ok=True)
    meta = json.load(open(META_FILE))
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(FASTA_FILE, "fasta")}
    db = json.load(open(DB_FINAL))
    db_by_qid = {d["qid"]: d for d in db["details"]}

    results = []
    if os.path.exists(RESULTS_FILE):
        results = json.load(open(RESULTS_FILE)).get("details", [])
    done = {r["qid"] for r in results}

    for i, m in enumerate(meta, 1):
        qid = m["qid"]
        if qid in done:
            print(f"  [{i}/{len(meta)}] Skipped (cached): {m['species'][:50]}", flush=True)
            continue

        seq = seqs.get(qid)
        if not seq:
            print(f"  [{i}/{len(meta)}] WARN: no sequence for {qid}", flush=True)
            continue

        print(f"  [{i}/{len(meta)}] BLASTing: {m['species'][:55]} ({m['category']}) ...", end=" ", flush=True)
        hit = None
        for attempt in range(MAX_ATTEMPTS):
            blast_result = submit_ncbi_blast(seq, qid)
            hit = extract_top_hit(blast_result)
            if not str(hit["description"]).startswith(("TIMEOUT", "POLL_", "SUBMIT_", "BLAST_FAILED", "NO_RID", "XML_")):
                break
            if attempt < MAX_ATTEMPTS - 1:
                print(f"attempt{attempt + 1} ", end="", flush=True)
                time.sleep(10)

        accession = hit.get("accession", "")
        ncbi_taxid = fetch_ncbi_taxid(accession)
        db_row = db_by_qid.get(qid, {})
        result = {
            "qid": qid,
            "taxid": m["taxid"],
            "species": m["species"],
            "category": m["category"],
            "accession": m["accession"],
            "db_correct": db_row.get("correct"),
            "db_top_taxid": db_row.get("top_taxid"),
            "db_pident": db_row.get("pident"),
            "ncbi_accession": accession,
            "ncbi_taxid": ncbi_taxid,
            "ncbi_correct": bool(ncbi_taxid) and ncbi_taxid == m["taxid"],
            "ncbi_pident": hit["pident"],
            "ncbi_top": hit["description"],
        }
        results.append(result)
        with open(RESULTS_FILE, "w") as f:
            json.dump({"details": results}, f, indent=2, ensure_ascii=False)
        status = "OK" if result["ncbi_correct"] else "FAIL"
        print(f"{status} (pident={hit['pident']:.1f}%)", flush=True)
        time.sleep(NCBI_DELAY)

    summary = {}
    for cat in CAT_ORDER:
        rows = [r for r in results if r["category"] == cat]
        n = len(rows)
        db_c = sum(1 for r in rows if r.get("db_correct"))
        ncbi_c = sum(1 for r in rows if r.get("ncbi_correct"))
        summary[cat] = {
            "n_queries": n,
            "db_correct": db_c,
            "db_accuracy": round(db_c / n * 100, 1) if n else 0.0,
            "ncbi_correct": ncbi_c,
            "ncbi_accuracy": round(ncbi_c / n * 100, 1) if n else 0.0,
        }
    n_total = len(results)
    db_total = sum(1 for r in results if r.get("db_correct"))
    ncbi_total = sum(1 for r in results if r.get("ncbi_correct"))
    summary["overall"] = {
        "n_queries": n_total,
        "db_correct": db_total,
        "db_accuracy": round(db_total / n_total * 100, 1) if n_total else 0.0,
        "ncbi_correct": ncbi_total,
        "ncbi_accuracy": round(ncbi_total / n_total * 100, 1) if n_total else 0.0,
    }

    final = {
        "validation_type": "external_ncbi_nt_comparison",
        "description": "Fair comparison on 2020-2026 external sequences: MaizePathogenDB vs NCBI-nt, species-level TaxID",
        "run_date": "2026-08-20",
        "total_external_queries": n_total,
        "results": summary,
        "details": results,
    }
    with open(RESULTS_FILE, "w") as f:
        json.dump(final, f, indent=2, ensure_ascii=False)

    print("\nFinal (species-level TaxID):")
    for cat in CAT_ORDER:
        s = summary[cat]
        print(f"  {cat}: OurDB {s['db_correct']}/{s['n_queries']} ({s['db_accuracy']:.1f}%) | "
              f"NCBI-nt {s['ncbi_correct']}/{s['n_queries']} ({s['ncbi_accuracy']:.1f}%)")
    o = summary["overall"]
    print(f"  overall: OurDB {o['db_correct']}/{o['n_queries']} ({o['db_accuracy']:.1f}%) | "
          f"NCBI-nt {o['ncbi_correct']}/{o['n_queries']} ({o['ncbi_accuracy']:.1f}%)")


if __name__ == "__main__":
    main()
