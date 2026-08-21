#!/usr/bin/env python3
"""Classify benchmark queries with MaizePathogenDB and NCBI-nt; compute classification metrics."""

import json
import os
import random
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from Bio import SeqIO

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
OUT = os.path.join(BASE, "docs", "validation", "classification_benchmark")
POS_FASTA = os.path.join(BASE, "sequences", "maize_pathogens_all.fasta")
NEG_FASTA = os.path.join(OUT, "negative_queries.fasta")
NEG_META = os.path.join(OUT, "negative_queries_meta.json")
BLAST_DB = os.path.join(BASE, "blast_db", "maize_pathogens_all")
BLAST_BIN = "/Users/gfgao/Desktop/blacksoil_metaG/tools/ncbi-blast-2.17.0+/bin/blastn"
BACKGROUND_FASTA = os.path.join(OUT, "negative_queries.fasta")
BACKGROUND_DB = os.path.join(OUT, "background_negatives")
CATALOG = json.load(open("/Users/gfgao/Desktop/blacksoil_metaG/Figshare/taxonomy.json"))
NCBI_150 = os.path.join(BASE, "docs", "validation", "ncbi_nt_comparison.json")

OURDB_FILE = os.path.join(OUT, "benchmark_ourdb_predictions.json")
NCBI_FILE = os.path.join(OUT, "benchmark_ncbi_predictions.json")
METRICS_FILE = os.path.join(OUT, "benchmark_metrics.json")

NCBI_DELAY = 2.0
NCBI_TIMEOUT = 300
NCBI_POLL = 5
NCBI_EMAIL = "maize_pathogen_db@example.com"
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MAX_ATTEMPTS = 1
WORKERS = 6

PID_THRESHOLD = 97.0
QCOV_THRESHOLD = 90.0
SPECIES_PID = 99.5
SPECIES_QCOV = 99.0
BACKGROUND_PID_GAP = 0.5
BACKGROUND_QCOV = 90.0
CAT_ORDER = ["bacteria", "viruses", "fungi", "oomycetes"]


def normalize_seq(seq):
    return re.sub(r"[^ACGTNacgtn]", "", seq).upper()


def parse_fasta(fpath):
    records = []
    current = None
    for line in open(fpath):
        if line.startswith(">"):
            if current:
                records.append(current)
            parts = line[1:].strip().split("|", 3)
            current = {"header": line[1:].strip(), "seq": "",
                       "taxid": parts[0], "species": parts[1] if len(parts) > 1 else "",
                       "category": parts[2] if len(parts) > 2 else "unknown"}
        elif current is not None:
            current["seq"] += line.strip()
    if current:
        records.append(current)
    return records


def group_of(rec, taxid_to_group):
    return taxid_to_group.get(rec["taxid"], rec["category"] if rec["category"] in ("bacteria", "viruses") else "fungi")


def load_queries():
    taxid_to_group = {}
    for r in CATALOG:
        g = "oomycetes" if "Oomycota" in (r.get("phylum") or "") else r["category"].lower()
        if g == "virus":
            g = "viruses"
        taxid_to_group[str(r["taxid"])] = g

    queries = []
    for rec in parse_fasta(POS_FASTA):
        queries.append({
            "qid": f"p{len(queries) + 1:04d}", "header": rec["header"], "seq": rec["seq"],
            "taxid": rec["taxid"], "species": rec["species"],
            "category": group_of(rec, taxid_to_group), "label": 1,
        })

    neg_meta = json.load(open(NEG_META))
    neg_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(NEG_FASTA, "fasta")}
    for m in neg_meta:
        queries.append({
            "qid": m["qid"], "header": m["qid"], "seq": neg_seqs[m["qid"]],
            "taxid": m["taxid"], "species": m["species"],
            "category": m["category"], "label": 0,
        })
    return queries


def run_ourdb(queries):
    qmap = {q["qid"]: q for q in queries}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for q in queries:
            f.write(f">{q['qid']}\n{q['seq']}\n")
        qfile = f.name
    cmd = (f'"{BLAST_BIN}" -query "{qfile}" -db "{BLAST_DB}" '
           f'-outfmt "6 qseqid sseqid pident length qcovs" -max_target_seqs 1 -num_threads 4')
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=1800)
    os.unlink(qfile)
    best = {}
    for line in proc.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        qid = parts[0]
        if qid not in best:
            best[qid] = {
                "sseqid": parts[1], "pident": float(parts[2]),
                "length": int(parts[3]), "qcovs": float(parts[4]),
            }
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for q in queries:
            f.write(f">{q['qid']}\n{q['seq']}\n")
        bgfile = f.name
    cmd_bg = (f'"{BLAST_BIN}" -query "{bgfile}" -db "{BACKGROUND_DB}" '
              f'-outfmt "6 qseqid sseqid pident length qcovs bitscore evalue" '
              f'-max_target_seqs 1 -num_threads 4')
    proc_bg = subprocess.run(cmd_bg, shell=True, capture_output=True, text=True, timeout=1800)
    os.unlink(bgfile)
    bg_best = {}
    for line in proc_bg.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        qid = parts[0]
        if qid not in bg_best:
            bg_best[qid] = {
                "sseqid": parts[1], "pident": float(parts[2]),
                "length": int(parts[3]), "qcovs": float(parts[4]),
                "bitscore": float(parts[5]), "evalue": parts[6],
            }
    results = []
    for q in queries:
        hit = best.get(q["qid"])
        bg = bg_best.get(q["qid"])
        candidate = bool(hit and hit["pident"] >= SPECIES_PID and hit["qcovs"] >= SPECIES_QCOV)
        rejected = bool(bg and hit and bg["pident"] >= hit["pident"] - BACKGROUND_PID_GAP
                        and bg["qcovs"] >= BACKGROUND_QCOV)
        pred = candidate and not rejected
        results.append({
            "qid": q["qid"], "label": q["label"], "category": q["category"],
            "predicted": pred, "hit": hit, "background_hit": bg,
        })
    with open(OURDB_FILE, "w") as f:
        json.dump({"results": results}, f, indent=2)
    return results


def reproduce_sample():
    random.seed(42)
    records = parse_fasta(POS_FASTA)
    all_recs = [r for r in records if r["category"] in ("bacteria", "viruses", "fungi")]
    by_cat_species = defaultdict(list)
    for r in all_recs:
        by_cat_species[(r["category"], r["species"])].append(r)
    total_species = sum(len({k for k in by_cat_species if k[0] == cat}) for cat in ["bacteria", "viruses", "fungi"])
    sample = {}
    for cat in ["bacteria", "viruses", "fungi"]:
        cat_species = [k for k in by_cat_species if k[0] == cat]
        prop = int(round(len(cat_species) / total_species * 150))
        chosen = random.sample(cat_species, min(prop, len(cat_species)))
        for k in chosen:
            sample[by_cat_species[k][0]["header"]] = random.choice(by_cat_species[k])
    return list(sample.values())


def submit_ncbi_blast(seq, seq_id="query"):
    params = {
        "CMD": "Put", "PROGRAM": "blastn", "DATABASE": "nt",
        "QUERY": f">{seq_id}\n{seq}",
        "HITLIST_SIZE": "1", "MEGABLAST": "on", "FILTER": "L",
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


def classify_one(q, catalog_taxids, precomputed):
    if q["qid"] in precomputed:
        return precomputed[q["qid"]]
    hit = None
    for attempt in range(MAX_ATTEMPTS):
        blast_result = submit_ncbi_blast(q["seq"], q["qid"])
        hit = extract_top_hit(blast_result)
        if not str(hit["description"]).startswith(("TIMEOUT", "POLL_", "SUBMIT_", "BLAST_FAILED", "NO_RID", "XML_")):
            break
        if attempt < MAX_ATTEMPTS - 1:
            time.sleep(10)
    accession = hit.get("accession", "")
    taxid = fetch_ncbi_taxid(accession)
    predicted = bool(taxid) and hit["pident"] >= SPECIES_PID and taxid in catalog_taxids
    return {
        "qid": q["qid"], "label": q["label"], "category": q["category"],
        "predicted": predicted, "accession": accession, "taxid": taxid,
        "pident": hit["pident"], "top": hit["description"],
    }


def run_ncbi(queries):
    catalog_taxids = set()
    for r in CATALOG:
        catalog_taxids.add(str(r["taxid"]))

    precomputed = {}
    if os.path.exists(NCBI_FILE):
        precomputed = {r["qid"]: r for r in json.load(open(NCBI_FILE))["results"]}
    if os.path.exists(NCBI_150):
        sample_list = reproduce_sample()
        old = json.load(open(NCBI_150))["results"]
        for rec, r in zip(sample_list, old):
            q = next((x for x in queries if x["header"] == rec["header"]), None)
            if q is not None:
                taxid = r.get("ncbi_taxid", "")
                precomputed[q["qid"]] = {
                    "qid": q["qid"], "label": q["label"], "category": q["category"],
                    "predicted": bool(taxid) and r.get("ncbi_pident", 0) >= SPECIES_PID and taxid in catalog_taxids,
                    "accession": r.get("ncbi_accession", ""), "taxid": taxid,
                    "pident": r.get("ncbi_pident", 0), "top": r.get("ncbi_top", ""),
                }

    todo = [q for q in queries if q["qid"] not in precomputed]
    print(f"NCBI-nt: {len(precomputed)} precomputed, {len(todo)} to run", flush=True)
    results = [precomputed[q["qid"]] for q in queries if q["qid"] in precomputed]
    done = set(precomputed.keys())

    def worker(q):
        return classify_one(q, catalog_taxids, {})

    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futures = {ex.submit(worker, q): q for q in todo}
        for i, fut in enumerate(as_completed(futures), 1):
            q = futures[fut]
            try:
                res = fut.result()
                results.append(res)
                done.add(q["qid"])
                status = "OK" if res["predicted"] else "NEG"
                print(f"  [{i}/{len(todo)}] {q['species'][:55]} ({q['category']}) {status}", flush=True)
            except Exception as e:
                print(f"  [{i}/{len(todo)}] {q['species'][:55]} ERROR {e}", flush=True)
            if i % 10 == 0:
                with open(NCBI_FILE, "w") as f:
                    json.dump({"results": results}, f, indent=2, ensure_ascii=False)

    with open(NCBI_FILE, "w") as f:
        json.dump({"results": results}, f, indent=2, ensure_ascii=False)
    return results


def metrics(predictions, label):
    tp = sum(1 for r in predictions if r["label"] == 1 and r["predicted"])
    fn = sum(1 for r in predictions if r["label"] == 1 and not r["predicted"])
    fp = sum(1 for r in predictions if r["label"] == 0 and r["predicted"])
    tn = sum(1 for r in predictions if r["label"] == 0 and not r["predicted"])
    sensitivity = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else 0.0
    balanced = (sensitivity + specificity) / 2
    return {
        "n_pos": tp + fn, "n_neg": tn + fp, "tp": tp, "fn": fn, "fp": fp, "tn": tn,
        "sensitivity": round(sensitivity * 100, 1),
        "specificity": round(specificity * 100, 1),
        "f1": round(f1 * 100, 1),
        "balanced_accuracy": round(balanced * 100, 1),
    }


def compute_metrics():
    ourdb = json.load(open(OURDB_FILE))["results"]
    ncbi = json.load(open(NCBI_FILE))["results"]
    catalog_taxids = {str(r["taxid"]) for r in CATALOG}
    for r in ncbi:
        r["predicted"] = bool(r.get("taxid")) and r.get("pident", 0) >= SPECIES_PID and r.get("taxid") in catalog_taxids
    out = {"ourdb": {}, "ncbi_nt": {}}
    for name, rows in (("ourdb", ourdb), ("ncbi_nt", ncbi)):
        out[name]["overall"] = metrics(rows, name)
        for cat in CAT_ORDER:
            cat_rows = [r for r in rows if r["category"] == cat]
            out[name][cat] = metrics(cat_rows, name)
    with open(METRICS_FILE, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    for name in ("ourdb", "ncbi_nt"):
        print(f"\n{name}:")
        for cat in ["overall"] + CAT_ORDER:
            m = out[name][cat]
            print(f"  {cat}: Sens {m['sensitivity']:.1f}% Spec {m['specificity']:.1f}% "
                  f"F1 {m['f1']:.1f}% Balanced {m['balanced_accuracy']:.1f}%")


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "metrics"
    os.makedirs(OUT, exist_ok=True)
    queries = load_queries()
    if phase == "local":
        run_ourdb(queries)
    elif phase == "ncbi":
        run_ncbi(queries)
    elif phase == "metrics":
        compute_metrics()
