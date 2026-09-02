#!/usr/bin/env python3
"""Classify fungal and oomycete queries against fixed NCBI ITS databases."""

import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_validation as rv

ROOT = rv.ROOT
QUERY_DIR = rv.QUERY_DIR
BLAST_DIR = os.environ.get("BLAST_DIR", rv.BLAST_DIR)
RESULT_DIR = os.environ.get("RESULT_DIR", rv.RESULT_DIR)
ITS_EUK_DB = os.path.join(ROOT, "260samples_fungi", "analysis", "db",
                          "ncbi_ITS_eukaryote", "ITS_eukaryote_sequences")
ITS_REFSEQ_DB = os.path.join(ROOT, "260samples_fungi", "analysis", "db",
                             "ncbi_ITS_RefSeq", "ITS_RefSeq_Fungi")

SPECIES_TH = (99.0, 90.0)
GENUS_TH = (95.0, 70.0)


def load_queries(fasta, meta_path):
    queries = []
    meta = {}
    with open(meta_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["qid"]] = row
    for record in rv.read_fasta(fasta):
        info = rv.parse_query_header(record["header"])
        m = meta.get(info["qid"], {})
        queries.append({**info, "species": m.get("species", info.get("species", ""))})
    return queries


def run_blast_taxids(query_file, db):
    cmd = [
        rv.BLASTN, "-query", query_file, "-db", db,
        "-outfmt", "6 qseqid sseqid pident qcovs staxids",
        "-max_target_seqs", "1", "-evalue", "1e-5", "-num_threads", "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    hits = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        qid = parts[0].split("|", 1)[0]
        staxid = parts[4].strip()
        if not staxid or staxid == "N/A":
            continue
        hits[qid] = {"pident": float(parts[2]), "qcovs": float(parts[3]),
                     "staxid": staxid}
    return hits


def classify(pos_hits, neg_hits, pos, neg, pid, qcov, level,
             catalog_taxids, nodes, names):
    def get_hit_taxid(hit):
        if "staxid" in hit:
            return hit["staxid"]
        return rv.subject_taxid(hit.get("sseqid", ""))

    tp = fp = fn = tn = 0
    for record in pos:
        hit = pos_hits.get(record["qid"])
        if not hit or hit["pident"] < pid or hit["qcovs"] < qcov:
            fn += 1
            continue
        hit_taxid = get_hit_taxid(hit)
        if level == "species":
            ok = rv.current_taxid(hit_taxid) == rv.current_taxid(record["taxid"])
        else:
            ok = (rv.norm_taxon(rv.genus_of(hit_taxid, nodes, names))
                  == rv.norm_taxon(rv.genus_of(record["taxid"], nodes, names)))
        if ok:
            tp += 1
        else:
            fn += 1
    for record in neg:
        hit = neg_hits.get(record["qid"])
        if not hit or hit["pident"] < pid or hit["qcovs"] < qcov:
            tn += 1
            continue
        hit_taxid = rv.current_taxid(get_hit_taxid(hit))
        if level == "species":
            positive = hit_taxid in catalog_taxids
        else:
            hit_genus = rv.norm_taxon(rv.genus_of(hit_taxid, nodes, names))
            positive = any(hit_genus == rv.norm_taxon(rv.genus_of(t, nodes, names))
                           for t in catalog_taxids)
        if positive:
            fp += 1
        else:
            tn += 1
    return rv.metrics(tp, fp, fn, tn)


def main():
    rv.load_merged()
    nodes, names = rv.load_taxdump()
    catalog_taxids = {rv.current_taxid(str(r.get("taxid"))) for r in rv.load_taxonomy()}
    pos_file = os.path.join(QUERY_DIR, "external_positives.fasta")
    neg_file = os.path.join(QUERY_DIR, "negatives.fasta")
    pos = [p for p in load_queries(pos_file, os.path.join(QUERY_DIR, "external_positives_meta.tsv"))
           if p["category"] in ("fungi", "oomycetes")]
    neg = [n for n in load_queries(neg_file, os.path.join(QUERY_DIR, "negatives_meta.tsv"))
           if n["category"] in ("fungi", "oomycetes")]

    # MPDB comparison on the same subset
    mpdb_hits = rv.run_blast(pos_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    mpdb_neg = rv.run_blast(neg_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    db_results = {"MPDB": {"hits_pos": mpdb_hits, "hits_neg": mpdb_neg}}
    for name, db in (("NCBI_ITS_eukaryote", ITS_EUK_DB),
                     ("NCBI_ITS_RefSeq_Fungi", ITS_REFSEQ_DB)):
        db_results[name] = {
            "hits_pos": run_blast_taxids(pos_file, db),
            "hits_neg": run_blast_taxids(neg_file, db),
        }

    results = {"categories": ["fungi", "oomycetes"], "databases": {}}
    for db_name, hits in db_results.items():
        entry = {"per_category": {}, "overall": {}}
        for category in ("fungi", "oomycetes"):
            pos_sub = [p for p in pos if p["category"] == category]
            neg_sub = [n for n in neg if n["category"] == category]
            entry["per_category"][category] = {
                "species": classify(hits["hits_pos"], hits["hits_neg"],
                                    pos_sub, neg_sub, *SPECIES_TH, "species",
                                    catalog_taxids, nodes, names),
                "genus": classify(hits["hits_pos"], hits["hits_neg"],
                                  pos_sub, neg_sub, *GENUS_TH, "genus",
                                  catalog_taxids, nodes, names),
            }
        for level in ("species", "genus"):
            pid, qcov = SPECIES_TH if level == "species" else GENUS_TH
            entry["overall"][level] = classify(
                hits["hits_pos"], hits["hits_neg"], pos, neg, pid, qcov,
                level, catalog_taxids, nodes, names)
        results["databases"][db_name] = entry

    os.makedirs(RESULT_DIR, exist_ok=True)
    out = os.path.join(RESULT_DIR, "ncbi_its_comparison.json")
    json.dump(results, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("wrote", out)
    for db_name, entry in results["databases"].items():
        for level, m in entry["overall"].items():
            print(db_name, level, {k: m[k] for k in
                                   ("sensitivity", "specificity", "f1", "balanced_accuracy")})


if __name__ == "__main__":
    main()
