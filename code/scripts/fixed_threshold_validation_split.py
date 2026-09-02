#!/usr/bin/env python3
"""Evaluate fixed recommended thresholds on a TaxID-level train/validation split.

No threshold selection happens here. Species calls use pident>=99, qcovs>=90;
genus calls use pident>=95, qcovs>=70. Both halves are reported so the
validation half is an unbiased check of the fixed thresholds.
"""

import csv
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_validation as rv

ROOT = rv.ROOT
QUERY_DIR = rv.QUERY_DIR
BLAST_DIR = os.environ.get("BLAST_DIR", rv.BLAST_DIR)
RESULT_DIR = os.environ.get("RESULT_DIR", rv.RESULT_DIR)
SEED = 42

LEVEL_THRESHOLDS = {
    "species": (99.0, 90.0),
    "genus": (95.0, 70.0),
}


def load_queries(fasta, meta_path):
    queries = []
    meta = {}
    with open(meta_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["qid"]] = row
    for record in rv.read_fasta(fasta):
        info = rv.parse_query_header(record["header"])
        m = meta.get(info["qid"], {})
        queries.append({
            **info,
            "seq": record["seq"],
            "species": m.get("species", info.get("species", "")),
        })
    return queries


def split_by_taxid(items, rng):
    taxids = sorted({item["taxid"] for item in items})
    rng.shuffle(taxids)
    half = (len(taxids) + 1) // 2
    tune_taxids = set(taxids[:half])
    tune = [x for x in items if x["taxid"] in tune_taxids]
    val = [x for x in items if x["taxid"] not in tune_taxids]
    return tune, val


def classify(pos_hits, neg_hits, pos, neg, pid, qcov, level,
             catalog_taxids, nodes, names):
    tp = fp = fn = tn = 0
    for record in pos:
        hit = pos_hits.get(record["qid"])
        if not hit or hit["pident"] < pid or hit["qcovs"] < qcov:
            fn += 1
            continue
        hit_taxid = rv.subject_taxid(hit["sseqid"])
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
        hit_taxid = rv.current_taxid(rv.subject_taxid(hit["sseqid"]))
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
    pos = load_queries(
        os.path.join(QUERY_DIR, "external_positives.fasta"),
        os.path.join(QUERY_DIR, "external_positives_meta.tsv"),
    )
    neg = load_queries(
        os.path.join(QUERY_DIR, "negatives.fasta"),
        os.path.join(QUERY_DIR, "negatives_meta.tsv"),
    )
    db = os.path.join(BLAST_DIR, "maize_pathogens_all")
    pos_hits = rv.run_blast(os.path.join(QUERY_DIR, "external_positives.fasta"), db)
    neg_hits = rv.run_blast(os.path.join(QUERY_DIR, "negatives.fasta"), db)

    rng = random.Random(SEED)
    results = {
        "seed": SEED,
        "split": "by_taxid_stratified_by_category",
        "thresholds": LEVEL_THRESHOLDS,
        "note": "fixed recommended thresholds; no threshold selection",
        "categories": {},
        "overall_validation": {},
    }
    val_sums = {"species": [0, 0, 0, 0], "genus": [0, 0, 0, 0]}
    for category in rv.CATEGORY_ORDER:
        pos_cat = [p for p in pos if p["category"] == category]
        neg_cat = [n for n in neg if n["category"] == category]
        pos_tune, pos_val = split_by_taxid(pos_cat, rng)
        neg_tune, neg_val = split_by_taxid(neg_cat, rng)
        entry = {
            "category": category,
            "n_pos_tune": len(pos_tune), "n_pos_val": len(pos_val),
            "n_neg_tune": len(neg_tune), "n_neg_val": len(neg_val),
            "levels": {},
        }
        for level, (pid, qcov) in LEVEL_THRESHOLDS.items():
            tune_m = classify(pos_hits, neg_hits, pos_tune, neg_tune,
                              pid, qcov, level, catalog_taxids, nodes, names)
            val_m = classify(pos_hits, neg_hits, pos_val, neg_val,
                             pid, qcov, level, catalog_taxids, nodes, names)
            entry["levels"][level] = {
                "thresholds": {"pident": pid, "qcovs": qcov},
                "tuning_metrics": {k: tune_m[k] for k in
                                   ("tp", "fp", "fn", "tn", "sensitivity", "specificity",
                                    "precision", "f1", "balanced_accuracy")},
                "validation_metrics": {k: val_m[k] for k in
                                       ("tp", "fp", "fn", "tn", "sensitivity", "specificity",
                                        "precision", "f1", "balanced_accuracy")},
            }
            for i, key in enumerate(("tp", "fp", "fn", "tn")):
                val_sums[level][i] += val_m[key]
        results["categories"][category] = entry

    for level, sums in val_sums.items():
        tp, fp, fn, tn = sums
        results["overall_validation"][level] = rv.metrics(tp, fp, fn, tn)

    os.makedirs(RESULT_DIR, exist_ok=True)
    out = os.path.join(RESULT_DIR, "fixed_threshold_validation_split.json")
    json.dump(results, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("wrote", out)
    for level, overall in results["overall_validation"].items():
        print(level, "overall validation:",
              {k: overall[k] for k in ("tp", "fp", "fn", "tn", "sensitivity",
                                       "specificity", "f1", "balanced_accuracy")})


if __name__ == "__main__":
    main()
