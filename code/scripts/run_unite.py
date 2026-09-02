#!/usr/bin/env python3
"""Classify fungal queries with the fixed UNITE QIIME2 classifier.

UNITE dynamic fungi is fungi-only, so oomycetes are excluded from the metric
calculation and reported separately as not applicable. Confidence threshold is
0.7, matching the project protocol and the 260-sample UNITE usage.
"""

import csv
import json
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_validation as rv

ROOT = rv.ROOT
QUERY_DIR = rv.QUERY_DIR
RESULT_DIR = os.environ.get(
    "RESULT_DIR",
    os.path.join(ROOT, "docs", "validation", "results"),
)
QIIME = os.environ.get(
    "QIIME",
    os.path.expanduser("~/miniconda3/envs/rachis-qiime2-2026.7/bin/qiime"),
)
UNITE_QZA = os.path.join(
    ROOT, "260samples_fungi", "analysis", "db",
    "unite_ver2025-02-19_dynamic_fungi-Q2-2026.4.qza",
)
CONF_THRESHOLD = 0.7
WORK = os.path.join(RESULT_DIR, "unite_work")
OUT_TAXONOMY = os.path.join(RESULT_DIR, "unite_taxonomy.tsv")
OUT_QUERIES = os.path.join(RESULT_DIR, "unite_per_query.tsv")
OUT_JSON = os.path.join(RESULT_DIR, "unite_comparison.json")
OUT_MD = os.path.join(RESULT_DIR, "UNITE_COMPARISON.md")


def run(cmd, timeout=7200):
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        print(result.stdout[-4000:], file=sys.stderr)
        print(result.stderr[-4000:], file=sys.stderr)
        raise SystemExit(f"command failed: {' '.join(cmd)}")
    return time.time() - t0


def build_query_files():
    pos_file = os.path.join(QUERY_DIR, "external_positives.fasta")
    neg_file = os.path.join(QUERY_DIR, "negatives.fasta")
    pos_records = {rv.parse_query_header(r["header"])["qid"]: r["seq"]
                   for r in rv.read_fasta(pos_file)}
    neg_records = {rv.parse_query_header(r["header"])["qid"]: r["seq"]
                   for r in rv.read_fasta(neg_file)}
    rows = []
    for label, fasta, meta_path in (
        ("pos", pos_file, os.path.join(QUERY_DIR, "external_positives_meta.tsv")),
        ("neg", neg_file, os.path.join(QUERY_DIR, "negatives_meta.tsv")),
    ):
        meta = {}
        with open(meta_path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = row["qid"].split("|", 1)[0]
                meta[key] = row
        for record in rv.read_fasta(fasta):
            info = rv.parse_query_header(record["header"])
            if info["category"] != "fungi":
                continue
            m = meta[info["qid"]]
            rows.append({
                "qid": info["qid"],
                "set": label,
                "taxid": info["taxid"],
                "species": m.get("species", info["species"]),
                "category": info["category"],
                "accession": m.get("accession", info.get("accession", "")),
            })
    with open(os.path.join(WORK, "queries.fasta"), "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(f">{row['qid']}\n")
            seq = pos_records[row["qid"]] if row["set"] == "pos" else neg_records[row["qid"]]
            fh.write(seq + "\n")
    with open(os.path.join(WORK, "queries_meta.tsv"), "w", encoding="utf-8") as fh:
        cols = ["qid", "set", "taxid", "species", "category", "accession"]
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    return rows


def load_unite_results():
    taxonomy = {}
    path = os.path.join(WORK, "taxonomy_export", "taxonomy.tsv")
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Feature ID"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            taxonomy[parts[0]] = {
                "taxonomy": parts[1],
                "confidence": float(parts[2]),
            }
    return taxonomy


def parse_unite_taxonomy(tax_string):
    out = {}
    for field in tax_string.split(";"):
        field = field.strip()
        if len(field) >= 3 and field[1:3] == "__":
            out[field[0]] = field[3:].replace("_", " ").strip()
    return out


def norm(value):
    return rv.norm_taxon(value)


def load_names_by_taxid():
    names_by_taxid = {}
    names_path = os.path.join(rv.TAXDUMP, "names.dmp")
    with open(names_path, encoding="utf-8") as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) >= 4 and fields[0].isdigit():
                taxid = rv.current_taxid(fields[0])
                if fields[3] in ("scientific name", "synonym", "genbank common name"):
                    names_by_taxid.setdefault(taxid, set()).add(norm(fields[1]))
    return names_by_taxid


def load_catalog_name_maps(names_by_taxid, nodes, names):
    """Map normalized NCBI names/synonyms to catalog TaxIDs and genera."""
    species_to_taxid = {}
    genus_set = set()
    with open(os.path.join(ROOT, "data", "species_list.tsv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            taxid = rv.current_taxid(str(row["taxid"]))
            names_to_add = list(names_by_taxid.get(taxid, set())) + [
                row["species"], row.get("full_name", ""), row.get("synonyms", "")
            ]
            for value in names_to_add:
                for name in value.split(";"):
                    species_to_taxid[norm(name)] = taxid
            genus = rv.genus_of(taxid, nodes, names)
            if genus:
                genus_set.add(norm(genus))
    return species_to_taxid, genus_set


def query_species_names(taxid, meta_species, names_by_taxid):
    """All normalized names that should count as a species match."""
    values = {norm(meta_species)}
    values.update(names_by_taxid.get(taxid, set()))
    return values


def evaluate(rows, taxonomy, catalog_map, catalog_genus_set, names_by_taxid,
             nodes, names):
    per_query = []
    stats = {"tp": 0, "fp": 0, "fn": 0, "tn": 0,
             "tp_genus": 0, "fp_genus": 0, "fn_genus": 0, "tn_genus": 0}
    for row in rows:
        call = taxonomy.get(row["qid"])
        taxid = rv.current_taxid(str(row["taxid"]))
        query_species = query_species_names(taxid, row["species"], names_by_taxid)
        query_genus = norm(rv.genus_of(taxid, nodes, names))
        is_pos = row["set"] == "pos"
        if not call or call["confidence"] < CONF_THRESHOLD:
            parsed = {}
        else:
            parsed = parse_unite_taxonomy(call["taxonomy"])
        unite_species = norm(parsed.get("s", ""))
        unite_genus = norm(parsed.get("g", ""))

        species_ok = bool(unite_species) and unite_species in query_species
        genus_ok = bool(unite_genus) and unite_genus == query_genus
        catalog_species_hit = bool(unite_species) and unite_species in catalog_map
        catalog_genus_hit = bool(unite_genus) and unite_genus in catalog_genus_set

        if is_pos:
            stats["tp" if species_ok else "fn"] += 1
            stats["tp_genus" if genus_ok else "fn_genus"] += 1
        else:
            stats["tn" if not catalog_species_hit else "fp"] += 1
            stats["tn_genus" if not catalog_genus_hit else "fp_genus"] += 1

        per_query.append({
            "qid": row["qid"],
            "set": row["set"],
            "taxid": row["taxid"],
            "species": row["species"],
            "accession": row["accession"],
            "unite_taxonomy": call["taxonomy"] if call else "",
            "confidence": round(call["confidence"], 6) if call else 0.0,
            "confidence_ok": bool(call and call["confidence"] >= CONF_THRESHOLD),
            "unite_genus": parsed.get("g", ""),
            "unite_species": parsed.get("s", ""),
            "species_ok": species_ok,
            "genus_ok": genus_ok,
            "catalog_species_hit": catalog_species_hit,
            "catalog_genus_hit": catalog_genus_hit,
        })
    overall_species = rv.metrics(stats["tp"], stats["fp"], stats["fn"], stats["tn"])
    overall_genus = rv.metrics(
        stats["tp_genus"], stats["fp_genus"], stats["fn_genus"], stats["tn_genus"]
    )
    return per_query, {"species": overall_species, "genus": overall_genus}


def write_md(meta, metrics):
    import json as _json
    ncbi_path = os.path.join(RESULT_DIR, "ncbi_its_comparison.json")
    compare = None
    if os.path.exists(ncbi_path):
        with open(ncbi_path, encoding="utf-8") as fh:
            data = _json.load(fh)
        compare = data.get("databases", {})
    lines = [
        "# Classification vs UNITE (QIIME2)",
        "",
        f"**Data**: release; positive/negative fungal queries only (UNITE dynamic fungi has "
        f"no oomycete classes).",
        f"**Classifier**: UNITE v10.0 dynamic fungi, release 2025-02-19, "
        f"`unite_ver2025-02-19_dynamic_fungi-Q2-2026.4.qza`.",
        f"**QIIME2**: {meta['qiime_version']}.",
        f"**Confidence threshold**: >= {CONF_THRESHOLD}; calls below threshold "
        f"are treated as no call.",
        f"**Queries**: {meta['n_queries']} (positive {meta['n_pos']}, negative "
        f"{meta['n_neg']}).",
        "",
        "| Level | TP | FP | FN | TN | Sensitivity | Specificity | F1 | Balanced |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for level in ("species", "genus"):
        m = metrics[level]
        lines.append(
            f"| {level} | {m['tp']} | {m['fp']} | {m['fn']} | {m['tn']} | "
            f"{m['sensitivity']}% | {m['specificity']}% | {m['f1']} | "
            f"{m['balanced_accuracy']}% |"
        )
    lines += [
        "",
        "## Notes",
        "",
        "- Oomycetes are excluded because the UNITE dynamic fungi classifier "
        "contains no oomycete reference classes; forcing oomycete queries "
        "through it would report them as Fungi and produce misleading metrics.",
        "- Species matching uses NCBI taxonomy names and synonyms for each "
        "query TaxID; UNITE species names that are not present in NCBI "
        "taxonomy are treated as mismatches and should be reviewed manually.",
        "- Full per-query calls are not included in the public package.",
        "",
    ]
    if compare:
        rows = []
        for db, label in (
            ("MPDB", "MaizePathogenDB (BLAST 99/90, 95/70)"),
            ("NCBI_ITS_eukaryote", "NCBI ITS_eukaryote (BLAST 99/90, 95/70)"),
        ):
            fungi = compare.get(db, {}).get("per_category", {}).get("fungi", {})
            if not fungi:
                continue
            sp = fungi.get("species", {})
            ge = fungi.get("genus", {})
            rows.append(
                f"| {label} | species | {sp.get('sensitivity')}% | "
                f"{sp.get('specificity')}% | {sp.get('f1')} | "
                f"{sp.get('balanced_accuracy')}% |"
            )
            rows.append(
                f"| {label} | genus | {ge.get('sensitivity')}% | "
                f"{ge.get('specificity')}% | {ge.get('f1')} | "
                f"{ge.get('balanced_accuracy')}% |"
            )
        if rows:
            lines += [
                "## Comparison on fungal queries only (same query set)",
                "",
                "UNITE is a naive-Bayes classifier with confidence >= 0.7; BLAST "
                "databases use species 99/90 and genus 95/70. Metrics are not "
                "threshold-equivalent.",
                "",
                "| Method | Level | Sensitivity | Specificity | F1 | Balanced |",
                "|---|---|---:|---:|---:|---:|",
            ]
            lines += rows
            lines += ["", "Source for MPDB/NCBI rows: `ncbi_its_comparison.json`."]
    with open(OUT_MD, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def main():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(os.path.join(WORK, "taxonomy_export"), exist_ok=True)

    rows = build_query_files()
    t_classify = 0.0
    if os.environ.get("SKIP_QIIME") != "1":
        run([QIIME, "tools", "import", "--type", "FeatureData[Sequence]",
             "--input-path", os.path.join(WORK, "queries.fasta"),
             "--output-path", os.path.join(WORK, "queries.qza")])
        t_classify = run([QIIME, "feature-classifier", "classify-sklearn",
                          "--i-classifier", UNITE_QZA,
                          "--i-reads", os.path.join(WORK, "queries.qza"),
                          "--o-classification", os.path.join(WORK, "taxonomy.qza"),
                          "--p-n-jobs", "4"])
        run([QIIME, "tools", "export",
             "--input-path", os.path.join(WORK, "taxonomy.qza"),
             "--output-path", os.path.join(WORK, "taxonomy_export")])
        with open(os.path.join(WORK, "classify_seconds.txt"), "w") as fh:
            fh.write(str(round(t_classify, 1)))
    elif os.path.exists(os.path.join(WORK, "classify_seconds.txt")):
        with open(os.path.join(WORK, "classify_seconds.txt")) as fh:
            t_classify = float(fh.read().strip())

    taxonomy = load_unite_results()
    nodes, names = rv.load_taxdump()
    names_by_taxid = load_names_by_taxid()
    catalog_map, catalog_genus_set = load_catalog_name_maps(
        names_by_taxid, nodes, names
    )
    per_query, metrics = evaluate(
        rows, taxonomy, catalog_map, catalog_genus_set, names_by_taxid,
        nodes, names
    )

    qiime_info = subprocess.run([QIIME, "info"], capture_output=True, text=True)
    qiime_version = ""
    for line in (qiime_info.stdout or "").splitlines():
        if "rachis version" in line:
            qiime_version = "QIIME2 2026.7 (rachis " + line.split(":", 1)[1].strip() + ")"
    meta = {
        "validation": "unite_comparison",
        "data_version": "release",
        "classifier": "UNITE v10.0 dynamic fungi 2025-02-19 (Q2-2026.4 qza)",
        "qiime_version": qiime_version or "unknown",
        "confidence_threshold": CONF_THRESHOLD,
        "n_queries": len(rows),
        "n_pos": sum(1 for r in rows if r["set"] == "pos"),
        "n_neg": sum(1 for r in rows if r["set"] == "neg"),
        "classify_s": round(t_classify, 1),
    }

    with open(OUT_TAXONOMY, "w", encoding="utf-8") as fh:
        fh.write("Feature ID\tTaxon\tConfidence\n")
        for qid, call in taxonomy.items():
            fh.write(f"{qid}\t{call['taxonomy']}\t{call['confidence']}\n")
    with open(OUT_QUERIES, "w", encoding="utf-8") as fh:
        cols = list(per_query[0].keys())
        fh.write("\t".join(cols) + "\n")
        for row in per_query:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump({"meta": meta, "overall": metrics}, fh, ensure_ascii=False,
                  indent=2)
    write_md(meta, metrics)
    print(json.dumps({"meta": meta, "overall": metrics}, indent=2))


if __name__ == "__main__":
    main()
