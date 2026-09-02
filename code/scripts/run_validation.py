#!/usr/bin/env python3
"""Run the MaizePathogenDB validation suite per PROTOCOL.md."""

import csv
import json
import os
import re
import subprocess
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEQ_DIR = os.environ.get("SEQ_DIR", os.path.join(ROOT, "release", "sequences"))
BLAST_DIR = os.environ.get("BLAST_DIR", os.path.join(ROOT, "release", "blast_db"))
QUERY_DIR = os.path.join(ROOT, "docs", "validation", "query_sets")
RESULT_DIR = os.environ.get(
    "RESULT_DIR", os.path.join(ROOT, "docs", "validation", "results")
)
DB_VERSION = os.environ.get("DB_VERSION", "release")
TAXONOMY_JSON = os.path.join(ROOT, "Figshare", "taxonomy.json")
TAXDUMP = os.path.join(ROOT, "260samples_fungi", "analysis", "db", "taxdump")
MERGED_DMP = os.path.join(TAXDUMP, "merged.dmp")

BLASTN = os.environ.get("BLASTN", "blastn")
CATEGORY_ORDER = ["bacteria", "viruses", "fungi", "oomycetes"]

PRIMERS = {
    "27F": "AGAGTTTGATCMTGGCTCAG",
    "1492R": "GGTTACCTTGTTACGACTT",
    "338F": "ACTCCTACGGGAGGCAGCAG",
    "806R": "GGACTACHVGGGTWTCTAAT",
    "515F": "GTGCCAGCMGCCGCGGTAA",
    "1392R": "ACGGGCGGTGTGTRC",
    "926R": "CCGYCAATTYMTTTRAGTTT",
    "ITS1F": "CTTGGTCATTTAGAGGAAGTAA",
    "ITS4": "TCCTCCGCTTATTGATATGC",
    "ITS5": "GGAAGTAAAAGTCGTAACAAGG",
    "ITS1": "TCCGTAGGTGAACCTGCGG",
    "ITS86F": "GTGAATCATCGAATCTTTGAA",
    "fITS7": "GTGARTCATCGAATCTTTG",
    "ITS2": "GCTGCGTTCTTCATCGATGC",
    "ITS3": "GCATCGATGAAGAACGCAGC",
    "ITS9mun": "GTACACACCGCCCGTCG",
    "ITS4ngs": "TCCTSCGCTTATTGATATGC",
}

PRIMER_PAIRS = {
    "16S_V3V4": ("338F", "806R"),
    "16S_V4": ("515F", "806R"),
    "16S_V3V8": ("338F", "1392R"),
    "16S_full_27F_1492R": ("27F", "1492R"),
    "16S_V4V5": ("515F", "926R"),
    "ITS_ITS1F_ITS4": ("ITS1F", "ITS4"),
    "ITS_ITS5_ITS4": ("ITS5", "ITS4"),
    "ITS_ITS1_ITS4": ("ITS1", "ITS4"),
    "ITS_ITS86F_ITS4": ("ITS86F", "ITS4"),
    "ITS_fITS7_ITS4": ("fITS7", "ITS4"),
    "ITS_ITS1F_ITS2": ("ITS1F", "ITS2"),
    "ITS_ITS3_ITS4": ("ITS3", "ITS4"),
    "ITS_ITS9mun_ITS4ngs": ("ITS9mun", "ITS4ngs"),
}

AMBIG = {
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
MERGED = {}


def load_merged():
    global MERGED
    if MERGED:
        return MERGED
    if os.path.exists(MERGED_DMP):
        with open(MERGED_DMP, encoding="utf-8") as fh:
            for line in fh:
                fields = [f.strip() for f in line.split("|")]
                if len(fields) >= 2 and fields[0].isdigit() and fields[1].isdigit():
                    MERGED[fields[0]] = fields[1]
    return MERGED


def current_taxid(taxid):
    taxid = str(taxid)
    seen = set()
    while taxid in MERGED and taxid not in seen:
        seen.add(taxid)
        taxid = MERGED[taxid]
    return taxid


def allowed_bases(base):
    return set(AMBIG.get(base.upper(), base.upper()))


def primer_site(primer, seq, max_mismatch=2):
    p = primer.upper()
    seq = seq.upper()
    for start in range(0, len(seq) - len(p) + 1):
        mismatches = 0
        for i in range(len(p)):
            if seq[start + i] not in allowed_bases(p[i]):
                mismatches += 1
                if mismatches > max_mismatch:
                    break
        if mismatches <= max_mismatch:
            return True
    return False


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_taxonomy():
    with open(TAXONOMY_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def load_taxdump():
    nodes = {}
    names = {}
    nodes_path = os.path.join(TAXDUMP, "nodes.dmp")
    names_path = os.path.join(TAXDUMP, "names.dmp")
    with open(nodes_path, encoding="utf-8") as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) >= 3:
                nodes[fields[0]] = (fields[1], fields[2])
    with open(names_path, encoding="utf-8") as fh:
        for line in fh:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) >= 4 and fields[3] == "scientific name":
                names[fields[0]] = fields[1]
    return nodes, names


def genus_of(taxid, nodes, names):
    taxid = current_taxid(taxid)
    seen = set()
    current = str(taxid)
    while current and current not in seen:
        seen.add(current)
        node = nodes.get(current)
        if not node:
            return names.get(current, "")
        parent, rank = node
        if rank == "genus" or rank == "subgenus":
            return names.get(current, "")
        current = parent
    return ""


def norm_taxon(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def read_fasta(path):
    records = []
    header = None
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    records.append({"header": header, "seq": "".join(seq)})
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            records.append({"header": header, "seq": "".join(seq)})
    return records


def parse_query_header(header):
    parts = header.split("|")
    return {
        "qid": parts[0], "taxid": parts[1], "species": parts[2],
        "category": parts[3], "accession": parts[4] if len(parts) > 4 else "",
    }


def run_blast(query_file, db, max_targets=1):
    cmd = [
        BLASTN, "-query", query_file, "-db", db,
        "-outfmt", "6 qseqid sseqid pident qcovs",
        "-max_target_seqs", str(max_targets), "-evalue", "1e-5", "-num_threads", "4",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    hits = {}
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        qid = parts[0].split("|", 1)[0]
        if qid not in hits:
            hits[qid] = {
                "sseqid": parts[1], "pident": float(parts[2]), "qcovs": float(parts[3]),
            }
    return hits


def subject_taxid(sseqid):
    return sseqid.split("|")[1]


def accuracy_rows(query_records, hits, nodes, names, thresholds=None):
    rows = []
    for record in query_records:
        meta = parse_query_header(record["header"])
        hit = hits.get(meta["qid"])
        hit_taxid = subject_taxid(hit["sseqid"]) if hit else ""
        species_ok = bool(hit) and current_taxid(hit_taxid) == current_taxid(meta["taxid"])
        hit_genus = genus_of(hit_taxid, nodes, names) if hit else ""
        query_genus = genus_of(meta["taxid"], nodes, names)
        genus_ok = bool(hit) and bool(hit_genus) and (
            norm_taxon(hit_genus) == norm_taxon(query_genus)
        )
        row = {
            "qid": meta["qid"], "taxid": meta["taxid"], "species": meta["species"],
            "category": meta["category"], "hit_taxid": hit_taxid,
            "pident": hit["pident"] if hit else 0.0,
            "qcovs": hit["qcovs"] if hit else 0.0,
            "species_retrieval_ok": species_ok,
            "genus_retrieval_ok": genus_ok,
        }
        if thresholds:
            pid, qcov = thresholds["species"]
            row["species_class_ok"] = species_ok and hit["pident"] >= pid and hit["qcovs"] >= qcov
            pid_g, qcov_g = thresholds["genus"]
            row["genus_class_ok"] = genus_ok and hit["pident"] >= pid_g and hit["qcovs"] >= qcov_g
        rows.append(row)
    return rows


def summarize_accuracy(rows):
    summary = {}
    for category in CATEGORY_ORDER:
        subset = [r for r in rows if r["category"] == category]
        if not subset:
            continue
        summary[category] = {}
        for key in [k for k in subset[0] if k.endswith("_ok")]:
            correct = sum(1 for r in subset if r[key])
            summary[category][key] = {
                "n": len(subset), "correct": correct,
                "accuracy": round(correct / len(subset) * 100, 1),
            }
    return summary


def write_json(path, data):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=2)


def internal_completeness():
    query_file = os.path.join(SEQ_DIR, "maize_pathogens_all.fasta")
    queries = read_fasta(query_file)
    hits = run_blast(query_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    per_category = defaultdict(list)
    for record in queries:
        meta = parse_query_header(record["header"])
        hit = hits.get(meta["qid"])
        ok = bool(hit) and current_taxid(subject_taxid(hit["sseqid"])) == current_taxid(meta["taxid"])
        per_category[meta["category"]].append({
            "qid": meta["qid"], "taxid": meta["taxid"], "hit": bool(hit), "ok": ok,
        })
    results = {}
    for category in CATEGORY_ORDER:
        details = per_category.get(category, [])
        correct = sum(1 for d in details if d["ok"])
        results[category] = {
            "n": len(details), "correct": correct,
            "accuracy": round(correct / len(details) * 100, 1) if details else None,
            "details": details,
        }
    return {"validation": "internal_completeness", "results": results}


def primer_coverage():
    sequences = defaultdict(list)
    for record in read_fasta(os.path.join(SEQ_DIR, "maize_pathogens_all.fasta")):
        meta = parse_query_header(record["header"])
        sequences[meta["category"]].append(record["seq"])
    coverage = {}
    for pair_name, (fwd, rev) in PRIMER_PAIRS.items():
        if pair_name.startswith("16S"):
            seqs = sequences.get("bacteria", [])
            category_key = "bacteria"
        else:
            seqs = sequences.get("fungi", []) + sequences.get("oomycetes", [])
            category_key = "fungi+oomycetes"
        reverse_primers = {rev} if rev.endswith("R") or rev in ("ITS2", "ITS4", "ITS4ngs") else set()
        covered = 0
        for seq in seqs:
            fwd_ok = primer_site(PRIMERS[fwd], seq)
            rev_seq = revcomp(seq) if rev in reverse_primers else seq
            rev_ok = primer_site(PRIMERS[rev], rev_seq)
            covered += int(fwd_ok and rev_ok)
        coverage[pair_name] = {
            "category": category_key, "n": len(seqs), "covered": covered,
            "coverage_pct": round(covered / len(seqs) * 100, 1) if seqs else None,
        }
    return {"validation": "primer_coverage", "results": coverage}


def external_retrieval():
    nodes, names = load_taxdump()
    query_file = os.path.join(QUERY_DIR, "external_positives.fasta")
    queries = read_fasta(query_file)
    hits = run_blast(query_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    thresholds = {"species": (99.0, 90.0), "genus": (95.0, 70.0)}
    rows = accuracy_rows(queries, hits, nodes, names, thresholds)
    return {
        "validation": "external_retrieval_mpdb",
        "thresholds": thresholds,
        "summary": summarize_accuracy(rows),
        "details": rows,
    }


def metrics(tp, fp, fn, tn):
    sensitivity = tp / (tp + fn) if (tp + fn) else None
    specificity = tn / (tn + fp) if (tn + fp) else None
    precision = tp / (tp + fp) if (tp + fp) else None
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision and sensitivity else None
    balanced = (sensitivity + specificity) / 2 if sensitivity is not None and specificity is not None else None
    return {
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "sensitivity": round(sensitivity * 100, 1) if sensitivity is not None else None,
        "specificity": round(specificity * 100, 1) if specificity is not None else None,
        "precision": round(precision * 100, 1) if precision is not None else None,
        "f1": round(f1 * 100, 1) if f1 is not None else None,
        "balanced_accuracy": round(balanced * 100, 1) if balanced is not None else None,
    }


def classification_benchmark():
    nodes, names = load_taxdump()
    catalog_taxids = {current_taxid(str(r.get("taxid"))) for r in load_taxonomy()}
    pos_file = os.path.join(QUERY_DIR, "external_positives.fasta")
    neg_file = os.path.join(QUERY_DIR, "negatives.fasta")
    pos_queries = read_fasta(pos_file)
    neg_queries = read_fasta(neg_file)
    pos_hits = run_blast(pos_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    neg_hits = run_blast(neg_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    pid, qcov = 99.0, 90.0
    per_category = {}
    for category in CATEGORY_ORDER:
        tp = fp = fn = tn = 0
        for record in pos_queries:
            meta = parse_query_header(record["header"])
            if meta["category"] != category:
                continue
            hit = pos_hits.get(meta["qid"])
            if (hit and hit["pident"] >= pid and hit["qcovs"] >= qcov
                    and current_taxid(subject_taxid(hit["sseqid"])) == current_taxid(meta["taxid"])):
                tp += 1
            else:
                fn += 1
        for record in neg_queries:
            meta = parse_query_header(record["header"])
            if meta["category"] != category:
                continue
            hit = neg_hits.get(meta["qid"])
            if (hit and hit["pident"] >= pid and hit["qcovs"] >= qcov
                    and current_taxid(subject_taxid(hit["sseqid"])) in catalog_taxids):
                fp += 1
            else:
                tn += 1
        per_category[category] = metrics(tp, fp, fn, tn)
    overall_tp = sum(v["tp"] for v in per_category.values())
    overall_fp = sum(v["fp"] for v in per_category.values())
    overall_fn = sum(v["fn"] for v in per_category.values())
    overall_tn = sum(v["tn"] for v in per_category.values())
    return {
        "validation": "classification_benchmark_mpdb",
        "thresholds": {"species_pident": pid, "species_qcov": qcov},
        "per_category": per_category,
        "overall": metrics(overall_tp, overall_fp, overall_fn, overall_tn),
    }


def cross_database_consistency():
    nodes, names = load_taxdump()
    query_file = os.path.join(QUERY_DIR, "cross_db_queries.fasta")
    queries = read_fasta(query_file)
    hits = run_blast(query_file, os.path.join(BLAST_DIR, "maize_pathogens_all"))
    rows = accuracy_rows(queries, hits, nodes, names)
    return {
        "validation": "cross_database_consistency",
        "summary": summarize_accuracy(rows),
        "details": rows,
    }


def write_summary(completeness, coverage, retrieval, benchmark, cross_db):
    def simple_acc_table(summary):
        lines = ["| Category | n | Correct | Accuracy |", "|---|---:|---:|---:|"]
        for category in CATEGORY_ORDER:
            row = summary.get(category)
            if not row or row.get("n") is None:
                continue
            lines.append(f"| {category} | {row['n']} | {row['correct']} | {row['accuracy']} |")
        return "\n".join(lines)

    def acc_table(summary, keys):
        lines = ["| Category | n | Correct | Accuracy |", "|---|---:|---:|---:|"]
        for category in CATEGORY_ORDER:
            row = summary.get(category)
            if not row:
                continue
            value = row.get(keys)
            if not value:
                continue
            lines.append(f"| {category} | {value['n']} | {value['correct']} | {value['accuracy']} |")
        return "\n".join(lines)

    lines = [
        f"# MaizePathogenDB {DB_VERSION} Validation Results",
        "",
        "- Protocol: `docs/validation/PROTOCOL.md`",
        "- Run: 2026-08-26",
        "- Independent positive queries: 675; negative queries: 500; cross-database queries: 565",
        "",
        "## Internal completeness",
        "",
        simple_acc_table(completeness["results"]),
        "",
        "## Primer / region coverage",
        "",
        "| Primer pair | Category | n | Covered | Coverage % |",
        "|---|---:|---:|---:|---:|",
    ]
    for pair, value in coverage["results"].items():
        lines.append(
            f"| {pair} | {value['category']} | {value['n']} | {value['covered']} | {value['coverage_pct']} |"
        )
    lines += [
        "",
        "## External retrieval against MPDB",
        "",
        "### Species-level retrieval (top-1, no threshold)",
        "",
        acc_table(retrieval["summary"], "species_retrieval_ok"),
        "",
        "### Genus-level retrieval",
        "",
        acc_table(retrieval["summary"], "genus_retrieval_ok"),
        "",
        "### Species-level classification (pident>=99, qcovs>=90)",
        "",
        acc_table(retrieval["summary"], "species_class_ok"),
        "",
        "### Genus-level classification (pident>=95, qcovs>=70)",
        "",
        acc_table(retrieval["summary"], "genus_class_ok"),
        "",
        "## Classification benchmark (species 99/90)",
        "",
        "| Category | n_pos | n_neg | TP | FP | FN | TN | Sensitivity | Specificity | Precision | F1 | Balanced |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for category in CATEGORY_ORDER:
        row = benchmark["per_category"][category]
        lines.append(
            f"| {category} | {row['tp'] + row['fn']} | {row['fp'] + row['tn']} | {row['tp']} | {row['fp']} | "
            f"{row['fn']} | {row['tn']} | {row['sensitivity']} | {row['specificity']} | "
            f"{row['precision']} | {row['f1']} | {row['balanced_accuracy']} |"
        )
    row = benchmark["overall"]
    lines.append(
        f"| Overall | {row['tp'] + row['fn']} | {row['fp'] + row['tn']} | {row['tp']} | {row['fp']} | "
        f"{row['fn']} | {row['tn']} | {row['sensitivity']} | {row['specificity']} | "
        f"{row['precision']} | {row['f1']} | {row['balanced_accuracy']} |"
    )
    lines += [
        "",
        "## Cross-database consistency",
        "",
        "### Species-level retrieval",
        "",
        acc_table(cross_db["summary"], "species_retrieval_ok"),
        "",
        "### Genus-level retrieval",
        "",
        acc_table(cross_db["summary"], "genus_retrieval_ok"),
        "",
        "## Additional comparisons",
        "",
        "- NCBI-nt head-to-head: `NCBI_NT_COMPARISON.md`.",
        "- NCBI ITS_eukaryote / ITS_RefSeq comparison: `NCBI_ITS_COMPARISON.md` and `ncbi_its_comparison.json`.",
        "- UNITE comparison: `UNITE_COMPARISON.md` and `unite_comparison.json`.",
        "- Performance (Usage Notes only): `PERFORMANCE.md` and `performance.json`.",
    ]
    with open(os.path.join(RESULT_DIR, "SUMMARY.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    load_merged()
    completeness = internal_completeness()
    coverage = primer_coverage()
    retrieval = external_retrieval()
    benchmark = classification_benchmark()
    cross_db = cross_database_consistency()
    write_json(os.path.join(RESULT_DIR, "internal_completeness.json"), completeness)
    write_json(os.path.join(RESULT_DIR, "primer_coverage.json"), coverage)
    write_json(os.path.join(RESULT_DIR, "external_retrieval_mpdb.json"), retrieval)
    write_json(os.path.join(RESULT_DIR, "classification_benchmark.json"), benchmark)
    write_json(os.path.join(RESULT_DIR, "cross_database_consistency.json"), cross_db)
    results = {
        "internal_completeness": completeness,
        "primer_coverage": coverage,
        "external_retrieval": retrieval,
        "classification_benchmark": benchmark,
        "cross_database_consistency": cross_db,
    }
    write_summary(completeness, coverage, retrieval, benchmark, cross_db)
    print("internal completeness:", {k: (v["accuracy"], v["n"]) for k, v in completeness["results"].items()})
    print("primer coverage pairs:", {k: v["coverage_pct"] for k, v in coverage["results"].items()})
    print("external retrieval:", {k: v.get("species_retrieval_ok", {}) for k, v in retrieval["summary"].items()})
    print("classification benchmark overall:", benchmark["overall"])
    print("cross-database consistency:", {k: v.get("species_retrieval_ok", {}) for k, v in cross_db["summary"].items()})
    print("results written to", RESULT_DIR)


if __name__ == "__main__":
    main()
