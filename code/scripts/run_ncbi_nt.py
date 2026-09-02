#!/usr/bin/env python3
"""NCBI-nt head-to-head against a fixed local snapshot.

Independent positive and negative queries are searched against NCBI-nt with the same thresholds used for MPDB
(species: pident>=99, qcovs>=90; genus: pident>=95, qcovs>=70).  Searches
are restricted to each query category's NCBI lineage via -taxidlist so that a
fungal query is not searched against the entire nt database; the restriction
is applied identically to positives and negatives and is recorded in the
output metadata.
"""

import argparse
import concurrent.futures
import csv
import json
import os
import subprocess
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QUERY_DIR = os.path.join(ROOT, "docs", "validation", "query_sets")
RESULT_DIR = os.path.join(ROOT, "docs", "validation", "results")
WORK_DIR = os.path.join(RESULT_DIR, "ncbi_nt")
BLASTN = os.environ.get("BLASTN", "blastn")
TAXDUMP = os.path.join(ROOT, "260samples_fungi", "analysis", "db", "taxdump")
TAXONOMY_JSON = os.path.join(ROOT, "Figshare", "taxonomy.json")
NT_DB = os.environ.get("NT_DB", os.path.join(ROOT, "nt_snapshot", "nt"))
BLASTDB = os.environ.get("BLASTDB", os.path.dirname(NT_DB))

CATEGORY_ORDER = ["bacteria", "viruses", "fungi", "oomycetes"]
THRESHOLDS = {"species": (99.0, 90.0), "genus": (95.0, 70.0)}


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


def parse_header(header):
    parts = header.split("|")
    return {
        "qid": parts[0],
        "taxid": parts[1],
        "species": parts[2] if len(parts) > 2 else "",
        "category": parts[3] if len(parts) > 3 else "",
    }


def load_taxdump():
    nodes = {}
    names = {}
    with open(os.path.join(TAXDUMP, "nodes.dmp"), encoding="utf-8") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 3:
                nodes[f[0]] = (f[1], f[2])
    with open(os.path.join(TAXDUMP, "names.dmp"), encoding="utf-8") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 4 and f[3] == "scientific name":
                names[f[0]] = f[1]
    merged = {}
    with open(os.path.join(TAXDUMP, "merged.dmp"), encoding="utf-8") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 2 and f[0].isdigit() and f[1].isdigit():
                merged[f[0]] = f[1]
    return nodes, names, merged


def current_taxid(taxid, merged):
    taxid = str(taxid)
    seen = set()
    while taxid in merged and taxid not in seen:
        seen.add(taxid)
        taxid = merged[taxid]
    return taxid


def lineage(taxid, nodes, merged):
    """Return taxid chain from node to root (current IDs)."""
    taxid = current_taxid(taxid, merged)
    chain = []
    seen = set()
    cur = taxid
    while cur and cur not in seen:
        seen.add(cur)
        chain.append(cur)
        node = nodes.get(cur)
        if not node:
            break
        cur = node[0]
    return chain


def rank_taxid(taxid, nodes, merged, ranks=("species", "genus", "subgenus")):
    chain = lineage(taxid, nodes, merged)
    # NCBI sometimes stores a virus taxid as 'no rank' with its parent
    # carrying the 'species' rank (e.g. Cucumovirus CMV).  Treat the parent as
    # the species-level node in that case.
    if ranks and ranks[0] == "species" and chain:
        node = nodes.get(chain[0])
        if node and node[1] not in ranks and len(chain) > 1:
            parent = nodes.get(chain[1])
            if parent and parent[1] in ranks:
                return chain[1]
    for cur in chain:
        node = nodes.get(cur)
        if node and node[1] in ranks:
            return cur
    # Fall back to the first node when no explicit rank was found.
    if ranks and ranks[0] == "species" and chain:
        return chain[0]
    return None


def descendants(root, children):
    stack = [root]
    out = set()
    while stack:
        cur = stack.pop()
        if cur in out:
            continue
        out.add(cur)
        stack.extend(children.get(cur, ()))
    return out


def load_catalog(nodes, merged):
    data = json.load(open(TAXONOMY_JSON, encoding="utf-8"))
    catalog = []
    for row in data:
        taxid = current_taxid(row["taxid"], merged)
        genus = rank_taxid(taxid, nodes, merged, ("genus", "subgenus"))
        catalog.append({"taxid": taxid, "genus": genus, "species": row["species"]})
    return catalog


def build_taxidlists(query_records, nodes, merged):
    children = defaultdict(list)
    for child, (parent, _rank) in nodes.items():
        children[parent].append(child)
    per_cat = defaultdict(set)
    roots = defaultdict(set)
    for rec in query_records:
        meta = parse_header(rec["header"])
        cat = meta["category"]
        per_cat[cat].add(current_taxid(meta["taxid"], merged))
        roots[cat].update(lineage(meta["taxid"], nodes, merged)[-1:])
    out = {}
    for cat in CATEGORY_ORDER:
        if not per_cat.get(cat):
            continue
        taxids = set()
        for root in roots.get(cat, ()):
            taxids.update(descendants(root, children))
        if not taxids:
            taxids = per_cat[cat]
        out[cat] = taxids
    return out


def write_taxidlist(cat, taxids):
    os.makedirs(WORK_DIR, exist_ok=True)
    path = os.path.join(WORK_DIR, f"taxidlist_{cat}.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(sorted(taxids)) + "\n")
    return path


def split_fasta(records, batch_size):
    for i in range(0, len(records), batch_size):
        chunk = records[i : i + batch_size]
        path = os.path.join(WORK_DIR, f"batch_{i // batch_size:03d}.fasta")
        with open(path, "w", encoding="utf-8") as fh:
            for rec in chunk:
                fh.write(f">{rec['header']}\n{rec['seq']}\n")
        yield i // batch_size, path


def blast_batch(batch_id, fasta, cat, taxidlist=None):
    out_path = os.path.join(WORK_DIR, f"hits_{cat}_{batch_id:03d}.tsv")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return batch_id, out_path
    cmd = [
        BLASTN, "-query", fasta, "-db", NT_DB,
        "-outfmt", "6 qseqid sseqid staxids pident qcovs",
        "-max_target_seqs", "5", "-evalue", "1e-5",
        "-num_threads", "4", "-out", out_path,
    ]
    if taxidlist:
        cmd[3:3] = ["-taxidlist", taxidlist]
    env = dict(os.environ, BLASTDB=BLASTDB)
    subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
    return batch_id, out_path


def summarize_group(records, nodes, names, merged, catalog):
    """Evaluate one set of records (all positives or all negatives)."""
    catalog_taxids = {c["taxid"] for c in catalog}
    catalog_desc = set()
    children = defaultdict(list)
    for child, (parent, _rank) in nodes.items():
        children[parent].append(child)
    for taxid in catalog_taxids:
        catalog_desc.update(descendants(taxid, children))
    catalog_desc.update(catalog_taxids)

    rows = []
    for rec in records:
        meta = parse_header(rec["header"])
        hit = rec.get("hit")
        query_taxid = current_taxid(meta["taxid"], merged)
        query_species = rank_taxid(query_taxid, nodes, merged, ("species",))
        query_genus = rank_taxid(query_taxid, nodes, merged, ("genus", "subgenus"))
        hit_taxid = ""
        hit_species = None
        hit_genus = None
        if hit and hit.get("staxid"):
            hit_taxid = current_taxid(hit["staxid"], merged)
            hit_species = rank_taxid(hit_taxid, nodes, merged, ("species",))
            hit_genus = rank_taxid(hit_taxid, nodes, merged, ("genus", "subgenus"))
        species_retrieval = bool(hit) and hit_species is not None and query_species is not None and hit_species == query_species
        genus_retrieval = bool(hit) and hit_genus is not None and query_genus is not None and hit_genus == query_genus
        pid = hit["pident"] if hit else 0.0
        qcov = hit["qcovs"] if hit else 0.0
        spid, sqcov = THRESHOLDS["species"]
        gpid, gqcov = THRESHOLDS["genus"]
        species_class = species_retrieval and pid >= spid and qcov >= sqcov
        genus_class = genus_retrieval and pid >= gpid and qcov >= gqcov
        fp = bool(hit) and hit_taxid in catalog_desc
        rows.append({
            "qid": meta["qid"], "taxid": meta["taxid"], "species": meta["species"],
            "category": meta["category"], "hit_taxid": hit_taxid,
            "hit_species_taxid": hit_species, "hit_genus_taxid": hit_genus,
            "pident": round(pid, 1), "qcovs": round(qcov, 1),
            "species_retrieval_ok": species_retrieval,
            "genus_retrieval_ok": genus_retrieval,
            "species_class_ok": species_class,
            "genus_class_ok": genus_class,
            "catalog_false_positive": fp,
        })
    summary = {}
    for cat in CATEGORY_ORDER:
        subset = [r for r in rows if r["category"] == cat]
        if not subset:
            continue
        summary[cat] = {}
        for key in ("species_retrieval_ok", "genus_retrieval_ok", "species_class_ok", "genus_class_ok"):
            correct = sum(1 for r in subset if r[key])
            summary[cat][key] = {
                "n": len(subset), "correct": correct,
                "accuracy": round(correct / len(subset) * 100, 1),
            }
        catalog_hits = sum(1 for r in subset if r["catalog_false_positive"])
        summary[cat]["catalog_hit_rate"] = {
            "n": len(subset), "correct": catalog_hits,
            "accuracy": round(catalog_hits / len(subset) * 100, 1),
        }
    return rows, summary


def summarize(hits_by_cat, nodes, names, merged, catalog):
    """Keep the previous interface for callers that pass category buckets."""
    records = []
    for cat, items in hits_by_cat.items():
        records.extend(items)
    return summarize_group(records, nodes, names, merged, catalog)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--categories", default="bacteria,viruses,fungi,oomycetes")
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--taxidlist", action="store_true", help="restrict each category to its NCBI lineage")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--analyze-tsv", help="post-process an existing BLAST outfmt-6 TSV instead of running BLAST")
    args = ap.parse_args()
    categories = [c for c in args.categories.split(",") if c in CATEGORY_ORDER]

    nodes, names, merged = load_taxdump()
    catalog = load_catalog(nodes, merged)

    pos_records = [r for r in read_fasta(os.path.join(QUERY_DIR, "external_positives.fasta"))
                   if parse_header(r["header"])["category"] in categories]
    neg_records = [r for r in read_fasta(os.path.join(QUERY_DIR, "negatives.fasta"))
                   if parse_header(r["header"])["category"] in categories]
    all_records = pos_records + neg_records

    if args.analyze_tsv:
        records = []
        hit_by_qid = {}
        with open(args.analyze_tsv, encoding="utf-8") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                qid = f[0].split("|", 1)[0]
                staxids = f[2]
                first_taxid = staxids.split(";")[0] if staxids else ""
                if qid not in hit_by_qid:
                    hit_by_qid[qid] = {
                        "qid": qid,
                        "hit": {
                            "sseqid": f[1], "staxid": first_taxid,
                            "pident": float(f[3]), "qcovs": float(f[4]),
                        },
                    }
        for rec in pos_records + neg_records:
            meta = parse_header(rec["header"])
            records.append({
                "header": rec["header"], "qid": meta["qid"],
                "hit": (hit_by_qid.get(meta["qid"]) or {}).get("hit"),
            })
        pos_rows, pos_summary = summarize_group(
            [r for r in records if r["qid"].startswith("e")], nodes, names, merged, catalog
        )
        neg_rows, neg_summary = summarize_group(
            [r for r in records if r["qid"].startswith("n")], nodes, names, merged, catalog
        )
        rows, summary = summarize_group(records, nodes, names, merged, catalog)
        result = {
            "validation": "ncbi_nt_comparison",
            "database": "local NCBI-nt snapshot (2026-08-23)",
            "snapshot_date": "2026-08-23",
            "thresholds": THRESHOLDS,
            "search_restriction": "full NCBI-nt (no lineage restriction)",
            "source": os.path.relpath(args.analyze_tsv, RESULT_DIR),
            "positives": {"summary": pos_summary, "details": pos_rows},
            "negatives": {"summary": neg_summary, "details": neg_rows},
            "combined_summary": summary,
            "combined_details": rows,
        }
        out_path = os.path.join(RESULT_DIR, "ncbi_nt_comparison.json")
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(result, fh, ensure_ascii=False, indent=2)
        print("summary:", {k: v for k, v in summary.items()})
        print("wrote", out_path)
        return

    taxidlists = build_taxidlists(all_records, nodes, merged) if args.taxidlist else {}

    jobs = []
    for cat in categories:
        taxidlist = write_taxidlist(cat, taxidlists[cat]) if args.taxidlist else None
        cat_records = [r for r in pos_records + neg_records if parse_header(r["header"])["category"] == cat]
        for batch_id, fasta in split_fasta(cat_records, args.batch_size):
            jobs.append((cat, batch_id, fasta, taxidlist))
    if args.dry_run:
        print(f"dry run: {len(jobs)} batches across {categories}")
        for cat, bid, fasta, taxidlist in jobs[:5]:
            print(cat, bid, fasta)
        return

    done = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(blast_batch, bid, fasta, cat, taxidlist) for cat, bid, fasta, taxidlist in jobs]
        for fut in concurrent.futures.as_completed(futures):
            bid, out_path = fut.result()
            done[bid] = out_path
            print(f"done batch {bid}: {out_path}", flush=True)

    # Load hits by category from all output files and attach them to every
    # query record so queries without a hit are still counted.
    hit_by_qid = {}
    for cat in categories:
        for path in sorted(os.listdir(WORK_DIR)):
            if not path.startswith(f"hits_{cat}_"):
                continue
            full = os.path.join(WORK_DIR, path)
            if not os.path.getsize(full):
                continue
            with open(full, encoding="utf-8") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 5:
                        continue
                    qid = f[0].split("|", 1)[0]
                    staxids = f[2]
                    first_taxid = staxids.split(";")[0] if staxids else ""
                    hit_by_qid[qid] = {
                        "qid": qid,
                        "hit": {
                            "sseqid": f[1], "staxid": first_taxid,
                            "pident": float(f[3]), "qcovs": float(f[4]),
                        },
                    }

    hits_by_cat = defaultdict(list)
    for rec in pos_records + neg_records:
        meta = parse_header(rec["header"])
        hits_by_cat[meta["category"]].append({
            "header": rec["header"],
            "qid": meta["qid"],
            "hit": (hit_by_qid.get(meta["qid"]) or {}).get("hit"),
        })

    rows, summary = summarize(hits_by_cat, nodes, names, merged, catalog)
    result = {
        "validation": "ncbi_nt_comparison",
        "database": "local NCBI-nt snapshot (2026-08-23)",
        "snapshot_date": "2026-08-23",
        "thresholds": THRESHOLDS,
        "search_restriction": "per-category NCBI lineage via -taxidlist" if args.taxidlist else "none (full NCBI-nt)",
        "summary": summary,
        "details": rows,
    }
    with open(os.path.join(RESULT_DIR, "ncbi_nt_comparison.json"), "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False, indent=2)
    print("summary:", {k: v for k, v in summary.items()})
    print("wrote", os.path.join(RESULT_DIR, "ncbi_nt_comparison.json"))


if __name__ == "__main__":
    sys.exit(main())
