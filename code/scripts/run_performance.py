#!/usr/bin/env python3
"""Performance measurement on fixed ASV query sets (Usage Notes only)."""

import json
import os
import random
import re
import subprocess
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASV_FASTA = os.path.join(ROOT, "260samples_fungi", "rep-seqs.fasta")
OUT_DIR = os.path.join(ROOT, "data", "performance_query_sets")
RESULT_DIR = os.environ.get(
    "RESULT_DIR", os.path.join(ROOT, "docs", "validation", "results")
)
MPDB_DB = os.path.join(ROOT, "release", "blast_db", "maize_pathogens_all")
ITS_EUK_DB = os.path.join(
    ROOT, "260samples_fungi", "analysis", "db",
    "ncbi_ITS_eukaryote", "ITS_eukaryote_sequences",
)
DB_PATHS = {
    "MPDB": MPDB_DB,
    "NCBI_ITS_eukaryote": ITS_EUK_DB,
}
DB_LABELS = {
    "MPDB": "release/blast_db/maize_pathogens_all",
    "NCBI_ITS_eukaryote": "260samples_fungi/analysis/db/ncbi_ITS_eukaryote/ITS_eukaryote_sequences",
}
BLASTN = os.environ.get("BLASTN", "blastn")
SEED = 42
SIZES = [100, 500, 1000, 2000, 5000, 10000, 19276]


def load_asvs():
    records = []
    header = None
    seq = []
    with open(ASV_FASTA, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header:
            records.append((header, "".join(seq)))
    return records


def measure(query_file, db):
    cmd = [
        "/usr/bin/time", "-l",
        BLASTN, "-query", query_file, "-db", db,
        "-outfmt", "6 qseqid sseqid pident", "-max_target_seqs", "1",
        "-evalue", "1e-5", "-num_threads", "4",
    ]
    start = time.perf_counter()
    result = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                            text=True, timeout=7200)
    elapsed = time.perf_counter() - start
    match = re.search(r"(\d+)\s+maximum resident set size", result.stderr)
    maxrss_bytes = int(match.group(1)) if match else 0
    return {
        "wall_s": round(elapsed, 3),
        "maxrss_mb": round(maxrss_bytes / 1024 / 1024, 1),
    }


def main():
    records = load_asvs()
    rng = random.Random(SEED)
    rng.shuffle(records)
    os.makedirs(OUT_DIR, exist_ok=True)
    results = {
        "note": "Usage Notes only; fixed ASV query set from 260samples_fungi/rep-seqs.fasta; no pathogen detection analysis",
        "seed": SEED,
        "databases": DB_LABELS,
        "rows": [],
    }
    for size in SIZES:
        subset = records[:size]
        qfile = os.path.join(OUT_DIR, f"asv_query_{size}.fasta")
        with open(qfile, "w", encoding="utf-8") as fh:
            for header, seq in subset:
                fh.write(">" + header + "\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
        for db_name, db in DB_PATHS.items():
            metric = measure(qfile, db)
            results["rows"].append({
                "n_queries": size,
                "database": db_name,
                **metric,
                "per_query_ms": round(metric["wall_s"] * 1000 / size, 3),
            })
            print(size, db_name, metric, flush=True)

    os.makedirs(RESULT_DIR, exist_ok=True)
    out = os.path.join(RESULT_DIR, "performance.json")
    json.dump(results, open(out, "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
