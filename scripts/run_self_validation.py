#!/usr/bin/env python3
"""Internal self-hit validation for current MaizePathogenDB build (573 sequences)."""
import os, json, subprocess, tempfile, collections

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
BLAST = "/Users/gfgao/Desktop/blacksoil_metaG/tools/ncbi-blast-2.17.0+/bin/blastn"
FASTA = os.path.join(BASE, "sequences", "maize_pathogens_all.fasta")
OUT = os.path.join(BASE, "docs", "validation")

def parse_fasta(path):
    records = []
    cur = None
    for line in open(path):
        if line.startswith(">"):
            if cur:
                records.append(cur)
            parts = line[1:].strip().split("|", 3)
            cur = {"taxid": parts[0], "species": parts[1] if len(parts)>1 else "",
                   "cat": parts[2] if len(parts)>2 else "", "seq": "", "header": line[1:].strip()}
        elif cur is not None:
            cur["seq"] += line.strip()
    if cur:
        records.append(cur)
    return records

records = parse_fasta(FASTA)
print(f"Parsed {len(records)} sequences", flush=True)

db_map = {
    "bacteria": "maize_pathogens_bacteria",
    "viruses": "maize_pathogens_viruses",
    "fungi": "maize_pathogens_fungi",
}

# Group by category and write one query file per category
by_cat = collections.defaultdict(list)
for rec in records:
    by_cat[rec["cat"]].append(rec)

all_results = {}
for cat, db_name in db_map.items():
    recs = by_cat.get(cat, [])
    if not recs:
        continue
    qfile = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    for j, rec in enumerate(recs, 1):
        qfile.write(f">q{j}\n{rec['seq']}\n")
    qfile.close()
    db_path = os.path.join(BASE, "blast_db", db_name)
    cmd = f'"{BLAST}" -query "{qfile.name}" -db "{db_path}" -outfmt "6 qseqid sseqid pident length" -max_target_seqs 1 -num_threads 4'
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    os.unlink(qfile.name)

    hit_map = {}
    for line in proc.stdout.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 4:
            qid = parts[0]
            hit_id = parts[1].split("|")[0].split(".")[0]
            pident = float(parts[2])
            if qid not in hit_map:
                hit_map[qid] = (hit_id, pident)

    correct = 0
    details = []
    for j, rec in enumerate(recs, 1):
        qid = f"q{j}"
        hit = hit_map.get(qid)
        ok = bool(hit and hit[0] == rec["taxid"])
        if ok:
            correct += 1
        details.append({"query_taxid": rec["taxid"], "hit_taxid": hit[0] if hit else None,
                        "pident": hit[1] if hit else 0, "correct": ok})
    acc = correct / len(recs) * 100
    all_results[cat] = {"n": len(recs), "correct": correct, "accuracy": round(acc, 1), "details": details}
    print(f"{cat}: {correct}/{len(recs)} ({acc:.1f}%)", flush=True)

total_n = sum(v["n"] for v in all_results.values())
total_c = sum(v["correct"] for v in all_results.values())
all_results["overall"] = {"n": total_n, "correct": total_c, "accuracy": round(total_c/total_n*100, 1)}
print(f"Overall: {total_c}/{total_n} ({all_results['overall']['accuracy']:.1f}%)", flush=True)

out_path = os.path.join(OUT, "self_validation_current.json")
with open(out_path, "w") as f:
    json.dump({"validation_type": "internal_self_hit", "date": "2026-08-18",
               "sequences": total_n, "results": all_results}, f, indent=2, ensure_ascii=False)
print(f"Saved: {out_path}")
