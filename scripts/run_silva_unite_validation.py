#!/usr/bin/env python3
"""SILVA/UNITE cross-validation against current MaizePathogenDB build."""
import os, json, subprocess, tempfile, collections

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
BLAST = "/Users/gfgao/Desktop/blacksoil_metaG/tools/ncbi-blast-2.17.0+/bin/blastn"
QUERY = os.path.join(BASE, "docs", "validation", "external", "silva_unite_cross_queries.fasta")
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
                   "cat": (parts[2].split()[0] if len(parts)>2 else ""), "seq": "", "header": line[1:].strip()}
        elif cur is not None:
            cur["seq"] += line.strip()
    if cur:
        records.append(cur)
    return records

records = parse_fasta(QUERY)
print(f"Parsed {len(records)} SILVA/UNITE queries", flush=True)

# Deduplicate: keep one representative sequence per taxid
seen = set()
deduped = []
for rec in records:
    key = rec["taxid"]
    if key not in seen:
        seen.add(key)
        deduped.append(rec)
records = deduped
print(f"Deduplicated to {len(records)} queries (one per taxid)", flush=True)

# Only evaluate query taxa that have reference sequences in the current database
with open(os.path.join(BASE, "..", "maize_pathogen_web", "data", "seq_meta.json")) as f:
    seq_meta = json.load(f)
valid_taxids = set(seq_meta.keys())
before = len(records)
records = [r for r in records if r["taxid"] in valid_taxids]
print(f"Filtered to {len(records)} queries present in current database (removed {before-len(records)})", flush=True)

db_map = {"bacteria": "maize_pathogens_bacteria", "fungi": "maize_pathogens_fungi"}
by_cat = collections.defaultdict(list)
for rec in records:
    cat = rec["cat"]
    if cat in db_map:
        by_cat[cat].append(rec)

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
        details.append({"query_taxid": rec["taxid"], "query_species": rec["species"],
                        "hit_taxid": hit[0] if hit else None,
                        "pident": hit[1] if hit else 0, "correct": ok})
    acc = correct / len(recs) * 100
    all_results[cat] = {"n": len(recs), "correct": correct, "accuracy": round(acc, 1), "details": details}
    print(f"{cat}: {correct}/{len(recs)} ({acc:.1f}%)", flush=True)

total_n = sum(v["n"] for v in all_results.values())
total_c = sum(v["correct"] for v in all_results.values())
all_results["overall"] = {"n": total_n, "correct": total_c, "accuracy": round(total_c/total_n*100, 1)}
print(f"Overall: {total_c}/{total_n} ({all_results['overall']['accuracy']:.1f}%)", flush=True)

out_path = os.path.join(OUT, "cross_validation_current.json")
with open(out_path, "w") as f:
    json.dump({"validation_type": "cross_validation_silva_unite", "date": "2026-08-18",
               "queries": total_n, "results": all_results}, f, indent=2, ensure_ascii=False)
print(f"Saved: {out_path}")
