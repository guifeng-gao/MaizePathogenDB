#!/usr/bin/env python3
"""
MaizePathogenDB vs NCBI-nt comparison (v2).
Scale from 15 to ~100 sequences with stratified random sampling.
"""

import os, json, sys, time, subprocess, random, re, tempfile
from collections import defaultdict
import requests
from Bio.Blast import NCBIWWW
import matplotlib
matplotlib.use('Agg')
from Bio.Blast import NCBIWWW
import matplotlib.pyplot as plt
import numpy as np

# ── Config ──────────────────────────────────────────────────────────
BASE     = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
FASTA    = os.path.join(BASE, "sequences", "maize_pathogens_all.fasta")
BLAST_DB = os.path.join(BASE, "blast_db", "maize_pathogens_all")
BLAST_BIN = "/tmp/ncbi-blast-2.17.0+/bin/blastn"
OUT_DIR  = os.path.join(BASE, "docs", "validation")
os.makedirs(OUT_DIR, exist_ok=True)

RESULTS_FILE = os.path.join(OUT_DIR, "ncbi_nt_comparison_v2.json")
FIGURE_FILE   = os.path.join(OUT_DIR, "Fig_NCBI_nt_Comparison_v2.pdf")

SAMPLE_N      = 200
RANDOM_SEED   = 42
NCBI_DELAY    = 3.5          # seconds between NCBI web BLAST requests
NCBI_TIMEOUT  = 600          # max seconds to wait for a single BLAST job
NCBI_POLL     = 5            # poll interval (seconds)
NCBI_EMAIL    = "maize_pathogen_db@example.com"

CAT_ORDER = ["bacteria", "viruses", "fungi"]
CAT_LABELS = {"bacteria": "Bacteria (16S)", "viruses": "Viruses (Genome)", "fungi": "Fungi (ITS)"}
CAT_COLORS = {"bacteria": "#2F5496", "viruses": "#C00000", "fungi": "#548235"}

# ── Parse FASTA ─────────────────────────────────────────────────────
def parse_fasta(fpath):
    records = []
    current = None
    for line in open(fpath):
        if line.startswith(">"):
            if current:
                records.append(current)
            parts = line[1:].strip().split("|", 3)
            taxid = parts[0]
            species = parts[1] if len(parts) > 1 else ""
            category = parts[2] if len(parts) > 2 else "unknown"
            current = {"header": line[1:].strip(), "seq": "",
                       "taxid": taxid, "species": species,
                       "category": category}
        elif current is not None:
            current["seq"] += line.strip()
    if current:
        records.append(current)
    return records

# ── Local BLAST ─────────────────────────────────────────────────────
def run_local_blast(query_seq):
    """Run blastn against MaizePathogenDB, return top hit dict."""
    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    tf.write(">query\n" + query_seq + "\n")
    tf.close()
    cmd = f'"{BLAST_BIN}" -query "{tf.name}" -db "{BLAST_DB}" -outfmt "6 sseqid pident length" -max_target_seqs 1 -num_threads 2'
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
    os.unlink(tf.name)
    lines = [l.strip() for l in proc.stdout.split("\n") if l.strip()]
    if not lines:
        return None
    line = lines[-1]  # last non-empty line (skip BLAST warnings)
    parts = line.split("\t")
    if len(parts) < 3:
        print(f"  WARN: unexpected BLAST output: {repr(line)}", flush=True)
        return None
    return {"sseqid": parts[0], "pident": float(parts[1]), "length": int(parts[2])}

# ── NCBI Web BLAST ──────────────────────────────────────────────────
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
        # For XML format, check if the response contains BLAST results
        if r2.text.startswith("<?xml") or "<BlastOutput" in r2.text:
            return r2.text
    return "TIMEOUT"

def extract_top_hit(xml_data):
    """Parse top hit from NCBI BLAST XML."""
    if isinstance(xml_data, str) and (xml_data.startswith("BLAST_ERROR") or xml_data.startswith("NO_") or xml_data.startswith("TIMEOUT")):
        return {"description": xml_data, "pident": 0.0}
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(xml_data)
        hits = root.findall(".//Hit")
        if not hits:
            return {"description": "NO_HITS", "pident": 0.0}
        top = hits[0]
        desc = top.findtext("Hit_def", "UNKNOWN")
        hsps = top.findall(".//Hsp")
        if not hsps:
            return {"description": desc, "pident": 0.0}
        identity = int(hsps[0].findtext("Hsp_identity", "0"))
        align_len = int(hsps[0].findtext("Hsp_align-len", "1"))
        pident = identity / align_len * 100 if align_len else 0.0
        return {"description": desc, "pident": pident}
    except ET.ParseError as e:
        return {"description": f"XML_PARSE_ERROR: {e}", "pident": 0.0}
    except Exception as e:
        return {"description": f"XML_ERROR: {e}", "pident": 0.0}

def check_ncbi_species(description, expected_species):
    """Genus-level match check."""
    if not description or description.startswith("ERROR") or description.startswith("NO_"):
        return False
    expected_words = expected_species.split()
    if not expected_words:
        return False
    expected_genus = expected_words[0].lower()
    desc_lower = description.lower()
    return expected_genus in desc_lower

# ═══════════════════════════════════════════════════════════════════
print("=" * 60)
print("MaizePathogenDB vs NCBI-nt Comparison  v2")
print("=" * 60)

# ── Step 1: Parse & Sample ──────────────────────────────────────────
random.seed(RANDOM_SEED)
all_records = parse_fasta(FASTA)
print(f"\nTotal sequences: {len(all_records)}")

by_cat_species = defaultdict(list)
for r in all_records:
    if r["category"] in CAT_ORDER:
        by_cat_species[(r["category"], r["species"])].append(r)

for cat in CAT_ORDER:
    n = len([k for k in by_cat_species if k[0] == cat])
    total = len([r for r in all_records if r["category"] == cat])
    print(f"  {cat}: {n} species, {total} sequences")

total_species = sum(len([k for k in by_cat_species if k[0] == cat]) for cat in CAT_ORDER)
sample = {}
for cat in CAT_ORDER:
    cat_species = [k for k in by_cat_species if k[0] == cat]
    prop = max(15, int(round(len(cat_species) / total_species * SAMPLE_N)))
    chosen = random.sample(cat_species, min(prop, len(cat_species)))
    for k in chosen:
        sample[by_cat_species[k][0]["header"]] = random.choice(by_cat_species[k])
    print(f"  Sampled {len(chosen)}/{len(cat_species)} species from {cat}")

sample_list = list(sample.values())
total = len(sample_list)
print(f"\nTotal test queries: {total}")
for cat in CAT_ORDER:
    print(f"  {cat}: {sum(1 for r in sample_list if r['category'] == cat)}")

# ── Step 2: Local BLAST ─────────────────────────────────────────────
print(f"\n{'=' * 60}")
print("Phase 1: Local BLAST against MaizePathogenDB")
print("=" * 60)

for i, rec in enumerate(sample_list):
    result = run_local_blast(rec["seq"])
    if result:
        top_taxid = result["sseqid"].split("|")[0].split(".")[0]
        rec["db_correct"] = (top_taxid == rec["taxid"])
        rec["db_pident"] = result["pident"]
        rec["db_length"] = result["length"]
        rec["db_top"] = result["sseqid"]
    else:
        rec["db_correct"] = False
        rec["db_pident"] = 0.0
        rec["db_length"] = 0
        rec["db_top"] = "NO HITS"
    if (i+1) % 20 == 0:
        print(f"  [{i+1}/{total}] local BLAST done")

db_ok = sum(1 for r in sample_list if r.get("db_correct"))
print(f"  MaizePathogenDB: {db_ok}/{total} correct ({db_ok/total*100:.1f}%)")

# ── Step 3: NCBI Web BLAST ──────────────────────────────────────────
print(f"\n{'=' * 60}")
print("Phase 2: NCBI web BLAST against nt database")
print("(this will take ~5-15 minutes)")
print("=" * 60)

# Load cached results
done_headers = set()
cached = {"config": {"sample_n": total, "seed": RANDOM_SEED}, "results": []}
if os.path.exists(RESULTS_FILE):
    with open(RESULTS_FILE) as f:
        cached_data = json.load(f)
    if "results" in cached_data:
        done_headers = {c["header"] for c in cached_data["results"]}
        cached = cached_data
    print(f"  Loaded {len(done_headers)} cached results from {RESULTS_FILE}")

idx = 0
for rec in sample_list:
    idx += 1
    if rec["header"] in done_headers:
        print(f"  [{idx}/{total}] Skipped (cached): {rec['species'][:40]}")
        continue

    print(f"  [{idx}/{total}] BLASTing: {rec['species'][:50]} ({rec['category']}) ...", end=" ", flush=True)

    # Retry up to 3 times for XML errors
    hit = None
    for attempt in range(3):
        blast_result = submit_ncbi_blast(rec["seq"])
        hit = extract_top_hit(blast_result)
        if "XML_PARSE_ERROR" not in str(hit["description"]) and "XML_ERROR" not in str(hit["description"]):
            break
        if attempt < 2:
            print(f"retry({attempt+1}) ", end="", flush=True)
            time.sleep(10)  # longer wait before retry

    rec["ncbi_description"] = hit["description"]
    rec["ncbi_pident"] = hit["pident"]
    rec["ncbi_correct"] = check_ncbi_species(hit["description"], rec["species"])

    status = "✓" if rec["ncbi_correct"] else "✗"
    print(f"{status} (pident={hit['pident']:.1f}%)")

    cached["results"].append({
        "header": rec["header"],
        "species": rec["species"],
        "category": rec["category"],
        "db_correct": rec.get("db_correct"),
        "db_pident": rec.get("db_pident", 0),
        "ncbi_correct": rec.get("ncbi_correct"),
        "ncbi_pident": rec.get("ncbi_pident", 0),
        "ncbi_description": rec.get("ncbi_description", ""),
    })
    with open(RESULTS_FILE, "w") as f:
        json.dump(cached, f, indent=2, ensure_ascii=False)

    time.sleep(NCBI_DELAY)

# ── Step 4: Analysis ────────────────────────────────────────────────
print(f"\n{'=' * 60}")
print("Final Results")
print("=" * 60)

results_by_cat = defaultdict(list)
for r in sample_list:
    results_by_cat[r["category"]].append(r)

summary = {}
for cat in CAT_ORDER:
    data = results_by_cat[cat]
    n = len(data)
    db_c = sum(1 for r in data if r.get("db_correct"))
    ncbi_c = sum(1 for r in data if r.get("ncbi_correct"))
    db_acc = db_c / n * 100 if n else 0
    ncbi_acc = ncbi_c / n * 100 if n else 0
    summary[cat] = {"n": n, "db_correct": db_c, "ncbi_correct": ncbi_c,
                    "db_accuracy": round(db_acc, 1), "ncbi_accuracy": round(ncbi_acc, 1)}
    print(f"\n  {CAT_LABELS[cat]} (n={n}):")
    print(f"    MaizePathogenDB: {db_c}/{n} ({db_acc:.1f}%)")
    print(f"    NCBI-nt:         {ncbi_c}/{n} ({ncbi_acc:.1f}%)")

n_total = total
db_total = sum(1 for r in sample_list if r.get("db_correct"))
ncbi_total = sum(1 for r in sample_list if r.get("ncbi_correct"))
summary["overall"] = {"n": n_total, "db_correct": db_total, "ncbi_correct": ncbi_total,
                      "db_accuracy": round(db_total/n_total*100, 1) if n_total else 0,
                      "ncbi_accuracy": round(ncbi_total/n_total*100, 1) if n_total else 0}
print(f"\n  Overall (n={n_total}):")
print(f"    MaizePathogenDB: {db_total}/{n_total} ({db_total/n_total*100:.1f}%)")
print(f"    NCBI-nt:         {ncbi_total}/{n_total} ({ncbi_total/n_total*100:.1f}%)")

summary["results"] = [
    {"species": r["species"], "category": r["category"],
     "db_correct": r.get("db_correct"), "db_pident": r.get("db_pident", 0),
     "ncbi_correct": r.get("ncbi_correct"), "ncbi_pident": r.get("ncbi_pident", 0),
     "ncbi_top": r.get("ncbi_description", "")}
    for r in sample_list
]
with open(RESULTS_FILE, "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"\nResults saved to: {RESULTS_FILE}")

# ── Step 5: Generate Figure ─────────────────────────────────────────
print(f"\n{'=' * 60}")
print("Generating comparison figure...")
print("=" * 60)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 150, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

fig = plt.figure(figsize=(14, 8))

# Panel A: Grouped bar chart
ax1 = fig.add_subplot(2, 1, 1)
categories = [CAT_LABELS[c] for c in CAT_ORDER]
db_accs = [summary[c]["db_accuracy"] for c in CAT_ORDER]
ncbi_accs = [summary[c]["ncbi_accuracy"] for c in CAT_ORDER]

x = np.arange(len(categories))
width = 0.32

bars1 = ax1.bar(x - width/2, db_accs, width, label='MaizePathogenDB',
                color='#2F5496', edgecolor='white', linewidth=1.5)
bars2 = ax1.bar(x + width/2, ncbi_accs, width, label='NCBI-nt',
                color='#C00000', edgecolor='white', linewidth=1.5)

for bar, v in zip(bars1, db_accs):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
             f'{v:.1f}%', ha='center', fontsize=10, fontweight='bold', color='#2F5496')
for bar, v in zip(bars2, ncbi_accs):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
             f'{v:.1f}%', ha='center', fontsize=10, fontweight='bold', color='#C00000')

for i, cat in enumerate(CAT_ORDER):
    n = summary[cat]["n"]
    ax1.text(i, 2, f'n={n}', ha='center', fontsize=8, color='gray')

ax1.set_xticks(x)
ax1.set_xticklabels(categories, fontsize=10)
ax1.set_ylabel('Top-1 Accuracy (%)', fontweight='bold')
ax1.set_ylim(0, 108)
ax1.set_title('A. MaizePathogenDB vs NCBI-nt: Classification Accuracy', fontweight='bold', loc='left')
ax1.legend(loc='lower right', frameon=True, fancybox=True, framealpha=0.9)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.spines[['right', 'top']].set_visible(False)

# Panel B: Key findings summary table
ax2 = fig.add_subplot(2, 1, 2)
ax2.axis('off')

detail_data = [
    ['Category', 'Queries', 'MaizePathogenDB', 'NCBI-nt', 'Key Finding'],
]
for cat in CAT_ORDER:
    s = summary[cat]
    n = s["n"]
    db_a = s["db_accuracy"]
    ncbi_a = s["ncbi_accuracy"]
    detail_data.append([
        CAT_LABELS[cat], str(n), f'{db_a:.1f}%', f'{ncbi_a:.1f}%',
        f'MaizePathogenDB {db_a:.0f}% vs NCBI-nt {ncbi_a:.0f}%'
    ])
detail_data.append([
    'OVERALL', str(summary['overall']['n']),
    f'{summary["overall"]["db_accuracy"]:.1f}%',
    f'{summary["overall"]["ncbi_accuracy"]:.1f}%',
    f'MaizePathogenDB {summary["overall"]["db_accuracy"]:.0f}% vs NCBI-nt {summary["overall"]["ncbi_accuracy"]:.0f}%'
])

table = ax2.table(cellText=detail_data, cellLoc='center', loc='center',
                  colWidths=[0.18, 0.10, 0.20, 0.20, 0.32])
table.auto_set_font_size(False)
table.set_fontsize(10)
for i in range(len(detail_data)):
    for j in range(5):
        cell = table[i, j]
        if i == 0:
            cell.set_facecolor('#333333')
            cell.set_text_props(color='white', fontweight='bold')
        elif i == len(detail_data) - 1:
            cell.set_facecolor('#E8E8E8')
            cell.set_text_props(fontweight='bold')

ax2.set_title(f'B. Summary Comparison (n={summary["overall"]["n"]} queries, '
              f'NCBI-nt wrong: {ncbi_total}/{n_total}, MaizePathogenDB wrong: {db_total}/{n_total})',
              fontweight='bold', loc='left', fontsize=10)

fig.suptitle('MaizePathogenDB v1.0 vs NCBI-nt Database: Sequence Classification Accuracy',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(FIGURE_FILE, facecolor='white', edgecolor='none')
plt.close()
print(f"  ✓ {FIGURE_FILE}")

print(f"\n{'=' * 60}")
print("DONE")
print(f"{'=' * 60}")
print(f"Results: {RESULTS_FILE}")
print(f"Figure:  {FIGURE_FILE}")
print(f"\nKey numbers:")
print(f"  MaizePathogenDB: {summary['overall']['db_accuracy']:.1f}% ({db_total}/{n_total})")
print(f"  NCBI-nt:         {summary['overall']['ncbi_accuracy']:.1f}% ({ncbi_total}/{n_total})")
print(f"  NCBI-nt errors:  {n_total - ncbi_total}/{n_total}")
print(f"  DB errors:       {n_total - db_total}/{n_total}")
