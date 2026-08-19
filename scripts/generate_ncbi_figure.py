#!/usr/bin/env python3
"""Generate publication-quality NCBI-nt comparison figure from JSON results."""

import os, json, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = "/Users/gfgao/Desktop/blacksoil_metaG/maize_pathogen_db"
JSON_FILE = os.path.join(BASE, "docs", "validation", "ncbi_nt_comparison_v2.json")
OUT_DIR = os.path.join(BASE, "docs", "validation")

with open(JSON_FILE) as f:
    data = json.load(f)

C = {"bacteria": "#2F5496", "viruses": "#C00000", "fungi": "#548235",
     "oomycete": "#BF8F00", "dark_gray": "#333333", "medium_gray": "#666666", "light_gray": "#E8E8E8"}

CAT_ORDER = ["bacteria", "viruses", "fungi"]
CAT_LABELS = {"bacteria": "Bacteria (16S)", "viruses": "Viruses (Genome)", "fungi": "Fungi (ITS)"}

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 10,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 150, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

fig = plt.figure(figsize=(14, 9))

# Panel A: Grouped bar chart
ax1 = fig.add_subplot(2, 1, 1)
categories = [CAT_LABELS[c] for c in CAT_ORDER]
db_accs = [data[c]["db_accuracy"] for c in CAT_ORDER]
ncbi_accs = [data[c]["ncbi_accuracy"] for c in CAT_ORDER]
db_counts = [f"{data[c]['db_correct']}/{data[c]['n']}" for c in CAT_ORDER]
ncbi_counts = [f"{data[c]['ncbi_correct']}/{data[c]['n']}" for c in CAT_ORDER]

x = np.arange(len(categories))
width = 0.32

bars1 = ax1.bar(x - width/2, db_accs, width, label='MaizePathogenDB',
                color=C["bacteria"], edgecolor='white', linewidth=1.5)
bars2 = ax1.bar(x + width/2, ncbi_accs, width, label='NCBI-nt',
                color=C["viruses"], edgecolor='white', linewidth=1.5)

for bar, v, cnt in zip(bars1, db_accs, db_counts):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
             f'{v:.1f}%\n({cnt})', ha='center', fontsize=9, fontweight='bold', color=C["bacteria"])
for bar, v, cnt in zip(bars2, ncbi_accs, ncbi_counts):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
             f'{v:.1f}%\n({cnt})', ha='center', fontsize=9, fontweight='bold', color=C["viruses"])

ax1.set_xticks(x)
ax1.set_xticklabels(categories, fontsize=10)
ax1.set_ylabel('Top-1 Accuracy (%)', fontweight='bold')
ax1.set_ylim(0, 112)
ax1.set_title('A. MaizePathogenDB vs NCBI-nt: Classification Accuracy',
              fontweight='bold', loc='left')
ax1.legend(loc='lower right', frameon=True, fancybox=True, framealpha=0.9)
ax1.grid(axis='y', alpha=0.3, linestyle='--')
ax1.spines[['right', 'top']].set_visible(False)

# Panel B: Summary table
ax2 = fig.add_subplot(2, 1, 2)
ax2.axis('off')

n_total = data["overall"]["n"]
detail_data = [
    ['Category', 'n', 'MaizePathogenDB', 'NCBI-nt', 'Difference'],
]
for cat in CAT_ORDER:
    s = data[cat]
    diff = s["db_accuracy"] - s["ncbi_accuracy"]
    detail_data.append([
        CAT_LABELS[cat], str(s["n"]),
        f'{s["db_correct"]}/{s["n"]} ({s["db_accuracy"]:.1f}%)',
        f'{s["ncbi_correct"]}/{s["n"]} ({s["ncbi_accuracy"]:.1f}%)',
        f'+{diff:.1f}%' if diff > 0 else f'{diff:.1f}%',
    ])

s = data["overall"]
diff = s["db_accuracy"] - s["ncbi_accuracy"]
detail_data.append([
    'OVERALL', str(s["n"]),
    f'{s["db_correct"]}/{s["n"]} ({s["db_accuracy"]:.1f}%)',
    f'{s["ncbi_correct"]}/{s["n"]} ({s["ncbi_accuracy"]:.1f}%)',
    f'+{diff:.1f}%' if diff > 0 else f'{diff:.1f}%',
])

# Key findings note
failures = [r for r in data["results"] if not r["ncbi_correct"]]
fail_species = ", ".join(r["species"][:35] for r in failures)

col_widths = [0.22, 0.07, 0.22, 0.22, 0.22]
table = ax2.table(cellText=detail_data, cellLoc='center', loc='center',
                  colWidths=col_widths)
table.auto_set_font_size(False)
table.set_fontsize(9)
for i in range(len(detail_data)):
    for j in range(5):
        cell = table[i, j]
        if i == 0:
            cell.set_facecolor(C["dark_gray"])
            cell.set_text_props(color='white', fontweight='bold')
        elif i == len(detail_data) - 1:
            cell.set_facecolor(C["light_gray"])
            cell.set_text_props(fontweight='bold')

# Failure note below table
db_fail = sum(1 for r in data["results"] if not r["db_correct"])
ax2.text(0.5, 0.12,
         f'Data: stratified random sample (n={n_total}; 1 seq/species)  |  '
         f'NCBI-nt errors: {len(failures)}  |  '
         f'MaizePathogenDB errors: {db_fail}',
         ha='center', va='top', fontsize=7, color=C["medium_gray"],
         transform=ax2.transAxes, style='italic')

ax2.set_title('B. Summary Comparison', fontweight='bold', loc='left', fontsize=11)

fig.suptitle('MaizePathogenDB v1.0 vs NCBI-nt Database: Sequence Classification Accuracy',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()

OUTPUT = os.path.join(OUT_DIR, "Fig_NCBI_nt_Comparison_Final_v2.pdf")
fig.savefig(OUTPUT, facecolor='white', edgecolor='none')
plt.close()
print(f"Figure saved: {OUTPUT}")

# Copy to replace old figure references
for name in ["Fig_NCBI_nt_Comparison.pdf", "Fig_NCBI_nt_Comparison_Final.pdf"]:
    dest = os.path.join(OUT_DIR, name)
    shutil.copy2(OUTPUT, dest)
    print(f"  -> {dest}")
