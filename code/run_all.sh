#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python3}"
MAKEBLASTDB="${MAKEBLASTDB:-makeblastdb}"

echo "==> MaizePathogenDB release reproducibility pipeline"

echo "==> 1/8 species_list"
MAIZE_FASTA="$ROOT/release/sequences/maize_pathogens_all.fasta" \
  "$PY" "$ROOT/scripts/build_species_list.py"

echo "==> 2/8 BLAST databases"
mkdir -p "$ROOT/release/blast_db"
for name in all bacteria viruses fungi oomycetes; do
  "$MAKEBLASTDB" \
    -in "$ROOT/release/sequences/maize_pathogens_${name}.fasta" \
    -dbtype nucl \
    -out "$ROOT/release/blast_db/maize_pathogens_${name}" \
    -title "MaizePathogenDB release ${name}" >/dev/null
done

echo "==> 3/8 sequence QC"
SEQ_DIR="$ROOT/release/sequences" \
MANIFEST="$ROOT/data/sequence_manifest.tsv" \
OUT_TSV="$ROOT/data/sequence_qc_report.tsv" \
  "$PY" "$ROOT/scripts/qc_sequences.py"

echo "==> 4/8 validation"
SEQ_DIR="$ROOT/release/sequences" \
BLAST_DIR="$ROOT/release/blast_db" \
RESULT_DIR="$ROOT/docs/validation/results" \
DB_VERSION="release" \
  "$PY" "$ROOT/scripts/run_validation.py"

echo "==> 5/8 fixed-threshold validation split"
BLAST_DIR="$ROOT/release/blast_db" \
RESULT_DIR="$ROOT/docs/validation/results" \
  "$PY" "$ROOT/scripts/fixed_threshold_validation_split.py"

echo "==> 6/8 NCBI ITS comparison"
BLAST_DIR="$ROOT/release/blast_db" \
RESULT_DIR="$ROOT/docs/validation/results" \
  "$PY" "$ROOT/scripts/run_ncbi_its.py"

if [[ "${SKIP_PERFORMANCE:-0}" != "1" ]]; then
  echo "==> 7/8 performance (Usage Notes only)"
  RESULT_DIR="$ROOT/docs/validation/results" \
    "$PY" "$ROOT/scripts/run_performance.py"
else
  echo "==> 7/8 performance skipped (SKIP_PERFORMANCE=1)"
fi

QIIME_BIN="$HOME/miniconda3/envs/rachis-qiime2-2026.7/bin/qiime"
if [[ "${SKIP_UNITE:-0}" != "1" && -x "$QIIME_BIN" ]]; then
  echo "==> 8/8 UNITE comparison"
  RESULT_DIR="$ROOT/docs/validation/results" \
    QIIME="$QIIME_BIN" \
    "$PY" "$ROOT/scripts/run_unite.py"
else
  echo "==> 8/8 UNITE skipped (SKIP_UNITE=1 or QIIME2 not installed)"
fi

echo "==> done"
