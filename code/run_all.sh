#!/usr/bin/env bash
set -euo pipefail

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")/scripts" && pwd)"
REPO_ROOT="$(cd "$SCRIPTS/../.." && pwd)"
PY="${PYTHON:-python3}"
MAKEBLASTDB="${MAKEBLASTDB:-makeblastdb}"
BLASTN="${BLASTN:-blastn}"

export MPDB_ROOT="$REPO_ROOT"
export BLASTN
export SEQ_DIR="$REPO_ROOT/sequences"
export BLAST_DIR="$REPO_ROOT/blast_db"
export QUERY_DIR="$REPO_ROOT/validation/query_sets"
export RESULT_DIR="$REPO_ROOT/validation/results"

echo "==> 1/5 BLAST databases"
for name in all bacteria viruses fungi oomycetes; do
  "$MAKEBLASTDB" \
    -in "$SEQ_DIR/maize_pathogens_${name}.fasta" \
    -dbtype nucl \
    -out "$BLAST_DIR/maize_pathogens_${name}" \
    -title "MaizePathogenDB marker-clean ${name}" >/dev/null
done

echo "==> 2/5 core validation"
"$PY" "$SCRIPTS/run_validation.py"

echo "==> 3/5 fixed-threshold validation split"
"$PY" "$SCRIPTS/fixed_threshold_validation_split.py"

if [[ "${SKIP_PERFORMANCE:-0}" != "1" ]]; then
  echo "==> 4/5 performance"
  "$PY" "$SCRIPTS/run_performance.py"
else
  echo "==> 4/5 performance skipped"
fi

if [[ "${SKIP_EXTERNAL:-1}" != "1" ]]; then
  echo "==> 5/5 optional external database comparisons"
  "$PY" "$SCRIPTS/run_ncbi_its.py"
  QIIME_BIN="${QIIME:-$HOME/miniconda3/envs/rachis-qiime2-2026.7/bin/qiime}"
  if [[ -x "$QIIME_BIN" ]]; then
    QIIME="$QIIME_BIN" "$PY" "$SCRIPTS/run_unite.py"
  fi
else
  echo "==> 5/5 external comparisons skipped (set SKIP_EXTERNAL=0 after providing fixed databases)"
fi

echo "==> done"
