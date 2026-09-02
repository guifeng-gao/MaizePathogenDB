# QIIME2 2026.7 Environment (macOS arm64)

The V4 UNITE comparison uses QIIME 2 2026.7. On Apple Silicon, the official
distribution ships `osx-64` environments that run under Rosetta 2. The full
`rachis-qiime2` YAML install was attempted first, but `packages.qiime2.org`
was too slow from the current network (about 20 KB/s per connection), so the
environment was built as a minimal equivalent: the same pinned QIIME2 2026.7
packages required for `qiime feature-classifier classify-sklearn`
(`qiime2`, `q2-types`, `q2-feature-classifier`, `q2cli`), with the rest of
the dependencies resolved by conda-forge/bioconda.

## Install commands used on 2026-08-25

```bash
# 1. Install Miniconda (arm64) if absent
bash Miniconda3-latest-MacOSX-arm64.sh -b -p "$HOME/miniconda3"

# 2. Accept Anaconda channel terms of service (required by conda >= 26)
"$HOME/miniconda3/bin/conda" tos accept \
  --override-channels --channel https://repo.anaconda.com/pkgs/main
"$HOME/miniconda3/bin/conda" tos accept \
  --override-channels --channel https://repo.anaconda.com/pkgs/r

# 3. Create the official QIIME2 2026.7 environment (Rosetta/osx-64)
export CONDA_SUBDIR=osx-64
"$HOME/miniconda3/bin/conda" create -y \
  -n rachis-qiime2-2026.7 \
  -c https://packages.qiime2.org/qiime2/2026.7/qiime2/released \
  -c conda-forge -c bioconda \
  python=3.12 qiime2=2026.7.0 q2-feature-classifier=2026.7.0 \
  q2-types=2026.7.0 q2cli=2026.7.0

# 4. Verify
"$HOME/miniconda3/envs/rachis-qiime2-2026.7/bin/qiime" info
```

The full official YAML remains the reference distribution:
`https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-osx-64-conda.yml`.

## Why 2026.7 and osx-64

- 2026.7 is the latest QIIME2 release (2026-08-25).
- The full `rachis-qiime2` distribution for 2026.7 publishes `linux-64` and
  `osx-64` environment files. The minimal environment above installs only the
  packages needed for classification and pins the same `2026.7.0` versions.
- The pinned classifier artifact
  `unite_ver2025-02-19_dynamic_fungi-Q2-2026.4.qza` was trained with
  scikit-learn 1.7.1, which matches the 2026.7 environment (`scikit-learn=1.7.1`).

## Verification status

- Install and `qiime info`: verified during the V4 UNITE run (2026-08-25).
- Classifier artifact read by QIIME2 2026.7: verified during the V4 UNITE run.
