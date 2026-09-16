# MaizePathogenDB Validation

This directory contains the corrected marker-clean validation query sets and
results. Query sequences were trimmed to the same marker space as the
corrected reference database.

## Corrected validation numbers

- Positive queries: 675 (bacteria 32, viruses 30, fungi 537, oomycetes 76)
- Retained negative queries: 487 (bacteria 93, viruses 100, fungi 246, oomycetes 48)
- Cross-database queries: 560 (bacteria 26, viruses 20, fungi 434, oomycetes 80)
- External retrieval: species 76.1% (514/675); genus 96.7% (653/675)
- Species classification: sensitivity 61.6%; specificity 93.8%; precision
  93.3%; F1 74.2; balanced accuracy 77.7%
- Cross-database consistency: species 77.7% (435/560); genus 95.9% (537/560)
- Fixed-threshold validation half: species sensitivity 57.9%, specificity
  94.6%, F1 71.5; genus sensitivity 92.6%, specificity 91.9%, F1 93.3

## Results

- Full summary: `validation/results/SUMMARY.md`
- Fixed split: `validation/results/fixed_threshold_validation_split.json` and
  `validation/results/FIXED_THRESHOLD_VALIDATION_SPLIT.md`
- NCBI ITS comparison: `validation/results/ncbi_its_comparison.json` and
  `validation/results/NCBI_ITS_COMPARISON.md`
- UNITE comparison: `validation/results/unite_comparison.json` and
  `validation/results/UNITE_COMPARISON.md`
- Performance: `validation/results/performance.json` and
  `validation/results/PERFORMANCE.md`

## NCBI-nt

The previously used local NCBI-nt snapshot is not available in the corrected
release environment. The NCBI-nt head-to-head was therefore not rerun, and old
NCBI-nt values are not mixed with the marker-clean MPDB results.
