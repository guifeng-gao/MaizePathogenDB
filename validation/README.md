# MaizePathogenDB Validation

## Protocol and results

- Protocol: `docs/validation/PROTOCOL.md`
- Final results: `docs/validation/results/SUMMARY.md`
- Fixed-threshold validation split:
  `docs/validation/results/FIXED_THRESHOLD_VALIDATION_SPLIT.md`
- NCBI-nt head-to-head:
  `docs/validation/results/NCBI_NT_COMPARISON.md`
- NCBI ITS comparison:
  `docs/validation/results/NCBI_ITS_COMPARISON.md` and
  `docs/validation/results/ncbi_its_comparison.json`
- UNITE comparison:
  `docs/validation/results/UNITE_COMPARISON.md` and
  `docs/validation/results/unite_comparison.json`
- Performance (Usage Notes only): `docs/validation/results/PERFORMANCE.md`
  and `docs/validation/results/performance.json`
- Query sets: `docs/validation/query_sets/`

## Final validation numbers

- Independent positives: 675 (bacteria 32, viruses 30, fungi 537,
  oomycetes 76); negatives: 500; cross-database queries: 565
  (bacteria 30, viruses 20, fungi 435, oomycetes 80).
- External retrieval: species 78.8%, genus 97.0%; species classification
  (99/90) 60.9%, genus classification (95/70) 92.4%.
- Classification benchmark: sensitivity 60.9%, specificity 93.8%,
  F1 73.6, balanced accuracy 77.3%.
- Cross-database consistency: species retrieval 75.8%, genus retrieval
  95.0%.
- NCBI-nt head-to-head (260 non-self queries): species retrieval
  71.5% vs 63.1%; genus retrieval 96.2% vs 91.9%.
- Fixed-threshold validation half: species 58.2/94.6/71.8/76.4; genus
  91.6/91.9/92.7/91.7.
