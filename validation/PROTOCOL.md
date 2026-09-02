# MaizePathogenDB Validation Protocol

**Applies to**: MaizePathogenDB release (frozen content)  
**Freeze date**: 2026-08-26  
**Status**: in force  

## 1. Purpose

本协议定义 MaizePathogenDB release 的验证口径、查询集、阈值、外部数据库版本、输出格式和变更控制。论文中出现的所有验证数字必须来自本协议定义的运行，禁止混用历史 benchmark。

## 2. Definitions

- **Query TaxID**: 查询序列来源物种的 NCBI Taxonomy ID；对 catalog 内物种使用 `species_list.tsv` 中的 TaxID。
- **Species-level correct**: top-1 hit 的 TaxID 与 query TaxID 完全一致。
- **Genus-level correct**: top-1 hit 的 genus 与 query 的 genus 一致（以固定 NCBI taxonomy dump 判断）。
- **No-hit**: 无 BLAST hit 或 evalue 超过 1e-5，记为错误/漏检，不剔除。
- **Independent query**: 与 release 参考序列不 100% 全长度一致（global identity < 100% 或 qcov < 100%）。
- **TAXID_PENDING species**: 14 个无 NCBI TaxID 的 catalog 物种，不进入需要 TaxID 比对的验证，单独报告为“未验证目录条目”。

## 3. Reference build

- `release/sequences/maize_pathogens_all.fasta`：6,133 条。
- 分类库：bacteria 402 / viruses 430 / fungi 4,509 / oomycetes 792。
- `release/blast_db/`：由 BLAST+ 2.17.0 构建。
- `release/maize_pathogens_taxonomy_sintax.fasta`：6,133 条 SINTAX。
- `data/species_list.tsv`：225 个 catalog 条目；201 个有序列；14 个 TAXID_PENDING。
- 所有 SHA-256 见 `CHECKSUMS.sha256`。

## 4. Pinned external resources

| Resource | Version / file | Notes |
|---|---|---|
| NCBI Taxonomy | `taxdump.tar.gz` downloaded 2026-08-21 | 所有 lineage 和 genus 判断以此为准 |
| NCBI ITS_eukaryote_sequences | downloaded 2026-08-21, md5 存档 | 真菌/卵菌对照 |
| NCBI ITS_RefSeq_Fungi | downloaded 2026-08-21, md5 存档 | 真菌对照 |
| UNITE | v10.0 dynamic, release 2025-02-19, `unite_ver2025-02-19_dynamic_fungi-Q2-2026.4.qza` | 仅真菌，不含卵菌 |
| SILVA | SSU Ref NR99 138.2 (release 2024-07-11) | 下载后记录 md5；16S 对照 |
| NCBI nt | 必须使用本地快照并记录下载日期与 `update_blastdb.pl --showall` 输出 | 不用未存档的 web BLAST 作为论文数字 |

## 5. Query sets

### Q1 Internal (completeness only)
- 全部 6,133 条参考序列，对分类内 BLAST 库 top-1 自命中。
- 仅作完整性检查，不称为准确率。

### Q2 Independent external positives（主验证）
- 来源：NCBI GenBank，`[PDAT]` 在 2020-01-01 至 2026-08-25 之间。
- 每个有序列 catalog 物种最多 10 条 marker 序列；排除与参考序列 100% 全长一致的记录。
- 目标：真菌 ≥500 条，细菌/病毒/卵菌 ≥30 条；若类别物种数不足则使用全部可用物种。
- 记录每个物种的查询日期、检索式、命中数和排除数。
- 输出：`docs/validation/query_sets/external_positives.fasta` + `external_positives_meta.tsv`。

### Q3 Negatives
- 500 条非玉米病原序列：100 细菌 / 100 病毒 / 250 真菌 / 50 卵菌。
- 条件：TaxID 不在 catalog；与 Q2 和参考序列无 100% 全长一致；固定后不再变动。
- 输出：`docs/validation/query_sets/negatives.fasta` + `negatives_meta.tsv`。

### Q4 Cross-database consistency
- 对每个有序列 catalog 物种，从已存档的 SILVA/UNITE/NCBI RefSeq 查询集和
  NCBI 独立序列中取非 100% 一致的序列，每物种最多 5 条；卵菌通过 NCBI
  ITS 补齐。
- 直接使用固定资源：`Figshare/docs/validation/external/silva_unite_cross_queries.fasta`（SILVA/UNITE 查询）与 `Figshare/docs/validation/external/refseq_queries.fasta`（病毒 RefSeq 查询）。
- SILVA 138.2 本地完整数据库为待办项；下载并记录 md5 后，Q4 细菌部分可替换为直接从 SILVA 138.2 提取。
- 输出：`docs/validation/query_sets/cross_db_queries.fasta` + `cross_db_queries_meta.tsv`。

### Q5 Performance (informational)
- 固定 ASV 查询集，规模 5/10/20/40/60/80/100/150/200/260（或等价划分）。
- 记录 OS、CPU、内存、BLAST 版本。

## 6. Validation tasks

### Internal completeness
- BLAST：top-1，`blastn -evalue 1e-5`。
- 报告：每类正确数 / 总数；不进入论文准确率声明。

### Primer / region coverage
- 引物集合固定：细菌 27F/1492R、338F/806R、515F/806R、338F/1392R；真菌 ITS1F/ITS4、ITS5/ITS4、ITS1/ITS4、ITS86F/ITS4、fITS7/ITS4、ITS1F/ITS2、ITS3/ITS4、ITS9mun/ITS4ngs。
- 规则：结合位点 ≤2 mismatches 视为覆盖；区域覆盖按参考注释判断。

### Independent external retrieval against MPDB
- 查询集：Q2。
- 方法：MPDB 本地 BLAST top-1。
- 报告口径：
  - Retrieval accuracy：top-1，无阈值，种级 / 属级。
  - Classification accuracy：种级 `pident>=99, qcovs>=90`；属级 `pident>=95, qcovs>=70`。
  - 四类分别报告，不合并掩盖低分类准确率。

### Head-to-head against NCBI-nt
- 条件：存在固定日期下载的 NCBI-nt 本地快照，或存档全部 web BLAST RID/XML。
- 当前状态：已完成（NCBI-nt 快照 2026-08-23）。
- 公平对照口径：排除 415 条“top-1 命中即查询自身 accession”的自匹配后，
  使用剩余 260 条独立查询。结果：MPDB 种级检索 71.5% vs NCBI-nt 63.1%；
  属级检索 96.2% vs 91.9%；种级分类 66.2% vs 63.1%；属级分类 95.4% vs
  91.9%；负样本假阳性 0.6%（3/500）。详见 `results/NCBI_NT_COMPARISON.md`。

### Classification benchmark against general databases
- 阳性：Q2；阴性：Q3。
- 方法：MPDB 本地 BLAST、NCBI ITS_eukaryote / ITS_RefSeq 固定库、UNITE
  QIIME2 分类器（置信度 >= 0.7，仅真菌）；NCBI-nt 本地快照已完成。
- 指标：sensitivity、specificity、F1、balanced accuracy，按类报告。
- 严格档 `pident>=99.5, qcovs>=99` 只作为 strain-confirmation 说明，不作为主结果。
- 禁止再使用“从数据库自身抽样的 573 阳性 / 150 条 benchmark”作为准确率。

### Cross-database consistency
- 查询集：Q4；方法：top-1 种级 / 属级。
- 逐库报告：SILVA / UNITE / RefSeq。

### Performance (informational)
- 方法：Q5；指标：wall time、peak RSS、per-ASV time。
- 结果只能放 Usage Notes，不参与技术验证结论，也不是 Data Descriptor 必需内容。
- 若使用 260 样本 ASV 数据，仅作为固定性能查询集，不得包含病原检出、丰度比较或生态学结论。
- 固定 ASV 查询集已使用；结果只放 Usage Notes。

## 7. Output and reporting

- 输出目录：`docs/validation/results/`。
- 每次运行必须保存：
  - 查询集 FASTA + metadata；
  - BLAST raw TSV / XML；
  - per-query JSON；
  - summary JSON 和 summary MD；
  - 运行日志（日期、seed、脚本 hash、外部库版本）。
- 论文中每个验证数字必须能追溯到上述文件。

## 8. Change control

- 任何对 release 数据、查询集、阈值或外部库版本的修改，都需要重跑验证并更新结果。
- 协议修改记录在本文件。

## 9. 当前执行状态

| 任务 | 状态 |
|---|---|
| Internal completeness | 已完成（release） |
| Primer coverage | 已完成（release） |
| External retrieval against MPDB | 已完成（release） |
| NCBI-nt head-to-head | 已完成（2026-08-23 快照；`results/NCBI_NT_COMPARISON.md`） |
| Classification benchmark (MPDB + fixed NCBI ITS) | 已完成（release） |
| UNITE comparison | 已完成（QIIME2 2026.7，仅真菌；结果见 `results/UNITE_COMPARISON.md`） |
| Cross-database consistency | 已完成（release） |
| Performance | 已完成（固定 ASV 查询集，Usage Notes only） |

## 10. Threshold exploration (diagnostic only)

- 官方与推荐口径统一：种级 `pident>=99, qcovs>=90`；属级 `pident>=95, qcovs>=70`。
- 固定推荐阈值的无泄漏分层结果：`docs/validation/results/fixed_threshold_validation_split.json` 与 `FIXED_THRESHOLD_VALIDATION_SPLIT.md`。

## 11. First run checklist

- [ ] 确认 release FASTA/SINTAX/BLAST 库与 `CHECKSUMS.sha256` 一致。
- [ ] 下载 SILVA 138.2 并记录 md5。
- [ ] 确认 UNITE qza、NCBI ITS_eukaryote、RefSeq、taxdump 版本与第 4 节一致。
- [ ] 构建 Q2/Q3/Q4 固定查询集。
- [ ] 下载并固定 NCBI nt 本地快照（或确认可存档 web BLAST RID）。
- [ ] 运行完整验证，写入 `docs/validation/results/`。
- [ ] 更新 `docs/validation/README.md`。
- [ ] 将主验证数字填入论文 Technical Validation。
