# 当前计划（10x ONT + WGS 配对 CNV）

日期：2026-09-04（相对 2026-08-28 / 2026-09-01 复查后同步）

## 结论（先定死）

- 数据：ONT ~10x 已比对 BAM；Illumina WGS ~10x 原始 PE FASTQ；**同一人配对**。
- 问题：体质性/种系 CNV（妊娠丢失队列），不是肿瘤、不是 scRNA。
- 策略：断点（Sniffles2 / DELLY）+ 剂量（CNVpytor 100 kb / 500 kb）互补；**配对重叠当主过滤器**。
- 不用：inferCNV、CNVkit 默认 exome、GATK-gCNV、FACETS、ichorCNA。
- 原始 FASTQ / ONT 原始目录只读；结果写到独立 `cnv_work/`（`WORKDIR`）。

## 规则摘要（与 README 一致）

- **LARGE_HIGH**：必须 **ONT depth + WGS depth**、RO≥50%、≥100 kb、硬区 mask 外；≥1 Mb 有 500 kb 时须 **两侧都有 500 kb**；**>10 Mb 双边 depth，不强制 SV**。
- 跨平台但未双 depth（如 ONT SV+WGS depth）→ **MEDIUM**；ONT 独有 → **ONT_SV / ONT_DEPTH**（保留，不丢）。
- CNVpytor 进 merge 前默认 Q0/pN/e-val1 过滤；缺 hard mask → merge 失败。
- 默认丢 X/Y（`KEEP_SEX_CHROM=1` 才保留）。
- 队列负担用 `09_cluster_breakpoints`（跨样本断点聚类），不要按每人 BED 行计位点。

## 正式入口（按这个顺序）

```text
00_scan_fastq.py
  → 00c_ont_wgs_coverage.py  →  ${WORKDIR}/manifest.from_ont.tsv  （含 ont_bam）
prepare_ref_index.sh         →  bwa-mem2 index + samtools faidx（队列提交前只做一次）
ONLY=0002C submit / run_sample.sh
  → … merge → AnnotSV → touch done/${sample}.done
全队列 submit_per_sample.sh（跳过已有 done/*.done；FORCE=1 重跑）
09_cluster_breakpoints.sh    →  cohort/LARGE_HIGH 位点表
```

改 merge 规则后：对已有 callset 用 `remake_merge.sh`（也会刷新 `.done`）。

## 本目录内容

| 文件 | 作用 |
|------|------|
| `README.md` | 分辨率合同、置信度、怎么跑 |
| `config.sh` | 集群路径与 10x 参数 |
| `scripts/00_scan_fastq.py` … `00c_*.py` | 扫 FASTQ、配对 ONT BAM → 正式 manifest |
| `scripts/prepare_ref_index.sh` | **队列前**建参考索引（单样本不再建） |
| `scripts/prepare_cnvpytor_ref.sh` | CNVpytor GC/conf 一次性准备 |
| `scripts/01` … `08` | 比对 → QC → SV/depth → merge → AnnotSV |
| `scripts/09_cluster_breakpoints.*` | 跨样本 LARGE_HIGH 断点聚类 |
| `scripts/run_sample.sh` | 单样本全流程；成功写 `done/${id}.done` |
| `scripts/remake_merge.sh` | 只重跑 merge+annotate |
| `scripts/submit_per_sample.sh` | jsub 提交；按 `.done` SKIP |
| `examples/manifest.tsv` | 含 `ont_bam` 的示例表头 |
| `ref/hard_mask.grch38.refseq.bed` | LARGE_HIGH 硬区 |
| `env/` | conda（不含 Manta） |

## 服务器上机清单

1. 拷贝本模块到分析盘；`sed -i 's/\r$//' scripts/*.sh scripts/*.py config.sh`
2. `conda env create -f env/environment.yml && conda activate cnv10x`
3. `bash scripts/check_env.sh`（含参考索引是否存在）
4. 确认 `REF_FASTA` 与 ONT BAM `@SQ` 一致（RefSeq `NC_*`）
5. `bash scripts/prepare_ref_index.sh`（若 check 报 MISS）
6. `bash scripts/prepare_cnvpytor_ref.sh`（若尚未做 GC）
7. 生成 `manifest.from_ont.tsv`（见 README）
8. 先 `ONLY=0002C bash scripts/submit_per_sample.sh` 或 `run_sample.sh`，再全队列
9. 全员 merge 后：`bash scripts/09_cluster_breakpoints.sh`

## 已关闭的复查项（2026-09-01 / 最新版）

Spectre 1 kb、500 kb 进 merge、CNVpytor QC、>10 Mb 不强制 SV、hard mask 必填、mask union、primary/性染色体、`fixmate -m`、`.done`、BWA 预索引；**LARGE_HIGH 强制双边 depth、≥1 Mb 强制双边 500 kb**；AnnotSV `-svtBEDcol 4`。染色体 NC_*→1–22 映射按实测可用暂不改。
