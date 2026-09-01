# 10x paired ONT + Illumina WGS germline CNV pipeline

进度与下一步见 [`PLAN.md`](PLAN.md)。上机先配环境：`env/environment.yml`，然后 `bash scripts/check_env.sh`。

## 当前状态

**流水线与置信度规则已入库；队列 callset / LARGE_HIGH 汇总尚未作为结果表提交。** 下一步仍是服务器配环境 → `manifest.tsv` → 单样本通流程 → 队列。

妊娠丢失队列：每个样本同时有 **~10x ONT**（已比对 BAM）和 **~10x Illumina PE WGS**（原始 FASTQ）。这不是 30x WGS，也不是 scRNA inferCNV。单平台在 10x 都会漏小事件、涨假阳性；**配对交叉验证是这条流水线的主过滤器**。

## 10x 下能报什么、不能报什么

| 事件大小 | 主要证据 | 10x 下怎么处理 |
|----------|----------|----------------|
| 非整倍体 / ≥1 Mb | 两边的 depth（CNVpytor 500 kb） | 高可信，配对应一致 |
| 100 kb–1 Mb | 两边 depth（100 kb bin）± SV 断点 | 两边同向重叠 → HIGH；单边 → 降级 |
| 10–100 kb | ONT Sniffles2 为主，DELLY/Manta 辅助 | **不要用 depth 单独确认** |
| 50 bp–10 kb | 几乎只能靠 ONT Sniffles2 | WGS 10x 的 depth 和 split-read 都不够 |
| 插入（INS） | 几乎只能靠 ONT | Illumina 10x 基本看不见，不参与 CNV 目录 |

Depth 层 bin = **100 kb**（另跑 500 kb 大事件）。CBS 在 ≤5x 召回会垮；10x 用 CNVpytor 的 mean-shift，不要用 CNVkit 默认 CBS 追焦点。inferCNV / FACETS / ichorCNA / GATK-gCNV 都不作为本队列主 caller。

## 数据现状

| 平台 | 现在有什么 | 深度 | 分析入口 |
|------|------------|------|----------|
| ONT | `{sample}.merged.bam` + `.bai` | ~10x | 直接 SV + depth CNV |
| Illumina | `NGS_Rawdata/*/03.release/*_{1,2}.fq.gz` | ~10x | **先比对**，禁止改原始 FASTQ |
| 配对 | 同一人两套数据 | — | `manifest.tsv` 把 ONT ID 和 FASTQ ID 绑在一起 |

甲基化 TSV 不参与 calling，只在注释阶段叠到 CNV 区间。

## 流水线

```
manifest.tsv  (ont_id, wgs_r1, wgs_r2, sex)
        |
        +--[WGS]  fastp → BWA-MEM2 (同一套 GRCh38) → sort → markdup → BAM
        |
        +-- QC  mosdepth（qc.tsv 记 pass；低于 MIN_MEDIAN_DEPTH=8 只警告）
        |       仅当 median<3 或 WGS breadth<0.40 才中止（QC_ABORT）
        |
        +--------------------------+---------------------------+
        | ONT                      | WGS                       |
        | Sniffles2 + TR BED       | DELLY（或单独装的 Manta） |
        |   → 筛 DEL/DUP           |   → 筛 DEL/DUP            |
        | Spectre 或 CNVpytor      | CNVpytor 100kb + 500kb    |
        +--------------------------+---------------------------+
        |
        v
  统一 BED → 按样本配对合并（50% reciprocal overlap）
        |
        v
  置信度分层  LARGE_HIGH / SHARED_SV / MEDIUM / ONT_SV / ONT_DEPTH / WGS_DEPTH / MASKED
        |
        v
  AnnotSV (GRCh38) → 仅 LARGE_HIGH
```

## 置信度规则（配对的意义）

对每个样本，把 ONT 与 WGS 的 DEL/DUP 做成同一坐标系后：

| 标签 | 规则 | 含义 |
|------|------|------|
| **LARGE_HIGH** | 两边同类型，RO≥50%，长度 ≥100 kb，**不落在硬区 mask**。100 kb–1 Mb：100 kb 双边 depth；≥1 Mb：有 500 kb 时需 500 kb 支持；**>10 Mb / 非整倍体：双边 depth 即可，不强制 SV 断点** | **病例对照主表**；写入 `cnv.high.bed` |
| **SHARED_SV** | 两边都有 SV 断点，长度 &lt;100 kb | 小 SV 负担分析；**不要**和大片段混计 |
| **MEDIUM** | 跨平台同向重叠但不满足上两档 | 候选，需 IGV |
| **ONT_SV** | 仅 ONT Sniffles DEL/DUP，≥50 bp | 小 CNV 的预期来源 |
| **ONT_DEPTH** | 仅 ONT depth，≥100 kb，未 mask，≤10 Mb | 单平台，优先级低 |
| **WGS_DEPTH** | 仅 WGS depth，≥100 kb，未 mask，≤10 Mb | 单平台探索 |
| **MASKED** | 本可进 depth/LARGE，但命中硬区，或超长且证据不足，或 ≥1 Mb depth 无 500 kb 支持 | 审计用，不当阳性 |
| **DROP** | 仅 WGS 小 SV；单边 depth &lt;100 kb；**性染色体（默认）**；非 primary contig | 不报 |

硬区 BED：`ref/hard_mask.grch38.refseq.bed`（近端着丝粒短臂、着丝粒±2 Mb、1q12/9qh/16qh 等）。**缺 mask 时 merge 直接失败**（`--require-hard-mask`）。默认**丢弃性染色体**（chrX/Y）；需要时设 `KEEP_SEX_CHROM=1`。

CNVpytor 进 merge 前默认质量过滤：`Q0≤0.5`、`pN≤0.5`、`e-val1≤1e-4`（可用 `CNVPYTOR_QC=0` 关闭）。100 kb + 500 kb 双分辨率都会进 merge；≥1 Mb 且含 depth 证据时要求有 500 kb 支持。

**不要**再把旧版笼统的 `HIGH`（含上千条亚 kb SV∩SV）直接做富集。

## 目录（全部写到分析盘，不动原始数据）

```
cnv_work/
  manifest.from_ont.tsv          # 正式入口（含 ont_bam）
  wgs_bam/{sample}.markdup.bam
  qc/{sample}.{ont,wgs}.mosdepth.*
  ont_sv/{sample}.sniffles.vcf.gz
  ont_cnv/{sample}.cnvpytor.100000.tsv
  ont_cnv/{sample}.cnvpytor.500000.tsv
  wgs_sv/{sample}.wgs_sv.cnv.vcf.gz
  wgs_cnv/{sample}.cnvpytor.100000.tsv
  wgs_cnv/{sample}.cnvpytor.500000.tsv
  merged/{sample}.cnv.bed
  merged/{sample}.cnv.high.bed          # LARGE_HIGH only
  merged/{sample}.cnv.shared_sv.bed     # SHARED_SV (<100 kb both-SV)
  annot/{sample}.annotsv.tsv
```

## 怎么跑

1. 改 `config.sh` 里的 `REF_FASTA`（必须与 ONT BAM 的 `@SQ` 一致）和 `TR_BED`。
2. 生成**正式** manifest（必须含 `ont_bam`）。`00_scan_fastq.py` 只扫 WGS；配对 BAM 用 `00c`：

```bash
source config.sh
python3 scripts/00_scan_fastq.py --ngs-root "${NGS_RAWDATA}" -o "${WORKDIR}/manifest.wgs.tsv"
python3 scripts/00c_ont_wgs_coverage.py \
  --ont-root "${ONT_BAM_ROOT}" \
  --ont-extra "${ONT_RAW_MERGED_ROOT}" \
  --wgs-manifest "${WORKDIR}/manifest.wgs.tsv" \
  -o "${MANIFEST}"   # 默认 WORKDIR/manifest.from_ont.tsv，含 ont_bam
# 示例列见 examples/manifest.tsv
```

3. 单样本：

```bash
source config.sh
bash scripts/run_sample.sh 0002C
```

仅重跑合并（callset 已在盘上、改了 merge 规则之后）：

```bash
source config.sh
bash scripts/remake_merge.sh 0002C
# 同步后务必: sed -i 's/\r$//' scripts/*.sh scripts/*.py config.sh
```

4. 队列：对 `manifest.from_ont.tsv` 的 `ont_id` 循环提交（不要在 `NGS_Rawdata/` 里跑）。已有 `cnv.high.bed` 的样本会被 SKIP；改规则后用 `FORCE=1` 或先 `remake_merge.sh`。

依赖：`bwa-mem2` `samtools` `sambamba`（推荐）`fastp` `mosdepth` `sniffles` `delly` `cnvpytor` `bcftools` `bedtools`（**Manta 不要 conda 装进 cnv10x**，见 `env/README.md`）；可选 `spectre` `AnnotSV`。深度主 caller 是 **CNVpytor**；若环境有 Spectre 会额外跑（专用 mosdepth `--by 1000`，与 QC 的 100 kb 分开）。

## 和 ONT 甲基化的衔接

`0002C_methylome_final_correct.tsv` 只在有了 `merged/{sample}.cnv.bed` 之后做区间重叠。10x 深度 CNV 的断点是 bin 宽，不要用它去解释单个 CpG。
