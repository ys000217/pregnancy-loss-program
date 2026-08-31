# 当前计划（10x ONT + WGS 配对 CNV）

日期：2026-08-28

## 结论（先定死）

- 数据：ONT ~10x 已比对 BAM；Illumina WGS ~10x 原始 PE FASTQ；**同一人配对**。
- 问题：体质性/种系 CNV（妊娠丢失队列），不是肿瘤、不是 scRNA。
- 策略：断点（Sniffles2 / Manta）+ 剂量（CNVpytor 100 kb / 500 kb）互补；**配对重叠当主过滤器**。
- 不用：inferCNV、CNVkit 默认 exome、GATK-gCNV、FACETS、ichorCNA。
- 原始 FASTQ / ONT 原始目录只读；结果写到独立 `cnv_work/`。

## 已写到本目录的内容

| 文件 | 作用 |
|------|------|
| `README.md` | 分辨率合同、置信度规则、流水线说明 |
| `config.sh` | 集群路径和 10x 参数（上机后先改这里） |
| `scripts/00_scan_fastq.py` … `08_annotate.sh` | 逐步脚本 |
| `scripts/run_sample.sh` | 单样本全流程 |
| `scripts/run_cohort.sh` | 队列循环（作业系统需按景行再包一层） |
| `examples/manifest.tsv` | 配对表格式 |
| `env/environment.yml` | conda 环境（不含 manta） |
| `env/README.md` | 安装失败排查、可选 Manta 单独安装 |
| `scripts/check_env.sh` | 上机后检查工具是否齐 |

## 接下来做什么（按这个顺序）

1. **本机/仓库**：计划与脚本已在 `cnv/`，拷到服务器分析盘（不要拷进 `NGS_Rawdata/`）。
2. **服务器：配环境与工具** ← 你现在要做的这一步
   - 建 conda 环境（见 `env/environment.yml`）
   - 跑 `bash scripts/check_env.sh`
   - 确认 `REF_FASTA` 与 ONT BAM 的 `@SQ` 同一套 GRCh38
   - 放 tandem-repeat BED（Sniffles2 必需）
3. **编 `manifest.tsv`**：`00_scan_fastq.py` 扫 FASTQ，手工把 `ont_id` 对上 `0002C` 这类 ONT 目录名。
4. **先跑 1 个配对样本**（建议 `0002C`）通流程，再提交队列。

## 服务器上需要的工具

必装：

- python ≥3.10
- fastp
- bwa-mem2
- samtools ≥1.19
- sambamba 或 samtools markdup
- mosdepth ≥0.3
- sniffles ≥2.2
- delly（WGS 断点 SV；conda 默认）
- cnvpytor
- bcftools ≥1.19
- bedtools ≥2.31

可选：

- spectre（ONT 大片段剂量，没有就只用 CNVpytor）
- AnnotSV（LARGE_HIGH 结果注释）

参考数据（不是软件，但环境配完就要就位）：

- 与 ONT BAM 一致的 `GRCh38.fa` + bwa-mem2 索引
- `human_GRCh38_TR.bed`（或等价 tandem-repeat BED）

## 明确还没做的

- 还没有在服务器上安装/验证上述工具
- 还没有 `manifest.tsv` 实表（ONT ID ↔ FASTQ 路径）
- 还没有对照 ONT BAM header 锁定 `REF_FASTA`
- 还没有用真实样本跑通
- 景行作业脚本（`yhbatch` / 队列名）还没按站点包装

配完环境并 `check_env.sh` 全绿之后：

1. `bash scripts/discover_paths.sh 0002C` — 找 ONT BAM 真实路径和 @SQ
2. 改 `config.sh`（路径一律用 share：`SHARE_ROOT`，不要用 `~/5250028_songyang`），`mkdir -p $SHARE_ROOT/cnv_work`
3. 编 `manifest.tsv`，跑 `bash scripts/run_sample.sh 0002C`
