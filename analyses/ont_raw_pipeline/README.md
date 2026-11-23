# ONT 原始数据（RawData）处理 pipeline

**这是原始数据（RawData）处理 pipeline，不是下游甲基化 / SV / CNV 分析。**

本模块只做 ONT trio 原始 fast5 到合并 modBAM 的三步：fast5 → pod5 → 按批次 dorado modBAM → 合并为一个 BAM。甲基化、结构变异、CNV 在 `analyses/sv/`、`analyses/cnv/`、`analyses/highvar_cpg/`、`analyses/pmd/` 等模块，是另外的分析线，不要把它们写成这条 pipeline 的下一步。旧 SHARE 队列的甲基化 / SV 下游流程同样是另一条分析线。

脚本在 `scripts/`，文件名与集群家目录一致，正文未改。`run_dorado_modbam_gpu_by_sample.sh` 会调用 `/data/home/4250160/run_dorado.sh`，所以 `run_dorado.sh` 一并入库。

## 来源

2026-11-23 从集群家目录 `/data/home/4250160/` 原样拷贝。公开仓库只保留这些工作副本；脚本里是集群绝对路径（`/data/home/4250160/` 与实验室共享盘）。不要往本模块加凭证、密钥或样本名单。

| 脚本 | 字节 | md5 | 集群 mtime |
|------|------|-----|------------|
| `convert_fast5_to_pod5_by_sample_ont_trio.sh` | 1457 | `fc19dfc7f7c614e987706fb9b6e345a9` | 2026-11-13 |
| `run_dorado_modbam_gpu_by_sample.sh` | 1371 | `da975d6bd60c8506b53d637a1548ad8a` | 2025-12-16 |
| `run_dorado.sh` | 334 | `96f2654c00750370ab1e10913910bde9` | 2025-10-14 |
| `run_merge_modbam_by_sample.sh` | 4078 | `88ebdc3b3723b1a357a9f3d29eedaeeb` | 2026-11-23 |

## 三步

按这个顺序跑。三份主脚本都是 jsub 作业脚本，提交前必须先 `export SAMPLE`。

### 1. Convert：fast5 → pod5

脚本 `scripts/convert_fast5_to_pod5_by_sample_ont_trio.sh`。

- 读 `ONT_trio/02.ONT/01.RawData`。完整路径已写在脚本里：`/data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_Lab_shared_data/ONT_trio/02.ONT/01.RawData`。按 `Release-*/Data-*/${SAMPLE}` 找该样本，再找其中的 `*fast5_pass`。
- conda 环境 `modkit`，`pod5 convert fast5`，56 线程。
- 写出 `~/pod5/$SAMPLE/<batch>/`（脚本中为 `/data/home/4250160/pod5/${SAMPLE}/<batch>/`）。`<batch>` 是 `*fast5_pass` 上一级目录名。

### 2. Dorado：pod5 → 按批次 modBAM

脚本 `scripts/run_dorado_modbam_gpu_by_sample.sh`，内部调用 `/data/home/4250160/run_dorado.sh`（仓库副本为 `scripts/run_dorado.sh`）。

- dorado 1.1.1（`/data/home/4250160/dorado-1.1.1-linux-x64/bin/dorado`），`run_dorado.sh` 用 glibc 2.35 的动态链接器把二进制拉起来。
- 参考基因组 GRCh38.p14：`/data/home/4250160/refseq/grch38.p14/GRCh38.p14.fasta`。
- 模型：`dna_r10.4.1_e8.2_400bps_hac@v5.2.0`，修饰碱基模型 `…_hac@v5.2.0_5mC_5hmC@v2`。`--device cuda:0`。conda 环境 `dorado`，GPU 队列。
- 读 `~/pod5/$SAMPLE` 下每个 `*.pod5`。写出 `~/modbam/$SAMPLE/<batch>_mods.bam/`（脚本中为 `/data/home/4250160/modbam/${SAMPLE}/<batch>_mods.bam/`）。dorado 在这个目录里再套一层 run 目录，BAM 在 `bam_pass/*.bam`。

### 3. Merge：按批次 BAM → 一个样本一个 BAM

脚本 `scripts/run_merge_modbam_by_sample.sh`。

- 写出 `~/modbam/$SAMPLE/$SAMPLE.merged.bam`（脚本中为 `/data/home/4250160/modbam/$SAMPLE/$SAMPLE.merged.bam`），并建索引。
- micromamba 环境 `bam_merge`（可执行文件 `/data/home/4250160/bin/micromamba`）。`samtools` 经 `micromamba run -n bam_merge` 调用。
- 合并成功后不删除输入，也不搬走原始批次文件。日志写在 `~/modbam/merge_logs/$SAMPLE.log`。
- 若合并 BAM 与索引（`.bai` 或 `.csi`）都已存在且非空，直接跳过并退出 0。
- 要求该样本下每一个 `*_mods.bam` 批次目录的 mtime 都早于 60 分钟（`MIN_AGE_MINUTES=60`）。有任何一个不够旧就报错退出。
- 每个批次目录里至少要有一个 `bam_pass/*.bam`，否则报错退出。
- `NUM_BAMS=1` 时 `cp -p` 该 BAM 再 `samtools index`；多于一个时 `samtools merge -f -c -p --write-index`。
- 写出后做非空与索引检查，再跑 `samtools quickcheck`。

## 集群用法

在集群上对家目录里的脚本提交（与入库副本内容相同）。`SAMPLE` 必须先导出，脚本里是 `SAMPLE="${SAMPLE:?}"`，没设会直接失败。

```bash
export SAMPLE=样本ID
jsub < /data/home/4250160/convert_fast5_to_pod5_by_sample_ont_trio.sh
jsub < /data/home/4250160/run_dorado_modbam_gpu_by_sample.sh
jsub < /data/home/4250160/run_merge_modbam_by_sample.sh
```

先等 convert 的 pod5 落盘，再交 dorado；dorado 的批次目录都超过 60 分钟、且每个批次都有 `bam_pass/*.bam` 之后，再交 merge。

## 已知行为

下面这些是现有脚本的行为，入库时没有改脚本。

- 调度器显示 “Successfully completed”，或脚本自己 echo 了 “Finished”，都不能当作成功。要看输出文件在不在、全不全。
- convert 包装脚本不检查 `pod5 convert` 的退出码。某个批次转换失败，脚本仍可能走到最后的 Finished。
- dorado 包装脚本不检查 dorado 的退出码。循环结束后，只要在 `~/modbam/$SAMPLE` 下找到任意一个 `*/bam_pass/*.bam`，就会执行 `rm -rf ~/pod5/$SAMPLE`。部分批次失败时，已经产出的 BAM 也会触发这次删除，该样本的 pod5 会整目录被清掉。
- 不要用旧的扫描脚本 `run_gpu_batches_scan_pod5.sh`。它在 `~/modbam/{sample}` 已经是目录时会跳过，而且任务之间可能重叠。
- 不要对整个目录 jsub 旧的 unify merge。它只做一层 `find`，对不上这批数据里 dorado 套出来的 `bam_pass` 布局。
- 不要对本批数据用 `merge_one_sample.sh`。它的 `ROOT` 写死在共享盘的 `modbam_Dorodo`，不是 `~/modbam`。
- 旧 SHARE 队列的甲基化 / SV 下游流程是另一条分析线，不要写成这条原始数据 pipeline 的下一步。

## 目录

```
analyses/ont_raw_pipeline/
├── README.md
└── scripts/
    ├── convert_fast5_to_pod5_by_sample_ont_trio.sh
    ├── run_dorado_modbam_gpu_by_sample.sh
    ├── run_dorado.sh
    └── run_merge_modbam_by_sample.sh
```
