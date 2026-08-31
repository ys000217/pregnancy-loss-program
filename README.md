# Pregnancy Loss Program

妊娠丢失（pregnancy loss）研究项目的**主仓库**。

本仓库按分析模块组织；各模块目录内自带说明与脚本。原始测序数据、参考基因组、AnnotSV 注释库等大文件**不纳入**本仓库（见 `.gitignore`）。

## 仓库结构

```
pregnancy-loss-program/
├── README.md                 # 本说明（总项目）
├── .gitignore
└── analyses/
    ├── sv/                   # 结构变异（SV）相关分析
    ├── gwas/                 # RPL GWAS 再分析 × gsMap 定位
    ├── cnv/                  # 10x 配对 ONT + WGS germline CNV
    ├── highvar_cpg/          # ONT 高变 CpG 筛选与聚类
    └── external_data/        # 流产/RPL GWAS + 胎盘 meQTL 显著 hit
```

## 模块一览

| 模块 | 路径 | 说明 |
|------|------|------|
| **SV 分析** | [`analyses/sv/`](analyses/sv/) | ONT 结构变异、liftover、AnnotSV 注释、SV–甲基化关联等 |
| **GWAS × gsMap** | [`analyses/gwas/`](analyses/gwas/) | 中国 RPL GWAS suggestive 再分析、λ_GC、gsMap 细胞定位尝试 |
| **CNV（配对）** | [`analyses/cnv/`](analyses/cnv/) | ~10x ONT + Illumina WGS 配对 germline CNV；LARGE_HIGH / SHARED_SV 分层 |
| **高变 CpG** | [`analyses/highvar_cpg/`](analyses/highvar_cpg/) | 全队列分出 abnormal 枝；NC 敏感性确认 30 例 abnormal-like；主流样本仅 8/9/10 周层内高变 |
| **外部数据** | [`analyses/external_data/`](analyses/external_data/) | 流产/RPL GWAS 与胎盘 meQTL 显著 hit 表（带 PMID）；供本队列定向复现 |

后续其它模块（例如甲基化图谱、临床表型、空间组学等）可继续放在 `analyses/` 下各自子目录中。

## SV 模块快速入口

详细脚本说明见：[`analyses/sv/README.md`](analyses/sv/README.md)

主流程示例：

```bash
cd analyses/sv
python sv_methylation_pipeline.py        # 全部分块
python sv_methylation_pipeline.py 2 6    # 仅跑指定 chunk
```

## GWAS × gsMap 模块快速入口

详细说明见：[`analyses/gwas/README.md`](analyses/gwas/README.md)

```bash
cd analyses/gwas
python scripts/reanalyze_gwas_suggestive.py   # suggestive P<1e-4 + λ_GC
# WSL: bash scripts/run_gsmap_suggestive87.sh
```

## CNV 模块快速入口

详细说明见：[`analyses/cnv/README.md`](analyses/cnv/README.md)

```bash
cd analyses/cnv
source config.sh
bash scripts/remake_merge.sh 0002C          # 仅重跑合并过滤
# FORCE=1 ONLY=0002C bash scripts/submit_per_sample.sh
```

## 高变 CpG 模块快速入口

详细说明见：[`analyses/highvar_cpg/README.md`](analyses/highvar_cpg/README.md)

```bash
# 服务器：步骤 2 已完成（NC 敏感性）；步骤 3 仅 8/9/10 周层内重算
# jsub < analyses/highvar_cpg/scripts/10_GA_stratified_mainstream_highvar_cpg.sh
```

## 外部数据模块快速入口

详细说明见：[`analyses/external_data/README.md`](analyses/external_data/README.md)；文献：[`analyses/external_data/docs/REFERENCES.md`](analyses/external_data/docs/REFERENCES.md)

```bash
cd analyses/external_data
# metadata/gwas_hits.tsv   — 流产/RPL 显著 SNP（Laisk / Sonehara / 中国候选基因）
# metadata/meqtl_hits.tsv  — 胎盘 meQTL ∩ GWAS 窗（Delahaye 2018 S6）
```

## 贡献约定

- **代码与文档**进 Git；**原始数据 / 大结果矩阵**留在本地或对象存储。
- 新增分析线时：在 `analyses/<模块名>/` 下建目录，并在本 README 的模块表中补一行。
