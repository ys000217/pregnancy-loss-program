# Pregnancy Loss Program

妊娠丢失（pregnancy loss）研究项目的**主仓库**。

本仓库按分析模块组织；各模块目录内自带说明与脚本。原始测序数据、参考基因组、AnnotSV 注释库等大文件**不纳入**本仓库（见 `.gitignore`）。

## 仓库结构

```
pregnancy-loss-program/
├── README.md                 # 本说明（总项目）
├── .gitignore
└── analyses/
    ├── sv/                   # 结构变异（SV）× 甲基化、致病元件富集
    ├── gwas/                 # RPL GWAS 再分析 × gsMap 定位
    ├── cnv/                  # 10x 配对 ONT + WGS germline CNV
    ├── highvar_cpg/          # ONT 高变 CpG 筛选与聚类
    ├── burden/               # SNV/SV burden、AnnotSV P/LP、全基因组 SV 位点富集
    ├── external_data/        # 流产/RPL GWAS、妊娠表型 GWS、胎盘 meQTL
    └── cell_composition/     # 胎盘细胞比例：孕周趋势与三分组比较
```

## 模块一览

| 模块 | 路径 | 说明 |
|------|------|------|
| **SV 分析** | [`analyses/sv/`](analyses/sv/) | SV–甲基化关联；致病断点×染色质富集；abnormal×DMR 阴性；结果在 `results/` |
| **GWAS × gsMap** | [`analyses/gwas/`](analyses/gwas/) | 中国 RPL GWAS suggestive 再分析、λ_GC、gsMap 细胞定位尝试（阴性） |
| **CNV（配对）** | [`analyses/cnv/`](analyses/cnv/) | ~10x ONT + Illumina WGS 配对 germline CNV；流水线已写，队列结果待跑 |
| **高变 CpG** | [`analyses/highvar_cpg/`](analyses/highvar_cpg/) | abnormal 枝 + 30 例 abnormal-like；步骤 3（8/9/10 周）待跑 |
| **SNV/SV burden** | [`analyses/burden/`](analyses/burden/) | mut/Mb、P/LP 负担；PASS SV 位点富集；strict abnormal-specific = 0 |
| **外部数据** | [`analyses/external_data/`](analyses/external_data/) | 流产/RPL GWAS、Liu 2026 GWS 子集、胎盘 meQTL；外部 hit gsMap 探索 |
| **细胞组成** | [`analyses/cell_composition/`](analyses/cell_composition/) | 6–12 周细胞比例趋势；三分组比较（多不显著） |

## SV 模块快速入口

详细说明见：[`analyses/sv/README.md`](analyses/sv/README.md)

```bash
cd analyses/sv
python sv_methylation_pipeline.py        # 全部分块
python pathogenic_element_enrichment.py  # 致病断点 × Roadmap 7 类富集
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
bash scripts/remake_merge.sh 0002C
```

## 高变 CpG 模块快速入口

详细说明见：[`analyses/highvar_cpg/README.md`](analyses/highvar_cpg/README.md)

```bash
# jsub < analyses/highvar_cpg/scripts/10_GA_stratified_mainstream_highvar_cpg.sh
```

## SNV/SV burden

见 [`analyses/burden/README.md`](analyses/burden/README.md)

```bash
Rscript analyses/burden/scripts/run_all.R
```

## 外部数据

见 [`analyses/external_data/README.md`](analyses/external_data/README.md)

## 细胞组成

见 [`analyses/cell_composition/README.md`](analyses/cell_composition/README.md)

## 贡献约定

- **代码与文档**进 Git；**原始数据 / 大结果矩阵**留在本地或对象存储。
- 新增分析线时：在 `analyses/<模块名>/` 下建目录，并在本 README 的模块表中补一行。
- 勿在仓库根目录堆分析产物；汇总表/图放进对应模块的 `results/`。
