# 外部数据模块

> **模块位置：** `analyses/external_data/`（妊娠丢失主仓库子模块）  
> 总项目说明见仓库根目录 [README.md](../../README.md)。

本模块收集**流产 / RPL 相关 GWAS 与胎盘 meQTL（SNP–CpG）的显著结果**，供本队列定向复现。每条 hit 带 PMID/DOI；文献清单见 [`docs/REFERENCES.md`](docs/REFERENCES.md)。

## 当前交付

| 文件 | 内容 | n |
|------|------|---|
| [`metadata/gwas_hits.tsv`](metadata/gwas_hits.tsv) | 流产/RPL 显著或候选 SNP（**GRCh38**） | 12 |
| [`metadata/meqtl_hits.tsv`](metadata/meqtl_hits.tsv) | 胎盘 meQTL ∩ GWAS 窗（**GRCh38**，自 hg19 liftOver） | 84 |
| [`docs/REFERENCES.md`](docs/REFERENCES.md) | 参考文献与 PMID | — |
| [`docs/DATA_PATHS.md`](docs/DATA_PATHS.md) | 原始补充表本地路径 | — |

### 优先级

| 级别 | 含义 | 示例 |
|------|------|------|
| **P0** | 东亚正式 GWAS / 与东亚 MHC 直接相关的胎盘 meQTL | Sonehara `rs9263738`；MHC/FGF9 窗内 meQTL |
| **P1** | 欧洲大样本流产 GWS + 对应窗 meQTL | Laisk GWS 位点 |
| **P2** | 候选基因或 Firth 未过信号 | 中国 CTLA4/FOXP3/IL17；Laisk `rs138993181` |

## 复现用法（概要）

1. **GWAS**：用 `gwas_hits.tsv` 的 `rsid` / 坐标在本队列做定向 case–control（注意祖先与效应等位基因）。  
2. **meQTL**：用 `meqtl_hits.tsv` 的 GRCh38 `snp_id`–`cpg_id` 坐标；`kgp*` 需先映射 rs，再在 ONT beta 上复测效应方向。  
3. 大文件（完整 S6）在 `data/`（gitignore），不入库。

```bash
cd analyses/external_data
# 查看
# head metadata/gwas_hits.tsv
# head metadata/meqtl_hits.tsv
```

重建 meQTL 窗过滤：

```bash
# 需本机 R + readxl；原始 S6 路径见 docs/DATA_PATHS.md
Rscript scripts/rebuild_meqtl_window_hits.R
```

## 目录

```
analyses/external_data/
├── README.md
├── docs/
│   ├── REFERENCES.md
│   └── DATA_PATHS.md
├── metadata/
│   ├── gwas_hits.tsv
│   └── meqtl_hits.tsv
├── scripts/
│   └── rebuild_meqtl_window_hits.R
├── data/          # gitignore：原始补充表 / 全量 mQTL
└── results/
```

## 不入库

- 完整 GWAS/meQTL sumstats、原始基因型、甲基化矩阵  
- `data/` 下 Delahaye S6 等补充表全文
