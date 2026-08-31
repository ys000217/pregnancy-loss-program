# 外部数据 — 本地路径与来源

> 大文件不进 Git。小 hit 表在 `metadata/`。

## 路径约定

| 角色 | 位置 | 说明 |
|------|------|------|
| 下载缓存 | `analyses/external_data/data/raw/` | 期刊补充表 |
| 处理后全表 | `analyses/external_data/data/processed/` | 如全量 4342 mQTL TSV |
| 可提交清单 | `metadata/` | `gwas_hits.tsv` / `meqtl_hits.tsv` |

## 数据集登记表

| ID | 来源 | 数据类型 | 用途 | 本地路径 | 状态 |
|----|------|----------|------|----------|------|
| LAISK2020 | PMID 33239672 | GWAS GWS hits（手录自主文） | 定向 SNP 复现 | `metadata/gwas_hits.tsv` | 已登记 |
| SONEHARA2024 | PMID 39019884 | 日本 uRPL GWAS lead | 东亚优先复现 | 同上 | 已登记 |
| FAN2018 | PMID 29476189 | 中国 CTLA4/FOXP3 候选 SNP | P2 | 同上 | 已登记 |
| LIU2024 | PMID 38847697 | 中国 IL17 候选 SNP | P2 | 同上 | 已登记 |
| DELAHAYE2018_S6 | PMID 30452450 | 胎盘 cis-mQTL S6（4342） | 窗过滤 → meQTL hits | `data/raw/Delahaye2018_PLoSGenet_S6_mQTL.xlsx` | 已下载；子集已入库 |
| DELAHAYE2018_FULL_TSV | 同上 | S6 转 TSV | 本地查询 | `data/processed/delahaye2018_placenta_mqtl_4342.tsv` | 本地 |

### 下载备注

- Delahaye S6：PLOS Genetics 补充文件 `journal.pgen.1007785.s018`（表名 `T6_mQTLs`）。  
- Laisk 位点：主文统计量 + dbSNP 坐标；完整 Supplementary Data 可按需再下。
- liftOver：`D:\gsMap\tools\liftOver` + `D:\gsMap\tools\hg19ToHg38.over.chain.gz`（UCSC）。meQTL 已全部转为 GRCh38。

## 备注

- 新增数据集：先改本表与 `docs/REFERENCES.md`，再更新 `metadata/`。
