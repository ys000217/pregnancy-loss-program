# 参考文献（外部数据模块）

本模块所有 hit 均有文献出处；PMID / DOI 也写在各 TSV 行内。

## GWAS / 候选基因关联

1. **Laisk T** et al. The genetic architecture of sporadic and multiple consecutive miscarriage. *Nature Communications*. 2020;11:5980.  
   PMID: [33239672](https://pubmed.ncbi.nlm.nih.gov/33239672/) · DOI: [10.1038/s41467-020-19742-5](https://doi.org/10.1038/s41467-020-19742-5)  
   → 欧洲血统 sporadic miscarriage GWS：`rs146350366`（FGF9）；multiple consecutive miscarriage GWS：`rs7859844`、`rs143445068`、`rs183453668`（Firth 后）。`rs138993181` 仅 Pmeta GWS、Firth 未过，表中标为 P2。

2. **Sonehara K** et al. Common and rare genetic variants predisposing females to unexplained recurrent pregnancy loss. *Nature Communications*. 2024;15:5744.  
   PMID: [39019884](https://pubmed.ncbi.nlm.nih.gov/39019884/) · DOI: [10.1038/s41467-024-49993-5](https://doi.org/10.1038/s41467-024-49993-5)  
   → 日本 uRPL GWAS lead：`rs9263738`（MHC）；HLA-C\*12:02–B\*52:01–DRB1\*15:02 保护性单倍型。

3. **Fan Q** et al. The synergic effects of CTLA-4/Foxp3-related genotypes and chromosomal aberrations on the risk of recurrent spontaneous abortion among a Chinese Han population. *Journal of Human Genetics*. 2018;63:579–587.  
   PMID: [29476189](https://pubmed.ncbi.nlm.nih.gov/29476189/) · DOI: [10.1038/s10038-018-0414-2](https://doi.org/10.1038/s10038-018-0414-2)  
   → 中国汉族候选基因：`rs231775`、`rs3087243`（CTLA4）；`rs2232365`、`rs2232368`（FOXP3）。优先级 P2（非全基因组扫描）。

4. **Liu** et al. Genetic polymorphism of IL-17 influences susceptibility to recurrent pregnancy loss in a Chinese population. *Medicine (Baltimore)*. 2024;103(23):e38333.  
   PMID: [38847697](https://pubmed.ncbi.nlm.nih.gov/38847697/) · DOI: [10.1097/MD.0000000000038333](https://doi.org/10.1097/MD.0000000000038333)  
   → 中国小样本候选基因：`rs2275913`（IL17A）、`rs763780`（IL17F）。优先级 P2。

## 胎盘 meQTL（SNP–CpG）

5. **Delahaye F** et al. Genetic variants influence on the placenta regulatory landscape. *PLoS Genetics*. 2018;14(11):e1007785.  
   PMID: [30452450](https://pubmed.ncbi.nlm.nih.gov/30452450/) · DOI: [10.1371/journal.pgen.1007785](https://doi.org/10.1371/journal.pgen.1007785)  
   → S6 Table：permutation 通过的 4,342 个胎盘 cis-mQTL。本模块提取与流产/RPL GWAS lead ±1 Mb（hg19）重叠的子集入 `meqtl_hits.tsv`。

## 登记但未入库全表的资源（后续可扩展）

| 资源 | 文献 / 入口 | 用途 |
|------|-------------|------|
| INMA 胎盘 cis-mQTL 浏览器 | Cilleros-Portet et al. *Nat Commun* 2025; DOI [10.1038/s41467-025-57760-3](https://doi.org/10.1038/s41467-025-57760-3)；https://irlab.shinyapps.io/shiny_mqtl_placenta/ | 按坐标查询更多胎盘 meQTL |
| PACE 胎盘 mQTL meta | 预印本 DOI [10.64898/2026.07.07.26357471](https://doi.org/10.64898/2026.07.07.26357471)；https://pace-placenta-mqtl.streamlit.app/ | 更大样本 lead mQTL |
| FinnGen O15_ABORT_SPONTAN | https://risteys.finngen.fi/endpoints/O15_ABORT_SPONTAN | 当前公开页未列稳定 GWS hit 表；有 sumstats 时可再筛 |

## 坐标说明

- **统一构建：GRCh38**（`genome_build=GRCh38`）。
- `gwas_hits.tsv`：主坐标 `chrom`/`pos` 来自 NCBI dbSNP GRCh38 `chrpos`；`pos_hg19_legacy` 仅审计用。
- `meqtl_hits.tsv`：Delahaye S6 原文为 hg19；已用 UCSC `liftOver` + `hg19ToHg38.over.chain.gz` 全部替换为主坐标 `snp_pos` / `cpg_*`；`*_hg19_legacy` 保留原文。
- 复现时请直接使用 GRCh38 列，勿再使用 legacy 列。
