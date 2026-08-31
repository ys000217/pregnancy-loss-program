# SNV / SV 突变负荷与全基因�?SV 位点富集

无偏 SNV（编码区非同�?mut/Mb）、SV 计数负荷、AnnotSV P/LP SV（含罕见子集），以及全基因组 PASS SV 逐位�?Fisher 富集�?

大文件（VCF、AnnotSV、ANNOVAR）不入库；默认从 `D:/ONT/figure2` �?`D:/ONT` 读取（可用环境变量覆盖）�?

## 输入（本地，不入库）

| 文件 | 用�?|
|------|------|
| `figure2/sample_phenotype_648.tsv` | 648 样本 Condition（abnormal / normal / control�?|
| `figure2/WGS_ONT_Intersection_648samples.vcf.gz` | 无偏 SNV callset |
| `figure2/S956.snp.annovar.hg38_multianno.txt.gz` | SNV ANNOVAR 注释 |
| `clinical_649.GRCh38.correct.vcf` | 无偏 SV callset |
| `clinical_649.GRCh38.annotsv.tsv` | AnnotSV（ACMG、人�?AF�?|

环境变量�?

- `BURDEN_ROOT`：表型与 SNV 数据根目录（默认 `D:/ONT/figure2`�?
- `BURDEN_OUT` / `BURDEN_PLOT`：输出表/图目录（默认本模�?`tables/`、`plots/`�?

## 主要结果（摘要）

- **SNV mut/Mb**：按 Condition / case–control 比较（见 `tables/group_comparison_stats.tsv`�?
- **SV 总负�?/ P/LP / 罕见 P/LP**：每样本计数与组间检�?
- **全基因组 PASS SV 位点**�?5,388 位点；三组两两独�?Fisher + 各自 FDR
- **Abnormal 特异**（严格）：`fdr(ab vs ctrl)<0.05` �?`ab_rate>ctrl` �?`ab_rate>norm` �?`fdr(ab vs norm)<0.05` �?**0 位点**
- ab vs ctrl FDR&lt;0.05 �?63 个，多数�?`case_vs_control` 模式（normal 也升高），不�?abnormal 特异

完整逐位点表 `sv_locus_enrichment_all_pass.tsv`（约 14 MB）留在本�?`figure2/burden_analysis/tables/`，不入库�?

## 运行

```bash
# 需能访问上述大文件
Rscript analyses/burden/scripts/run_all.R
# 或分�?
python analyses/burden/scripts/compute_burden.py
python analyses/burden/scripts/compute_sv_locus_enrichment.py
Rscript analyses/burden/scripts/analyze_and_plot.R
```

工作副本也可保留�?`figure2/burden_analysis/`（`figure*/` �?`.gitignore` 忽略）�?
