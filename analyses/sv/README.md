> **模块位置：** 本目录是妊娠丢失主仓库中的 **SV 分析子模块**（nalyses/sv/）。  
> 总项目说明见仓库根目录 [README.md](../../README.md)。
# 项目脚本索引（README）

> 本文件梳理本次分析涉及的脚本、数据与产出，便于后续查阅与复现。

## 一、主流程脚本

| 脚本 | 功能 |
|---|---|
| `sv_methylation_pipeline.py` | **核心主流程**。胎盘 ONT 甲基化 × 种系结构变异(SV) 关联分析，分 8 个 chunk（0-7），可单独或按区间运行 |

**运行方式**

```bash
python sv_methylation_pipeline.py        # 全部 chunk 0-7
python sv_methylation_pipeline.py 2      # 只跑 chunk 2
python sv_methylation_pipeline.py 2 6    # 跑 chunk 2-6
```

**关键配置（脚本顶部）**

- `EXCLUDE_ABNORMAL`：`True`=主分析（剔除 47 个 abnormal，601 样本）；`False`=敏感性分析（保留 648 样本 + abnormal 协变量）
- M 值转换用 `log2`（`mtransform` 函数），严格按 Du et al. 2010
- 模型：`M ~ n_SV + Gestational_Week`（主效应）；交互另加 `Status`、`n_SV×Status`

**chunk 说明**

| chunk | 输入 → 输出 |
|---|---|
| 0 | gencode GTF → `genes_grch38.tsv`（一次性，可跳过） |
| 1 | VCF → `sv_carriers.tsv`（SV 携带矩阵） |
| 2 | 基因表 → `genes_windows.tsv`（蛋白编码基因 4 窗口） |
| 3 | 窗口 × SV 断点 → `gene_window_patients.tsv`（每患者 SV 计数） |
| 4 | 11GB 矩阵全基因组回归 → `sv_methylation_results.tsv` + `sv_methylation_pvals.npy` |
| 5 | BH-FDR 筛选 → `sv_methylation_sig_pairs.tsv` + `sv_methylation_selected_sites.tsv` |
| 6 | SV × case/control 交互 → `interaction_results.tsv` + `interaction_significant.tsv` |
| 7 | 汇报图 + 汇总 → `figures/*.png` + `report_summary.tsv` |

## 二、本次新建的诊断/分析脚本

| 脚本 | 做了什么 | 结论 |
|---|---|---|
| `diagnose_beta_m_dist.py` | β vs M 分布 + 均值-方差图 | M 值适合线性回归 |
| `check_invariant_cpg.py` | 量化不变量 CpG 占比 | 0%（无退化风险） |
| `diagnose_var_vs_nsv.py` | 方差 vs 自变量 n_SV 的同方差诊断 | 近似同方差 |
| `check_sv_index_alignment.py` | SV 编号对齐校验（第1版，方法有误，废弃） | — |
| `check_sv_index_alignment2.py` | SV 编号对齐校验（正确版，按 SVTYPE） | 0% 不一致 |
| `summary_abnormal_dmr.py` | DMR 文件规模/方向摘要 | — |
| `dmr_sv_colocalization.py` | 阶段1：abnormal 特异 SV × hyper-DMR 共定位 | 非特异（65% 共定位，无区分度） |
| `dmr_sv_enrichment_test.py` | 受控富集检验（富集 SV vs 背景 SV 的 DMR 邻近度） | 阴性（方向甚至相反） |
| `dmr_sv_quantitative.py` | 阶段2：读 11GB 矩阵，验证 SV 携带是否预测 DMR 甲基化 | 阴性（FDR=0） |

## 三、关键数据文件

| 文件 | 说明 |
|---|---|
| `clinical_649.tsv` | 临床表（Group1=case/control，Group2=SPL/RPL，Group4=abnormal，孕周） |
| `matrix_covariates.tsv` | 样本协变量（Status/FID/Age/Gestational_Week/细胞比例/PC） |
| `matrix_fid.txt` | 648 个样本 ID（矩阵列顺序） |
| `sv_carriers.tsv` | 每个 SV 的携带者列表 |
| `E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed` | SV 断点（GRCh38） |
| `E:\甲基化数据矩阵\甲基化CpG位点矩阵.txt` | 11GB 甲基化矩阵（位点×样本） |
| `figure2\abnormal_*_{spl,rpl}_g{8,9,10}\segments_genome.bed` | modkit DMR 结果（spl_g9 作废） |

## 四、主要产出

| 目录/文件 | 说明 |
|---|---|
| `results_main\` | 主分析结果（601 样本）+ `figures\` 6 张图 |
| `figures\` | 敏感性分析（648 样本）6 张图 |
| `report_summary.tsv` | 汇总表（当前为敏感性分析） |
| `sv_methylation_results.tsv` | 全量检验结果 |
| `sv_methylation_sig_pairs.tsv` / `selected_sites.tsv` | 显著结果 |
| `interaction_significant.tsv` | 显著交互 |
| `diagnosis_*.png` | 诊断图（分布/均值方差/同方差） |

## 五、最终结论速览

1. **SV × 甲基化主分析**：完成，主分析与敏感性结果高度一致（2,588 vs 2,711 位点）。
2. **M 值正确性**：log2 logit 是公认标准，经诊断确认适合线性回归。
3. **SV 编号对齐**：正确，0% 不一致。
4. **abnormal 成因（DMR × SV 联合）**：三条证据全部阴性，SV 无法解释 abnormal 的全局高甲基化。
