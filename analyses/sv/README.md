> **模块位置：** 本目录是妊娠丢失主仓库中的 **SV 分析子模块**（`analyses/sv/`）。  
> 总项目说明见仓库根目录 [README.md](../../README.md)。

# SV × 甲基化与致病元件富集

胎盘 ONT 甲基化 × 种系结构变异（SV）关联、abnormal×DMR 阴性检验、AnnotSV 致病断点相关 pair 的染色质元件富集。

## 一、主流程脚本

| 脚本 | 功能 |
|---|---|
| `sv_methylation_pipeline.py` | **核心主流程**。分 8 个 chunk（0–7） |
| `pathogenic_element_enrichment.py` | AnnotSV ACMG 4/5 断点 × 显著 pair → Roadmap 7 类 / 基因窗口富集（fig2a/b） |
| `enrichment_prom_enh.py` / `regen_fig2.py` | 转调致病元件富集脚本的兼容入口 |

```bash
cd analyses/sv
python sv_methylation_pipeline.py        # 全部 chunk 0-7
python sv_methylation_pipeline.py 2 6    # 仅指定区间
python pathogenic_element_enrichment.py  # 致病断点 × 染色质富集
```

**关键配置（主流程脚本顶部）**

- `EXCLUDE_ABNORMAL`：`True`=主分析（剔 47 abnormal，601 样本）；`False`=敏感性（648 + abnormal 协变量）
- M 值：`log2` logit（Du et al. 2010）
- 模型：`M ~ n_SV + Gestational_Week`；交互另加 `Status`、`n_SV×Status`

## 二、诊断与 abnormal×DMR 脚本

| 脚本 | 结论 |
|---|---|
| `diagnose_beta_m_dist.py` 等 | M 值适合线性回归；同方差近似成立 |
| `check_sv_index_alignment2.py` | SV 编号 0% 不一致 |
| `dmr_sv_colocalization.py` / `dmr_sv_enrichment_test.py` / `dmr_sv_quantitative.py` | abnormal 特异 SV **不能**解释 hyper-DMR（三条阴性） |

## 三、目录布局（本仓库）

```
analyses/sv/
├── README.md
├── *.py                      # 脚本
├── metadata/                 # 小临床/协变量表（可提交）
├── results/                  # 可提交的汇总与富集表
│   └── plots/                # 诊断图（勿放根目录 figures/，已被 gitignore）
└── docs/                     # 汇报稿
```

大矩阵 / 全量回归 / AnnotSV 全文 **不入库**（见根 `.gitignore`）。本地运行时大输入路径仍写在脚本顶部（如 `D:\ONT\figure2`、11GB 甲基化矩阵）。

## 四、主要结果（已入库摘要）

| 文件 | 内容 |
|---|---|
| `results/report_summary.tsv` | 主分析汇总（约 601 样本口径） |
| `results/sv_methylation_sig_pairs.tsv` | FDR 显著 SV–CpG pair |
| `results/sv_methylation_selected_sites.tsv` | 入选位点 |
| `results/interaction_significant.tsv` | 显著交互 |
| `results/pathogenic_sv_acmg45.tsv` | AnnotSV ACMG 4+5（full）目录 |
| `results/pair_7class_enrichment.tsv` | Roadmap 7 类：全部显著 vs 致病断点相关 |
| `results/pair_gene_window_enrichment.tsv` | up/dn/body 富集 |
| `results/sig_pairs_with_chromatin.tsv` | 显著 pair + 染色质类 + `patho_related` |
| `results/dmr_sv_*.tsv` | abnormal DMR×SV 共定位/定量摘要 |
| `results/plots/diagnosis_*.png` | β/M、均值方差、同方差诊断图 |
| `metadata/clinical_649.tsv` 等 | 临床与矩阵列顺序 |

### 结论速览

1. **SV × 甲基化主效应稳健**（约数千显著 pair / ~10³ 基因量级；主分析与敏感性一致）。
2. **稳健交互基因**：OR11H12、FMO2、MMEL1（对照中 SV→低甲基化，病例中效应减弱）。
3. **abnormal 全局高甲基化不能用种系 SV 解释**（DMR 共定位 / 富集 / 定量三条阴性）。
4. **致病断点相关显著 pair 的 Het-ZNF 富集**：统计上相对「致病相关非显著背景」过 FDR，但 **17/17 条 Het-ZNF pair 来自单一基因 HRNR（dn）**；去掉 HRNR 后无 Het-ZNF。宜按**基因级**解读，不宜写成「致病 pair 系统富集异染色质」。效应多为负（SV↑ → M↓），与 abnormal–PMD 高甲基化叙事方向不同，**不能**直接当作 abnormal–PMD–染色质假说的证据。

汇报材料：`docs/分析汇报_PPT.md`、`docs/分析汇报_主分析.docx`。
