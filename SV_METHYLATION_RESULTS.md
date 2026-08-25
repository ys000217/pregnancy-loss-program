# SV × 甲基化 分析结果（最终版：窄窗口 + 携带率5-95% + SV数量）

> 复现：*Global DNA methylation differences involving germline structural variation impact gene expression in pediatric brain tumors* (Nat Commun 2025)。
> 基因集：仅 protein_coding（20,065 个）。
> 本轮改进：① 去掉 ±1Mb 宽窗口，只留 上游100kb/下游100kb/基因体；② 窗口携带率限制在 5%–95%；③ 自变量改为「窗口内 SV 数量」（剂量-反应），替代二元「有无 SV」。

---

## 1. 输入与规模

| 项 | 值 |
|---|---|
| 样本 | 648（495 case / 153 control） |
| CpG 位点 | 3,059,870（22 条常染色体，GRCh38） |
| 非 TRA SV | 95,604 个，96,452 断点 |
| 分析基因 | protein_coding 20,065（有可回归窗口者 14,093） |
| 有效窗口 | 24,231（3 窗口 × 携带率 5–95%） |
| 总检验数 | 4,090,352 |

## 2. 方法

- 每个 CpG 映射到最近 protein_coding 基因（TSS ≤1Mb，超出的位点丢弃，丢弃率 27%）。
- 自变量 = 该患者在窗口内的**不同 SV 个数**（计数去重，DEL 的 L/R 两端只算一次）。
- 模型：`M值 ~ n_SV + Age + 孕周 + 4类细胞比例`。
- 窗口：上游100kb / 下游100kb / 基因体；携带率 5%–95%。
- BH-FDR<10%。

## 3. 主结果

- 显著（位点,窗口）对：**4,334**（FDR<10%，p 阈值 1.06e-4）。
- **入选位点：3,497**（0.11%），覆盖 **1,248 个基因**。
- 按窗口：下游100kb 1,696；上游100kb 1,531；基因体 1,107。
- λ_GC = 1.18。
- 效应方向：SV 携带者低甲基化 73% / 高甲基化 27%。

### Top 基因

| 基因 | 入选位点 |
|---|---|
| ZFP37 | 271 |
| MMEL1 | 122 |
| CWH43 | 82 |
| DEFB115 | 80 |
| PTPN20 | 53 |
| ZNF565 | 47 |
| ZDHHC11 | 41 |
| KCNJ18 | 38 |
| TMEM242 | 36 |
| VAMP3 / RRS1 / POLE4 | 32–30 |

## 4. SV × case/control 交互（校正饱和后）

- 对 4,334 显著对跑 `M ~ n_SV + Status + n_SV×Status + 协变量`，全部可估。
- **仅 2 个交互显著（FDR<10%）**（此前含 ±1Mb 宽窗口时为 116 个，多为饱和假象）：

| 基因 | 窗口 | SV主效应 | 交互效应 | 交互p | case携带 | control携带 |
|---|---|---|---|---|---|---|
| OR11H12 | 上游100kb | −1.16 | +0.96 | 3.3e-7 | 59 | 12 |
| FMO2 | 上游100kb | −3.51 | +2.71 | 3.3e-5 | 85 | 15 |

- 解读：这两个位点「SV 导致的低甲基化」效应，在病例中的衰减/方向与对照显著不同（交互为正 = case 中效应减弱）。携带者几十人，比之前饱和窗口可信，但仍仅 2 个，且均为上游窗口。

## 5. 产物（`D:\ONT\`）

- `sv_ewas_results.tsv` / `sv_ewas_sig_pairs.tsv` / `sv_ewas_selected_sites.tsv`
- `interaction_results.tsv` / `interaction_significant.tsv`（2 行）
- `figures\fig1~fig6.png`、`report_summary.tsv`

## 6. 局限

1. 甲基化比例（无覆盖数）→ M 值线性模型；2. 低覆盖位点方差大；3. n_SV 仍聚合不同 SV（未按 SVTYPE 分层）；4. 47 例离群样本已纳入；5. 交互仅 2 个，样本量受限。
