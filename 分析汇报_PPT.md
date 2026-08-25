# 胎盘 ONT 甲基化 × 种系结构变异(SV) 关联分析 —— 汇报

> 用途：制作 PPT。每节对应一组幻灯片，图片路径已列出，可直接插图。

---

## 一、研究背景与目标

- **生物学问题**：种系结构变异(SV) 是否通过影响局部甲基化，参与妊娠结局（自然流产 SPL / 复发性流产 RPL）。
- **复现参考**：Global DNA methylation differences involving germline structural variation impact gene expression in pediatric brain tumors（Nat Commun 2025, 16:4713）。
- **数据规模**：648 例样本（含 495 case + 153 control）、95,635 个 SV、约 300 万个 CpG 位点（11GB 甲基化矩阵）。

---

## 二、分析流程（8 个 chunk）

| 步骤 | 内容 |
|---|---|
| chunk0 | gencode GTF 解析 → 蛋白编码基因坐标 |
| chunk1 | VCF → SV 携带矩阵 |
| chunk2 | 蛋白编码基因 → 上游100kb / 下游100kb / 基因体 三窗口 |
| chunk3 | 基因窗口 × SV 断点 → 每患者 SV 计数 |
| chunk4 | 全基因组回归（M ~ n_SV + 孕周） |
| chunk5 | BH-FDR 筛选（FDR<10%） |
| chunk6 | SV × case/control 交互 |
| chunk7 | 汇报图 + 汇总表 |

**核心设计**

- 因变量：M 值 = `log2(β/(1−β))`（Du et al. 2010 标准）
- 自变量：窗口内"不同 SV 数量"（剂量-反应）
- 模型：`M ~ n_SV + Gestational_Week`
- 多重检验：BH-FDR < 10%

---

## 三、方法学验证（这部分保证结论可靠）

### 3.1 M 值 vs β 值：为什么用 M 值

- β 值被压缩在 [0,1]、双峰、**方差随均值变化 10.1 倍**（异方差）；
- M 值无界、方差跨度降至 **2.5 倍**（更同方差），更适合线性回归。

**关键图**：
- `diagnosis_meanvar.png`（均值-方差对比，最核心）
- `diagnosis_beta_dist.png`（β 分布）
- `diagnosis_M_dist.png`（M 分布）

### 3.2 同方差诊断（回归假设检验）

- 残差方差随自变量 n_SV 变化极小（β 1.08 倍 / M 1.16 倍），**近似同方差**。
- 基因组膨胀系数 λ 接近 1，假阳性控制良好。

**关键图**：`diagnosis_var_vs_nsv.png`、`figures/fig6_qq.png`

### 3.3 SV 编号对齐校验

- VCF（CN1 坐标）→ liftover → BED（GRCh38 坐标）。
- 校验：BED 的 SV<k> 与 VCF 记录顺序，**SVTYPE 69,691/69,691 全部一致（0% 不一致）**。
- 27% SV 在 liftover 中丢弃（已由 inner join 正确处理）。

---

## 四、主要结果

### 4.1 主分析 vs 敏感性分析

| 指标 | 主分析（剔 47 abnormal，601 样本） | 敏感性（648 样本 + abnormal 协变量） |
|---|---|---|
| 检验位点数 | 2,228,920 | 2,233,500 |
| 总检验数 | 4,079,157 | 4,090,352 |
| **显著对（FDR<10%）** | **3,193** | **3,419** |
| **入选位点** | **2,588** | **2,711** |
| **命中基因** | **1,033** | **1,068** |
| p 阈值 | 7.83×10⁻⁵ | 8.35×10⁻⁵ |
| λ | 1.061 | 1.066 |
| 高/低甲基化 | 46.6% / 53.4% | 44.9% / 55.1% |
| 显著交互 | 5 | 3 |

**结论**：两种口径结果高度一致，主效应稳健。

**关键图**：
- `figures/fig1_manhattan.png`（曼哈顿图）
- `figures/fig2_window_bar.png`（窗口分布：up/dn/body）
- `figures/fig3_top_genes.png`（Top 基因）

### 4.2 显著交互（SV × case/control）

| 基因 | 区域 | 主分析 int_p | 敏感性 int_p |
|---|---|---|---|
| OR11H12 | 上游启动子 | 5.1×10⁻⁷ | 4.3×10⁻⁷ |
| FMO2 | 上游启动子 | 6.8×10⁻⁵ | 1.4×10⁻⁵ |
| MMEL1 | 下游 3'调控 | 1.4×10⁻⁴ | 6.6×10⁻⁵ |
| SLC22A25 | 基因体/上游 | 1.2×10⁻⁴ | 消失 |

- 三个基因（OR11H12、FMO2、MMEL1）在两种分析中都稳健显著。
- 效应方向一致：对照组 SV 使甲基化降低，病例组该效应被显著减弱。

---

## 五、abnormal 样本成因探索（DMR × SV 联合）

### 5.1 背景

- 47 个 abnormal 样本（全部为 case，SPL 22 / RPL 25），表现为**全局高甲基化**。
- 目的：寻找导致 abnormal 离群的 SV。
- 数据：modkit DMR（abnormal vs control，分孕周 g8/g9/g10，spl_g9 作废）。

### 5.2 三条证据链（全部阴性）

| 检验 | 结果 |
|---|---|
| ① 单 SV 富集 | 2,971 个 p<0.05，**少于随机期望 3,483**；无任何 SV 在 abnormal 中高渗透 |
| ② 空间共定位 | abnormal 富集 SV **不**比背景更靠近 hyper-DMR（方向甚至相反） |
| ③ 定量因果 | 携带 SV **不**预测 DMR 甲基化（中位 p=0.50，BH-FDR=0） |

### 5.3 顺带发现

- 500 个最强 DMR 中，**77%（385 个）在甲基化矩阵里无 CpG 覆盖**——提示 DMR（全基因组 ONT）与矩阵（3M CpG）不是同一套 CpG 集合。

### 5.4 结论

> **SV 无法解释 abnormal 样本的全局高甲基化。** 三条独立证据交叉印证，该方向应终止。

---

## 六、总体结论

1. SV × 甲基化关联分析顺利完成，主效应稳健（2,588 / 2,711 位点，1,033 / 1,068 基因）。
2. 方法学严谨：M 值（log2 logit）标准、同方差、λ≈1、SV 编号 0% 错配。
3. 显著交互锁定 3 个稳健基因（OR11H12 / FMO2 / MMEL1）。
4. abnormal 的全局高甲基化**不由 SV 驱动**，需转向技术性排查或其他遗传/表观机制。

---

## 附：图片清单

| 图片 | 用途 |
|---|---|
| `diagnosis_meanvar.png` | M 值 vs β 值的方差稳定性（方法学核心） |
| `diagnosis_beta_dist.png` / `diagnosis_M_dist.png` | 分布对比 |
| `diagnosis_var_vs_nsv.png` | 同方差诊断 |
| `figures/fig1_manhattan.png` | 曼哈顿图 |
| `figures/fig2_window_bar.png` | 窗口效应分布 |
| `figures/fig3_top_genes.png` | Top 基因 |
| `figures/fig6_qq.png` | QQ 图（λ 膨胀） |
