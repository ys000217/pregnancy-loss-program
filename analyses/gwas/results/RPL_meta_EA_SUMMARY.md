# 东亚 RPL 整合分析 × gsMap 结果（中国 + 日本）

> 生成日期：2026-08-17。整合用户自己的中国 RPL GWAS 与东亚（日本）RPL GWAS，跑 gsMap 检验是否出现多基因信号与细胞类型定位。

## 一句话结论

**把中国 RPL（N≈608）与日本 RPL（Sonehara 2024，N=1,728/24,315）做 meta 整合后，全基因组 mean χ² = 0.965（仍 < 1，零假设水平），gsMap 无任何细胞类型显著（best p_cauchy = 0.362）。整合不能把流产变成多基因性状。**

---

## 1. 整合方案

| 数据集 | 群体 | 表型 | N(case/ctrl) | 基因组版本 | 变异数 |
|---|---|---|---|---|---|
| combine_gwas_result_v1 | 中国 | 复发性流产 RPL | ~608 总 | hg38 | 364,183 |
| Sonehara 2024 (hum0197) | 日本 | 复发性流产 RPL | 1,728 / 24,353 | hg19 | 8,717,430 |

- 用户数据 hg38 → liftOver 到 hg19（复用上一 session 的 `user_hg19.bed`，303,404 个变异映射成功）。
- 按 `chr:hg19pos` + 等位基因（含 strand-flip 补齐）匹配，**sample-size 加权固定效应 meta**（权重 = √N_eff）。
- Sonehara N_eff = 4×1728×24353/(1728+24353) ≈ 6,454；用户 N_eff ≈ 608（总 N 近似，case/ctrl 拆分未知）。
- 效应等位基因约定：用户 A1=ALT，Sonehara Allele2=ALT，z=BETA/SE；跨研究对齐符号。
- 输出映射到 HapMap3 rsID（LDSC weights 1,242,190 个），剔除 strand-ambiguous（A/T、C/G）。

## 2. Meta 结果（金标准 = mean χ²）

| 指标 | 整合后 |
|---|---|
| HapMap3 SNP 总数 | **993,327** |
| 两研究共同（matched） | 48,276 |
| 仅日本（Sonehara-only） | 923,090 |
| 仅中国（user-only） | 21,961 |
| **mean χ²** | **0.9653** |
| λ_GC | 1.02 |
| max χ² | 21.85 |
| **fuel（mean χ² − 1）** | **−0.035** |

判读：mean χ² ≈ 0.965 < 1，与单个数据集（用户 1.022、Sonehara 0.966）一致，**零假设水平**。S-LDSC 燃料 = mean χ²−1 ≈ 0（甚至略负），无多基因信号可归因。对比身高/智力 mean χ² ≈ 1.4–1.8。

## 3. OLA1 位点：小样本假阳性实锤

用户 top hits（hg38 chr2:174.09–174.19Mb = hg19 chr2:174.96–175.06Mb）在独立日本数据中**不重复、甚至方向相反**：

| 用户 top hit (hg38) | 用户 z | 日本 z（同变异） | meta 后 z |
|---|---|---|---|
| chr2:174095142 | −4.95 (p=7.3e-7) | +0.64 | **−2.07** |
| chr2:174193638 | −4.95 (p=7.4e-7) | +0.63 | −2.05 |
| chr2:174134460 | −4.90 (p=9.6e-7) | +0.46 | −2.02 |
| chr2:174094826 | −4.72 (p=2.4e-6) | +0.67 | −1.99 |
| chr2:174190584 | −4.68 (p=2.8e-6) | +0.64 | −1.92 |

meta 后 z 从 ≈−4.9 稀释到 ≈−2.0（p≈0.05，未过 GW 阈值），日本数据呈弱反向——强烈提示这些小样本位点是 winner's curse 假阳性。

## 4. gsMap 结果（250K HapMap3 子集，243,851 共同 SNP）

| 步骤 | 结果 |
|---|---|
| ④ spatial LDSC | 121,767 spots；per-spot median z² = 0.76（纯零）；min p = 9.5e-5（> Bonferroni 4.1e-7，无一 spot 显著） |
| ⑤ cauchy（细胞类型） | **best p_cauchy = 0.362（Spinal cord）**，其余均 > 0.82 |

**无任何细胞类型显著**。对照：阳性对照 IQ 的 Brain p_cauchy = 2.45e-24；上一 session 用户单数据 best p_cauchy = 0.084（同样阴性）。

## 5. 文件位置

- 整合脚本：`D:\ONT\figure3\gwas\meta_rpl_ea.py`
- 验证脚本：`D:\ONT\figure3\gwas\verify_meta.py`
- 子集脚本：`D:\ONT\figure3\gwas\subset_meta_250k.py`
- 全量整合 sumstats：`D:\gsMap\RPL_GWAS\RPL_meta_EA.sumstats.gz`（993,327 SNP）
- gsMap 用子集：`D:\gsMap\RPL_GWAS\RPL_meta_EA_250k.sumstats.gz`（250,000 SNP）
- meta 报告：`D:\gsMap\RPL_GWAS\RPL_meta_EA.report.txt`
- gsMap 输出：`.../E16.5_E1S1.MOSTA/spatial_ldsc/E16.5_E1S1.MOSTA_RPL_meta_EA.csv.gz`、`.../cauchy_combination/E16.5_E1S1.MOSTA_RPL_meta_EA.Cauchy.csv.gz`

## 6. 结论

1. **整合（中国+日本 RPL）不改变结论**：mean χ² = 0.965，仍无多基因信号。
2. **gsMap 无法定位**：best p_cauchy = 0.362，无细胞类型显著。
3. **用户 OLA1 top hit 是假阳性**：在独立日本数据中反向、meta 后消失。
4. 这再次确认 RPL 是**寡基因/罕见变异**性状（少数大效应罕见位点），不是 S-LDSC 能吃到的海量小效应多基因性状。
