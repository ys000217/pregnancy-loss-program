# RPL × gsMap 项目交接文档

> 生成日期：2026-08-17。本文件供另一个 agent 无缝接管本项目。
> 语言约定：正文中文，关键术语/数值保留英文，避免歧义。

---

## 1. 项目目标（一句话）

调查**复发性流产（RPL / 反复流产）**的遗传信号，评估是否能用 **gsMap**（GNN + S-LDSC 的空间转录组 × GWAS 定位工具）把 RPL 的遗传信号定位到细胞类型。

**核心已证实的结论（详见 §3）**：RPL 不是强多基因性状，所有公开的流产/RPL GWAS 都缺乏多基因信号，**gsMap/S-LDSC 在结构上无法定位这类性状**。这不是工具或内存问题，是性状本身的遗传架构问题。

---

## 2. 参考基因组（关键前提）

- 本项目 T2T 参考基因组是 **CN1v0.8.1（T2T-CN1）**，不是 GRCh38。（见 memory `ont-reference-genome-is-cn1`）
- 用户自己的 GWAS 数据是 **hg38**（PLINK `.glm.logistic.hybrid`），中国人群，N≈608。
- gsMap 的 1000G EUR LD 参考面板是 **hg19/GRCh37**，两者之间需 **liftOver**。

---

## 3. 已完成的探索与核心数字

### 3.1 自己的 RPL 数据（`D:\ONT\figure3\gwas\`）

| 文件 | 变体数 | mean χ² | GW 位点(p<5e-8) | 说明 |
|---|---|---|---|---|
| `combine_gwas_result_v1.Status.glm.logistic.hybrid` | 364,183 | 1.022 | 0 | gsMap 实际用的输入 |
| `gwas_result_RPL_v1.Status.glm.logistic.hybrid` | 649,536 | — | 2 | |
| `gwas_result_RPL_v3_simplified.Status.glm.logistic.hybrid` | 649,536 | — | 2 | |

- **所有文件的 top hits 都聚集在 chr2:174.07–174.22 Mb（hg38）= OLA1 基因座**。
- 具体：combine top = chr2:174095142 (p=7.3e-7)；RPL_v1/RPL_v3 top = chr2:174185614 (p=1.24e-8 / 3.1e-9)。**所有 z 为负（保护性，OR≈0.11–0.20）**。
- OLA1 座 hg38 chr2:174.07–174.22Mb ↔ hg19 chr2:174,914,728–175,114,728。

### 3.2 gsMap 直接跑 RPL 数据 → 纯阴性

- 用 `RPL_combine_gsmap.sumstats.gz`（204,641 SNPs，mean χ²=1.080）跑 gsMap quick_mode steps ④⑤。
- 结果：**无任何细胞类型显著**，best p_cauchy = 0.084，per-spot mean z² = 1.028（纯零）。
- 阳性对照 IQ：Brain p_cauchy = 2.45e-24（证明 pipeline 本身正常）。

### 3.3 OLA1 抬升实验 → 证明「单个强位点推不动 gsMap」

- 把 OLA1 座 62 个 SNP 抬到 Z=+8（χ²=64，因 chisq_max=80 过滤阈值，Z 不能超过 8.94）。
- 结果**仍为阴性**（best p_cauchy=0.166）。
- 诊断根因：
  - 62 个抬升 SNP 里只有 **21 个**在 HapMap3 参考面板内（其余被丢弃）；
  - 「燃料」= mean χ²−1 只从 0.028 → 0.048，而 IQ 的燃料 ≈ 0.6，差一个数量级；
  - 要凑到 IQ 级别燃料，需 **~508 个 SNP 同时 χ²=80**。
- **结论：S-LDSC 只吃「海量中等效应 SNP」的多基因信号，不认「少数强位点」。** 这是 gsMap 的结构性限制。

### 3.4 公开 GWAS 的多基因信号实测（本 session 新增，金标准 = mean χ²）

| 数据集 | 性状 | 病例/对照 | **mean χ²** | λ_GC | 最大χ² | 多基因信号 |
|---|---|---|---|---|---|---|
| GCST90018786 (Sakaue 2021) | Abortion — UKB+FinnGen+BBJ | 7,069/250,492 | **0.965** | 0.92 | 27.9 | ❌ 无，连 1 个 GW 位点都没有 |
| GCST011888 (Laisk 2020) | Sporadic miscarriage | 49,996/174,109 | **0.926** | 0.97 | 30.7 | ❌ 无，1 个罕见位点 |
| GCST011887 (Laisk 2020) | **Multiple consecutive miscarriage（复发性）** | 750/150,215 | **1.023** | 0.57 | >200 | ❌ 无，3 个罕见位点 |
| FinnGen O15_ABORT_SPONTAN | 自然流产 | ~200K | 1.028 | — | — | ❌ 无 |
| 日本 RPL (Sonehara 2024) | 复发性流产 | 1,728 | 0.966 | — | — | ❌ 无 |
| **自己的数据** | RPL | ~608 | 1.022 | — | — | ❌ 无 |

**判读**：mean χ² 全部 ≈ 0.93–1.03 = 零假设水平。S-LDSC 燃料 = mean χ²−1 ≈ 0。对比身高/智力 mean χ² ≈ 1.4–1.8。

### 3.5 生物样本库普查结论

- **UK Biobank**：有公开流产 GWAS（Laisk 2020、Sakaue 2021 均含 UKB），汇总统计可在 GWAS Catalog 下载，但都无多基因信号。
- **中国版 biobank（CKB，中国慢性病前瞻性研究）**：**没有任何公开的流产/RPL GWAS**，只有流行病学研究（流产与 CVD/糖尿病/死亡）。
- **东亚唯一**：日本 BBJ（含在 GCST90018786）+ 日本 RPL（Sonehara）。
- **23andMe 2024 转祖先 GWAS**（334K 散发 + 52K 复发病例，10+1 个位点）：汇总统计受限（需数据访问协议），且 52K 复发病例也只找到 1 个 GW 位点，仍非强多基因。
- Laisk 2020 SNP 遗传力 **仅 ~5%（SE 0.4）**，散发性 1 个 / 复发性 3 个位点，全是**罕见(MAF 0.5–6.4%)大效应(OR 1.4–3.8)**变异 = 寡基因架构。

---

## 4. 数据 / 文件位置

### 自己的数据（Windows `D:\` = WSL `/mnt/d/`）
- GWAS 原始文件：`D:\ONT\figure3\gwas\`（3 个 `.glm.logistic.hybrid`）
- 分析脚本（本 session 创建的）：
  - `investigate_top_hits.py` — 找 top hits
  - `build_boosted_sumstats.py` — OLA1 抬升 sumstats
  - `diagnose_boost.py` / `diagnose_ola1_gene.py` — 抬升诊断
  - `calc_mean_chisq_abortion.py` — 流式 mean χ²（读 beta/se，只适用 Sakaue 列序）
  - `calc_mean_chisq_pval.py` — 逐行 norm.ppf（慢，已弃用）
  - `calc_mean_chisq_fast.py` — **向量化版本（用这个）**，读 p_value 列，`scipy.special.ndtri`
  - `search_gwascat.py` — GWAS Catalog API 检索

### gsMap（WSL）
- gsMap 代码：`/home/administrator/gsmap_env/lib/python3.10/site-packages/gsMap/`
- 资源目录：`/mnt/d/gsMap/gsMap_resource/`（LD 面板、weights_hm3_no_hla、baseline、gtf、quick_mode 预计算）
- 演示工作目录：`/mnt/d/gsMap/example_quick_mode/Mouse_Embryo/E16.5_E1S1.MOSTA/`
  - 结果：`spatial_ldsc/*.csv.gz`（121,767 spots）、`cauchy_combination/*.csv.gz`
- RPL 输入 sumstats：`/mnt/d/gsMap/RPL_GWAS/`
  - `RPL_combine_gsmap.sumstats.gz`（204,641 SNPs）
  - `RPL_boost_OLA1.sumstats.gz`（62 SNP 抬升版）
- 拷贝出的交付物：`/mnt/d/ONT/figure3/gwas/gsmap_RPL_results/`

### 本 session 下载的公开汇总统计（大文件，可复现）
- `/mnt/d/ONT/figure3/gwas/laisk_recurrent.h.tsv.gz`（934MB，GCST011887 复发性流产）
- `/mnt/d/ONT/figure3/gwas/laisk_sporadic.h.tsv.gz`（1.2GB，GCST011888 散发性流产）
- （GCST90018786 未落盘，是流式算的）

---

## 5. 环境与操作要点（重要的坑）

- **运行环境是 WSL（Ubuntu）**，不是 Windows Git Bash。默认 Git Bash 里 `python3` 不存在。
- **调 Python 用**：`wsl.exe -e bash -c "python3 /mnt/d/..."`（脚本路径必须 `/mnt/d/`，且要包在 `bash -c "..."` 里，否则 MSYS 会做路径转换导致 `C:/Program Files/Git/mnt/...` 报错）。
- **gsMap 专用 python**（有 scipy 1.15.3 / pandas）：`/home/administrator/gsmap_env/bin/python`
- **网络**：curl/wget 都能通 EBI/GWAS Catalog；但 `WebFetch` 被企业策略拦截。EBI 大文件（>900MB）下载常被服务器断连（SSL unexpected eof），**用 `wget -c --tries=20 --timeout=60 --waitretry=5` 断点续传**。
- **GWAS Catalog harmonised 文件有两种列序**：
  - Sakaue 2021 文件：`chromosome base_pair_location ... beta standard_error ...`（beta 在第 4/5 列）
  - Laisk 2020 文件：`hm_variant_id hm_rsid ... p_value ... beta standard_error`（hm_ 列在前，beta/se 在第 22/23 列且为 NA）
  - **Laisk 文件只有 `p_value` 有值**，必须从 p 反推 χ² = (Φ⁻¹(1−p/2))²，用 `calc_mean_chisq_fast.py`（自动按列名找 `p_value`）。

---

## 6. gsMap 关键技术事实（做任何新分析前先读）

- **Pipeline 6 步**：① find_latent_representations(GNN) → ② latent_to_gene → ③ generate_ldscore → ④ run_spatial_ldsc(S-LDSC，**决定性**) → ⑤ run_cauchy_combination → ⑥ run_report(HTML，24GB 会 OOM)。
- **S-LDSC 核心方程**：`E[χ²_j] = 1 + N·Σ_c τ_c·ℓ(j,c)`，回归响应变量 = `Z²`（符号被丢弃）。**「燃料」= mean χ² − 1 = N·h²/M**。
- **`filter_sumstats_by_chisq`**：`chisq = Z**2`，默认 `chisq_max = max(0.001×N.max, 80) = 80`，删掉 chisq>80 的 SNP（即单 SNP 超额上限 79）。
- **ACAT cauchy 组合**：`cct_stat = Σ w_i·tan((0.5−p_i)π)`；`pval = 1 − cauchy.cdf(cct_stat)`（**单侧**）；per-spot p = `1−Φ(z)`（上尾）。
- **Z 的均匀缩放是 no-op**（S-LDSC t 统计量尺度不变），只有**选择性（非均匀）抬升**才改变显著性——而这在逻辑上是循环论证。
- **quick_mode 用预计算**：weights_hm3_no_hla、snp_gene_weight_matrix.h5ad、baseline、SNP_gene_pair、gencode.v46lift37 gtf（hg19）。
- **HapMap3 参考 ~1.2M SNP** 是 LDSC 参考；用户稀疏 GWAS SNP 大多不落在其中（这是稀疏位点被丢弃的主因）。
- 参考基因注释：marker score feather 里基因名在 `HUMAN_GENE_SYM` 列（16,331 基因 × 121,767 spots）；SNP-gene 权重矩阵 var 有 `Gene` 列（19,084 基因）。

---

## 7. 当前决策点（下一步选项）

已确认「公开数据 → gsMap」这条路走不通。摆在面前的选项（供新 agent / 用户选择）：

- **(a)** 把复发性流产 GCST011887 转成 gsMap sumstats 格式硬跑一次，作为**正式阴性对照**（预期仍为空，但能补齐证据链）。
- **(b)** 转向 **OLA1 的共定位 / 表达定位**路线（不走 S-LDSC，直接做 OLA1 在空转数据里的细胞类型表达定位 + 共定位分析）。
- **(c)** 换一个**有强多基因信号的代理性状**（如与 RPL 遗传相关的身高/BMI/激素水平等）跑 gsMap，看细胞类型定位结果作为间接证据。
- **(d)** 接受「RPL 无多基因信号」的结论，改从**罕见变异 / 家系 / 功能实验**角度研究 OLA1。

---

## 8. 关键术语速查

- **mean χ² / λ_GC**：多基因信号的金标准指标。mean χ² = mean(Z²)≈1 表示无信号；>1 表示多基因信号（S-LDSC 的燃料来源）。λ_GC = median(χ²)/0.4549。
- **fuel（燃料）**：mean χ² − 1 = N·h²/M，S-LDSC 能归因给注释的总超额。
- **S-LDSC**：stratified LD score regression，gsMap 第④步的核心。
- **ACAT / cauchy**：gsMap 第⑤步，把 per-spot 的 S-LDSC p 值组合成 per-cell-type 的显著性。
- **寡基因(oligogenic) vs 多基因(polygenic)**：前者=少数大效应位点，后者=海量小效应位点。gsMap 只吃后者。
