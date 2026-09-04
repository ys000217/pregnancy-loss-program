# 高变 CpG 筛选与聚类

从 ONT modkit `*.cpg.5mC.bed.gz` 构建位点矩阵，按方差筛选高变 CpG，对妊娠丢失队列做无监督聚类。

本地大矩阵、日志、临床表**不入库**；工作副本可放在仓库根目录的 `筛选高变CpG/`（已被 `.gitignore` 忽略数据文件）。

## 分析逻辑（当前）

### 1–2. 全队列 + NC：主要收获是「疑似 abnormal」

全队列高变 CpG 能分开 **clinical abnormal 富集枝**（约 46/47 abnormal），同枝上还有 **30 例 normal_case + 1 例 control**。去掉 clinical abnormal 后重算高变（job **09**），这 30 例仍成簇且与全队列枝完全重叠。

**叙事口径：** 这 30 例是**疑似 abnormal**（临床未标 abnormal，甲基化为 abnormal-like）。主要收获是用高变位点 + 敏感性分析，在 **7 / 11 / 12 周**找到这批疑似样本（**不在 8–10 周**）。名单：[`metadata/abnormal_like_normal_case_30.txt`](metadata/abnormal_like_normal_case_30.txt)。

全队列/NC 上 control 与其余 normal_case **整体仍分不开**；加大 k 无济于事。

### 3. 分周 8 / 9 / 10（job 11）：问 control vs normal_case

按周各算高变。两套：

| 集合 | 含义 | W8 | W9 | W10 |
|------|------|----|----|-----|
| **all** | 该周全部（含 abnormal） | 339 | 101 | 74 |
| **noabn** | 该周去掉 clinical abnormal | 308（72/236） | 93（29/64） | 66（10/56） |

因疑似 30 例不在 8–10 周，**noabn ≈ 干净主流 + 控孕周**。本步问题：

1. 周内 control 与 normal_case 能否分开？
2. 若有稳定亚群，跨周是否重复、是否随孕周变化？

本地分析脚本：`scripts/run_week8910_noabn_control_vs_normal.R`（Nature 出图规范）。

### 4. 8–10 周合并（job 12）

把 8+9+10 合成一层再算（`all` / `noabn`），与分周对照，看合并是否被某一周主导。

| 集合 | 预期 n |
|------|--------|
| **all** | **514** |
| **noabn** | **467** |

### 结果状态（入库文档口径）

| 步骤 | 状态 | 关键产物 |
|------|------|----------|
| 1 全队列高变聚类 | **完成** | 本地 `筛选高变CpG/`（大矩阵 gitignore）；名单与逻辑见上文 |
| 2 NC 敏感性 → 30 例 abnormal-like | **完成** | [`metadata/abnormal_like_normal_case_30.txt`](metadata/abnormal_like_normal_case_30.txt) |
| 3 8/9/10 周分周 all + noabn | **待跑/已跑看服务器** | `scripts/11_week8910_all_and_noabn_highvar.sh` |
| 4 8–10 周合并 all + noabn | **待跑** | `scripts/12_week8to10_pooled_all_and_noabn_highvar.sh` |

**生物学要点：** 高变空间主轴是 abnormal / abnormal-like vs 主流混合群，而不是 control vs normal_case 的干净分开；30 例 abnormal-like 与临床 abnormal 定义不符，但甲基化表型可复现。

## 服务器流水线（JSUB）

脚本在 `scripts/`。01–08 全队列准备；**09** 为 NC 敏感性测试；**10** 为主流样本孕周分层。

| 脚本 | 作用 |
|------|------|
| `01_make_filelist.sh` | 扫描 bed.gz → `sample_file_map.tsv` |
| `02_run_single_methylation.sh` | cov≥5 提取 beta → `CpG_long_tmp_serial/` |
| `03_count_cpg.sh` | 全队列位点覆盖计数 |
| `04_build_matrix_join.sh` | 按 `CpG_95pct.list` 拼全矩阵 |
| `05_compute_cpg_variance.sh` | 全队列方差（依赖 95% 白名单） |
| `06_select_high_variance_CpG.sh` | variance > 0.005 |
| `07_extract_RPL_study_matrix_648.sh` | meQTL 用 648 样本 bed 矩阵 |
| `08_build_top100k_matrix_for_pca.sh` | 全队列方差 top 100k 矩阵 |
| `09_NC_normal_control_highvar_cpg.sh` | 剔除 abnormal 后重算覆盖/方差（敏感性测试） |
| `10_GA_stratified_mainstream_highvar_cpg.sh` | 剔除 abnormal+like 后仅 8/9/10 周（与 11-noabn 样本集相同） |
| `11_week8910_all_and_noabn_highvar.sh` | **8/9/10 分周**：全样本 与 去 abnormal 各算组内高变 |
| `12_week8to10_pooled_all_and_noabn_highvar.sh` | **8–10 周合并**：全样本 与 去 abnormal 各算高变 |

09–12 注意：用 `bash script` 或 `jsub < script`（勿 `./script`，Windows CRLF 会 `/bin/bash^M`）；勿 `sort|head`（`pipefail` SIGPIPE）；临床表去掉 `\r`。

```bash
# 服务器：8–10 周合并 all + 去 abnormal
# jsub < analyses/highvar_cpg/scripts/12_week8to10_pooled_all_and_noabn_highvar.sh
```

产出：
- 分周（11）：`.../week8910_hvar/{all,noabn}/{W8,W9,W10}/`
- 合并（12）：`.../week8to10_pooled_hvar/{all,noabn}/W8_10/`

## 本地 R 分析

默认读取 `筛选高变CpG/` 下已下载的矩阵。

| 脚本 | 作用 |
|------|------|
| `run_clustering.R` | 步骤 1：全队列 top-N、Ward 树、Class3 对照 |
| `run_kscan_normal_subcluster.R` | 全队列 k 扫描 |
| `run_NC_matrix_cluster.R` | 步骤 2：`CpG_matrix_NC.tsv` 敏感性测试 |
| `check_cluster_overlap.R` | 全队列 abnormal 枝 vs NC 小簇重叠 |
| `check_nc_outputs.R` | 检查 09 产出、abnormal 泄漏 |
| `plot_gradient_scatter.R` | 全队列 / NC 的 PCA 与 PCoA 点图 |
| `summarize_nc_outlier_gw.R` | 30 例孕周分布 |
| `run_normal_control_hvar.R` | 旧路径（全队列 top100k 内重排方差，有位点池偏倚） |

```bash
Rscript analyses/highvar_cpg/scripts/run_NC_matrix_cluster.R
```
