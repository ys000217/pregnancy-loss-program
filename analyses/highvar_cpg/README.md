# 高变 CpG 筛选与聚类

从 ONT modkit `*.cpg.5mC.bed.gz` 构建位点矩阵，按方差筛选高变 CpG，对妊娠丢失队列做无监督聚类。

本地大矩阵、日志、临床表**不入库**；工作副本可放在仓库根目录的 `筛选高变CpG/`（已被 `.gitignore` 忽略数据文件）。

## 分析逻辑（当前）

三步嵌套，后一步都建立在前一步已经确认的结构上。

### 1. 全队列高变 CpG：abnormal 可分，并带出 30 例 normal_case

在**全体样本**上按方差取高变 CpG（Ward + Pearson）。主结构是两枝，不是 clinical 三分：

- **abnormal 富集枝**：约 46/47 临床 abnormal，外加 **30 例 normal_case + 1 例 control**
- **主流混合群**：其余 control 与 normal_case 混在一起

k 加大只会把主流群切碎，不会把 control / normal 整体切开。全队列脚本 `01–08` + 本地 `run_clustering.R` / `run_kscan_normal_subcluster.R` 已完成，不必重跑。

### 2. 敏感性测试：去掉 abnormal 后重算高变位点，这 30 例仍在

剔除临床 abnormal 后，只在 control + normal_case 上重算覆盖与方差（job **09**）。这 30 例 normal_case 仍单独成簇，且与全队列 abnormal 枝上的非 abnormal 成员 **100% 重叠**。

因此：这 30 例**不符合临床 abnormal 定义**，但在高变甲基化空间上是 **abnormal-like**。同簇 1 例 control（`zz230300`，11 周）记为同枝成员，不改临床标签。

这 30 例在 11–12 周富集，但 **82 例 11–12 周 normal_case 里只有 28 例进该簇**，晚孕周不是充分条件。名单：[`metadata/abnormal_like_normal_case_30.txt`](metadata/abnormal_like_normal_case_30.txt)。

job 09 与本地 `run_NC_matrix_cluster.R` / `check_cluster_overlap.R` 已完成，不必为改叙述而重跑。

### 3. 8 / 9 / 10 周组内高变（此前未做）

按整周、在人数足够的 8/9/10 周**各算一套**高变 CpG。两套样本集合（脚本 **11**）：

| 集合 | 含义 | W8 | W9 | W10 |
|------|------|----|----|-----|
| **all** | 该周全部样本（含 abnormal） | 339（abn 31 / ctrl 72 / nor 236） | 101（8 / 29 / 64） | 74（8 / 10 / 56） |
| **noabn** | 该周去掉 clinical abnormal | 308（72 / 236） | 93（29 / 64） | 66（10 / 56） |

脚本 **10**（再去掉 30 例 like）在这三周上与 **11-noabn 样本集相同**，因为 like 不在 8/9/10 周。可只跑 11。

30 例 abnormal-like 孕周：7 周 2 例、11 周 9 例、12 周 19 例（SPL 10 / RPL 20）；**8/9/10 周为 0**。

### 4. 8–10 周合并高变（此前未做；脚本 12）

把 8+9+10 周**合成一层**再算覆盖与方差（相对脚本 11 只改「是否分周」）：

| 集合 | 含义 | 预期 n |
|------|------|--------|
| **all** | 8–10 周全部样本（含 abnormal） | **514**（339+101+74） |
| **noabn** | 8–10 周去掉 clinical abnormal | **467**（308+93+66） |

产出一层目录：`W8_10/`。like 仍不在这三周，noabn 不必再剔 like。

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
