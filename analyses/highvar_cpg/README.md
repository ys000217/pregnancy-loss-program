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

### 3. 再剔除 abnormal-like 后，只在 8 / 9 / 10 周层内重算（下一步）

去掉 **47 例 abnormal + 30 例 abnormal-like** 后，按整周清点：5–7 周、11–12 周人数不足或几乎没有 control，**不进入本步**。只对人数够的三周，**各周单独**算覆盖与高变 CpG（不把不同周混在一层里）。

| 层 | Group3 | 预期 n（control / normal_case） |
|----|--------|--------------------------------|
| `W8` | 8 周 | 308（72 / 236） |
| `W9` | 9 周 | 93（29 / 64） |
| `W10` | 10 周 | 66（10 / 56）；control 偏少 |

合计约 **467**。服务器脚本：**10**（尚未跑）。

### 结果状态（入库文档口径）

| 步骤 | 状态 | 关键产物 |
|------|------|----------|
| 1 全队列高变聚类 | **完成** | 本地 `筛选高变CpG/`（大矩阵 gitignore）；名单与逻辑见上文 |
| 2 NC 敏感性 → 30 例 abnormal-like | **完成** | [`metadata/abnormal_like_normal_case_30.txt`](metadata/abnormal_like_normal_case_30.txt) |
| 3 8/9/10 周主流层内高变 | **待跑** | `scripts/10_GA_stratified_mainstream_highvar_cpg.sh` |

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
| `10_GA_stratified_mainstream_highvar_cpg.sh` | 再剔除 30 例 abnormal-like，仅 8/9/10 周层内重算 |

09/10 注意：`chmod +x` 或 `bash script`；勿 `sort|head`（`pipefail` 下 SIGPIPE / exit 141）；临床表去掉 Windows `\r`。

```bash
# 服务器（下一步）
# jsub < analyses/highvar_cpg/scripts/10_GA_stratified_mainstream_highvar_cpg.sh
```

产出目录（服务器）：`.../prepare_methylation/mainstream_by_GA/{W8,W9,W10}/`

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
