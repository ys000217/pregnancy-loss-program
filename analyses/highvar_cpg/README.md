# 高变 CpG 筛选与聚类

从 ONT modkit `*.cpg.5mC.bed.gz` 构建位点矩阵，按方差筛选高变 CpG，对妊娠丢失队列做无监督聚类。

本地大矩阵、日志、临床表**不入库**；工作副本可放在仓库根目录的 `筛选高变CpG/`（已被 `.gitignore` 忽略数据文件）。

## 生物学结论（当前）

全队列高变 CpG 主要分开 **abnormal-like 枝 vs 主流混合群**，而不是 clinical 三分（abnormal / normal_case / control）。

去掉 abnormal 后，在 normal_case + control 上重算方差，仍分不开两组整体；但会稳定抠出 **30 例 normal_case + 1 例 control**。这 31 人与全队列 k=3 **簇 3**（abnormal 富集枝）完全重叠。

嵌套标签建议：`abnormal` / `normal_outlier`（30 例）/ `mainstream`。

## 服务器流水线（JSUB）

脚本在 `scripts/`。01–08 为全队列准备；**09** 为 NC 专用（剔除 abnormal 后重算覆盖与方差）。

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
| `09_NC_normal_control_highvar_cpg.sh` | **NC 专用**覆盖 + 方差 + top 10000 矩阵 |

09 注意：`chmod +x` 后再 `jsub`；`sort\|head` 在 `pipefail` 下会 SIGPIPE（exit 141），脚本已改为完整 sort 再截断。临床表需去掉 Windows `\r`，否则 abnormal 剔除失败。

## 本地 R 分析

默认读取 `筛选高变CpG/` 下已下载的矩阵。

| 脚本 | 作用 |
|------|------|
| `run_clustering.R` | 全队列 top-N 梯度、FM、k=3 vs Class3、热图 |
| `run_kscan_normal_subcluster.R` | 全局 k 扫描 + normal_case 子集重算方差 |
| `run_normal_control_hvar.R` | 旧路径：在全队列 top100k 内对 NC 重排方差（有位点池偏倚） |
| `run_NC_matrix_cluster.R` | **主分析**：`CpG_matrix_NC.tsv` 上 Ward 聚类 |
| `check_nc_outputs.R` | 检查 09 产出行列、abnormal 泄漏 |
| `check_cluster_overlap.R` | 全队列簇 3 vs NC 簇 2 成员重叠 |

```bash
Rscript analyses/highvar_cpg/scripts/run_NC_matrix_cluster.R
```
