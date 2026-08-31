# GWAS × gsMap 分析模块

> **模块位置：** `analyses/gwas/`（妊娠丢失主仓库子模块）  
> 总项目说明见仓库根目录 [README.md](../../README.md)。

本模块整理中国 RPL 队列 GWAS 汇总统计的再分析，以及将 suggestive 位点纳入 **gsMap**（空间转录组 × GWAS）尝试细胞类型定位的流程与结果摘要。

## 一、分析概要

| 项目 | 内容 |
|---|---|
| 表型 | 复发性流产（RPL），中国队列，N≈608 |
| 原始输入（本地，不入库） | `combine_gwas_result_v1.Status.glm.logistic.hybrid`（hg38） |
| Suggestive 阈值 | **P &lt; 1×10⁻⁴**（小样本宽松线；全基因组线仍为 5×10⁻⁸） |
| 基因组膨胀 | λ_GC ≈ **1.015** |
| Suggestive 位点数 | **87**（0 个过 5×10⁻⁸） |
| 主信号区 | chr2:174.07–174.22 Mb（hg38，OLA1 附近，49/87） |

## 二、脚本

| 脚本 | 功能 |
|---|---|
| `scripts/reanalyze_gwas_suggestive.py` | 读 PLINK hybrid：λ_GC、曼哈顿/QQ、按 suggestive 线导出位点 |
| `scripts/build_suggestive87_sumstats.py` | 将 87 个命中 liftOver→rsID，写入 gsMap 格式 sumstats（\|Z\|=8 boost） |
| `scripts/run_gsmap_suggestive87.sh` | WSL：小鼠胚胎 quick_mode（示例 ST）× suggestive-87 |
| `scripts/run_gsmap_sug87_CS17.sh` | WSL：人胚胎 CS17 HESTA Stereo-seq × suggestive-87 |
| `scripts/summarize_cs17_sug87.py` | 汇总 CS17 Cauchy / spot 指标并拷贝小结果 |
| `scripts/gsmap_healthcheck.py` | mean χ² / λ_GC / 与 1000G EUR 重叠体检 |
| `scripts/meta_rpl_ea.py` 等 | 东亚 meta 相关脚本（可选；表型差异大时不作敏感性主依据） |

### 运行（示例）

```bash
# Windows：GWAS 再分析（需 pandas/numpy/scipy/matplotlib）
python scripts/reanalyze_gwas_suggestive.py

# WSL：构建 sumstats + gsMap（需 gsmap 环境与本地资源路径）
source /home/administrator/gsmap_env/bin/activate
python scripts/build_suggestive87_sumstats.py
bash scripts/run_gsmap_suggestive87.sh      # 小鼠胚胎示例
bash scripts/run_gsmap_sug87_CS17.sh       # 人胚胎 CS17 HESTA
```

本地大文件与 gsMap 资源路径见 `docs/DATA_PATHS.md`（不提交原始 hybrid / 大矩阵）。

## 三、主要结果文件（本目录 `results/`）

| 文件 | 说明 |
|---|---|
| `GWAS_suggestive_hits_1e-4.csv` | 87 个 suggestive 位点 |
| `GWAS_reanalysis_summary_1e-4.txt` | λ_GC 与计数摘要 |
| `combine_GWAS_Manhattan_suggestive_1e-4.png` / `*_QQ_*.png` | 曼哈顿图 / QQ 图 |
| `RPL_suggestive87_mapping_report.txt` | 87 位点 → hg19/rsID/sumstats 映射 |
| `RPL_cauchy_celltype_level.csv.gz` | 既往全基因组 RPL gsMap Cauchy |
| `RPL_OLA1boost_cauchy_celltype_level.csv.gz` | 既往 OLA1 抬升实验 Cauchy |
| `RPL_sug87_cauchy_celltype_level.csv.gz` | suggestive-87 × 小鼠胚胎（best p_cauchy=0.116，阴性） |
| `RPL_sug87_gsmap_summary.md` | 小鼠胚胎 gsMap 定位摘要 |
| `CS17_RPL_sug87_cauchy.csv.gz` | suggestive-87 × **人胚胎 CS17** Cauchy（best p_cauchy=0.146，Eye，阴性） |
| `CS17_RPL_sug87_summary.txt` | CS17 spot / Cauchy 摘要 |
| `GWAS_suggestive_hits_1e-3.csv` / `GWAS_reanalysis_summary_1e-3.txt` | 更松阈值 P&lt;1×10⁻³ 的探索性再分析（非正式主结论） |
| `CS17_RPL_sug1e3_*.csv.gz` / `*_summary.txt` | suggestive-1e-3 × CS17（探索；主叙述仍以 87 / 1e-4 为准） |

## 四、gsMap 解读要点

gsMap / S-LDSC 依赖**全基因组多基因信号**（mean χ² ≫ 1）。本队列 mean χ²≈1.0–1.08；即便将可映射的 suggestive SNP 抬到 \|Z\|=8，在小鼠胚胎示例与人胚胎 CS17 Stereo-seq 上均为阴性。换相关人胚胎片子仍不显著，进一步支持瓶颈在 GWAS 侧燃料不足。

**正式主结论**：suggestive **P&lt;1×10⁻⁴ → 87 位点**；0 个过 5×10⁻⁸；λ_GC≈1.015；主峰 chr2 OLA1 附近。

## 五、不入库内容

- `*.Status.glm.logistic.hybrid`、公开大型 sumstats（Laisk / FinnGen 等）
- gsMap 资源目录、`.h5ad`、完整 LD score、spot 级 spatial_ldsc 大表
- Python `.venv/`
