# 集群环境安装说明

## 为什么刚才 `conda env create` 失败

1. **Manta 不能和 Python 3 装进同一个 conda 环境**  
   Bioconda 里的 `manta` 仍依赖 Python 2.7，与 `python>=3.10` 冲突，solver 会直接失败。

2. **部分新版本包要求 glibc ≥ 2.28**  
   若集群是较老的 CentOS 7（glibc 2.17），solver 可能继续报 `__glibc >=2.28` 相关错误。

本仓库已改为：**conda 只装 Python 3 工具；WGS 断点 SV 默认用 DELLY**（与 Manta 同属 DEL/DUP 筛选用途）。若必须用 Manta，见下文「可选：单独装 Manta」。

---

## 标准安装（推荐）

```bash
cd /data/jhinno/appform/data/share/5250028/PrivateShareGroup/5250028_songyang/cnv

# 若之前创建失败，先删半成品环境
conda env remove -n cnv10x -y 2>/dev/null || true

conda env create -f env/environment.yml
conda activate cnv10x
bash scripts/check_env.sh
```

`check_env.sh` 要求：**DELLY 或 Manta 至少有一个**；其余 fastp / bwa-mem2 / samtools / mosdepth / sniffles / cnvpytor / bcftools / bedtools 必须存在。

### 脚本报 `$'\r': command not found`（Windows 换行）

从 Windows 拷来的 `.sh` 常是 CRLF。在 cnv 目录执行：

```bash
bash scripts/fix_crlf.sh
# 或：find . \( -name '*.sh' -o -name 'config.sh' \) -exec sed -i 's/\r$//' {} +
bash scripts/check_env.sh
```

---

## 若仍报 glibc / python 冲突

先看系统 glibc：

```bash
ldd --version | head -1
```

**方案 A：用 mamba 重试（有时比 conda 好解）**

```bash
conda install -n base -c conda-forge mamba
mamba env create -f env/environment.yml
```

**方案 B：集群 module（若有）**

```bash
module avail samtools bcftools bwa
# 加载站点已有模块后，只 conda 装 sniffles cnvpytor delly 等缺的
```

**方案 C：固定更老的 Python 3.10 build（glibc 2.17 集群）**

```bash
conda create -n cnv10x -c conda-forge -c bioconda \
  python=3.10.13 \
  fastp bwa-mem2 samtools sambamba mosdepth sniffles delly cnvpytor bcftools bedtools htslib
```

---

## 可选：单独装 Manta（不放进 cnv10x）

Manta 自带 Python 2 运行时，应**独立目录**安装，不要 `conda install manta` 进 cnv10x：

```bash
# 示例：解压官方 Linux 包（版本与路径按你下载的文件改）
mkdir -p ~/tools/manta-1.6.0
tar -xzf manta-1.6.0.centos6_x86_64.tar.bz2 -C ~/tools/manta-1.6.0 --strip-components=1

export PATH=~/tools/manta-1.6.0/bin:$PATH
configManta.py --help   # 能跑即可

# 写入 ~/.bashrc 或作业脚本：
# export PATH=~/tools/manta-1.6.0/bin:$PATH
```

装好后 `05_wgs_sv.sh` 会**优先用 Manta**；没有 Manta 时自动用 DELLY。

---

## 可选工具

| 工具 | 作用 | 缺了怎么办 |
|------|------|------------|
| spectre | ONT 大片段 depth CNV | 只用 CNVpytor |
| AnnotSV | HIGH CNV 注释 | 只输出 BED，不注释 |
| Manta | WGS 断点 SV | 用 DELLY（默认） |

---

## 装完还要准备的参考数据

- `REF_FASTA`：与 ONT BAM `@SQ` 一致的 GRCh38
- `TR_BED`：Sniffles2 tandem-repeat BED（如 `human_GRCh38_TR.bed`）

在 `config.sh` 里改路径，不要拷原始 FASTQ/BAM。
