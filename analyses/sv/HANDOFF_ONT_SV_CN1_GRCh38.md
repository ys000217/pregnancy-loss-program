# ONT SV 项目交接文档（供下一个 agent 无缝接手）

> 本文件记录从项目启动至今的**全部上下文、已做工作、关键结论、当前状态、下一步任务**。
> 下一个 agent 读完本文档即可直接继续，无需重新探索。

---

## 1. 项目目标

复现一篇论文的方法：**从结构变异(SV)的 breakpoint 富集，寻找 DNA 甲基化变化的区域**。

因此核心需求是：**拿到 SV breakpoint 的坐标信息**，并且这些坐标要在 **GRCh38 参考**上，才能和甲基化数据（在 GRCh38 上）做交集/富集分析。

---

## 2. 数据文件清单（均在 `E:\genotype_data\`）

| 文件 | 大小 | 说明 |
|---|---|---|
| `S957.1kGPont.merged.vcf` | 3.6 GB | 原始合并 VCF，1984 样本（957 临床 + 外部 1kGP），SURVIVOR 合并，95,635 条记录 |
| `S957.Sniffles2.merged.vcf` | 9.8 GB | 957 样本，Sniffles2 合并，551,637 条，含 **202,886 条 BND**（全为染色体内） |
| `S957.hc.merged.addmiss.vcf` | 2.95 GB | 957 样本，119,798 条，**0 BND、31 TRA** |
| `S957.1kGPont.merged.649clin_1027kgp.vcf` | 3.16 GB | **★ 筛选后的最终文件**：1676 样本 = 649 临床 + 1027 个 1kGP，95,635 条 |
| `临床信息表0.2.xlsx` | — | 临床表，649 样本，列：Sample_ID/Group1/Group2/Gender/Age/Gravida/Loss/Para/Gestational_Week/Group3/CRL/Complications/Group4 |
| `CN1_v0.8.genome.fasta.gz` | 876,039,413 字节 | **★ 参考基因组 fasta（已下载并核对）**，见 §5 |
| `S957.Sniffles2.merged.anno.txt` | — | ANNOVAR 注释 |
| `S957.hc.merged.addmiss.anno.txt` | — | ANNOVAR 注释 |

工作目录：`d:\ONT\`（含筛选脚本 `filter_vcf.awk`、`clean_info.awk`、`exclude_308.txt`）。

---

## 3. 已完成的工作（按时间顺序）

1. **样本筛选**：从 `S957.1kGPont.merged.vcf`（1984 样本）中，只保留 **649 个临床表样本 + 1027 个 1kGP 样本**，剔除 308 个非临床 S957 样本，生成 `S957.1kGPont.merged.649clin_1027kgp.vcf`。只做筛选，保持 VCF 格式合法。
2. **INFO 清理**：删掉会误导的 `AC/AN/SUPP/SUPP_VEC`（这些反映原始 1984 样本的统计，筛选后不准确）。保留的 INFO 字段：`CIEND, CIPOS, CHR2, END, MAPQ, RE, IMPRECISE, PRECISE, SVLEN, SVMETHOD, SVTYPE, STRANDS`。FORMAT = `GT:IS:OT:DV:DR`。
3. **变异类型统计**（用于复现论文）：三份 VCF 的 SVTYPE 已统计。BND 结论：`S957.Sniffles2.merged.vcf` 有 202,886 条 BND；SURVIVOR 合并文件里 BND 被重命名为 **TRA**（`S957.hc.merged.addmiss.vcf` 有 31 条 TRA）。
4. **参考基因组识别**（本阶段核心，见 §4、§5）：确认三份 VCF 用的参考是 **CN1_v0.8**，不是 GRCh38。

---

## 4. 关键结论（最重要的部分）

### 4.1 参考基因组 = CN1_v0.8（T2T-CN1）

- 三份 SV VCF 全部基于 **CN1_v0.8** —— 汉族男性完整单倍体 T2T 参考基因组。
- 出处：*Cell Research* 2023, Yang C. et al., "The complete and fully-phased diploid genome of a male Han Chinese"。DOI: `10.1038/s41422-023-00849-5`。
- **Novogene 交付时把它叫 "T2T.v0.8"**（= "T2T" 前缀 + 官方版本号 "v0.8"），这就是用户当初看到的名字。
- **不是 GRCh38**。用户此前误以为是 GRCh38，已纠正。

### 4.2 排除过的其它参考（均不匹配）

GRCh38、T2T-CHM13v2.0(hs1)、T2T-YAO、Han1(HG00621)——都已通过染色体长度比对排除。

### 4.3 最终铁证：25/25 染色体长度逐碱基一致

已下载 CN1_v0.8 fasta 并解压统计 25 条序列长度，与 VCF 头部 `##contig` **完全一致**（见 §5 表）。

---

## 5. 参考基因组 CN1_v0.8 详情

### 5.1 官方信息（CNCB/GWH API 返回）

```
assemblyName    : "CN1_v0.8"
assemblyLevel   : "Complete"
accession       : GWHCBHP00000000
bioproject      : PRJCA016397
biosample       : SAMC1215858
submitter       : Medicine School, Zhejiang University（浙江大学医学院）
publication     : Cell Research 2023, 33(10):745-61
```

- 元数据 API：`https://ngdc.cncb.ac.cn/gwh/api/public/assembly/GWHCBHP00000000`
- 下载直链（**已验证 HTTP 200，无需登录**）：
  `https://download.cncb.ac.cn/gwh/Animals/Homo_sapiens_CN1_v0.8_GWHCBHP00000000/GWHCBHP00000000.genome.fasta.gz`
- 其它存放地：CNGB accession `CNA0069006`（BioProject `CNP0004252`）；CRGD 浏览器 `genome.zju.edu.cn`（只放 v0.6，另含 `CHM13_to_CN1.paf`、`CHM13_vs_CN1.SV.bb` 比对文件）；GitHub `T2T-CN1/CN1`（只有代码+BED，无 fasta）。

### 5.2 25 条染色体长度（= VCF 的 contig 长度，逐碱基一致）

fasta 内序列名为 `GWHCBHP00000001`…`GWHCBHP00000025`，顺序 = chr1…chr22、chrX、chrY、chrM。Novogene 把名字换成了 chr 形式，但顺序和长度完全一致。

| # | 染色体 | 长度 | # | 染色体 | 长度 |
|---|---|---|---|---|---|
| 1 | chr1 | 254,717,395 | 14 | chr14 | 98,522,907 |
| 2 | chr2 | 242,636,919 | 15 | chr15 | 94,829,848 |
| 3 | chr3 | 199,957,529 | 16 | chr16 | 90,698,698 |
| 4 | chr4 | 191,432,471 | 17 | chr17 | 83,813,795 |
| 5 | chr5 | 183,988,138 | 18 | chr18 | 80,489,117 |
| 6 | chr6 | 171,967,578 | 19 | chr19 | 62,832,331 |
| 7 | chr7 | 158,872,745 | 20 | chr20 | 66,744,663 |
| 8 | chr8 | 145,796,002 | 21 | chr21 | 45,247,884 |
| 9 | chr9 | 137,877,317 | 22 | chr22 | 50,276,510 |
| 10 | chr10 | 134,863,434 | X | chrX | 157,365,100 |
| 11 | chr11 | 135,971,833 | Y | chrY | 65,593,761 |
| 12 | chr12 | 133,770,623 | M | chrM | 16,571 |
| 13 | chr13 | 105,970,862 | | **合计** | **3,094,254,031** |

---

## 6. 已生成/已下载文件

| 文件 | 位置 | 说明 |
|---|---|---|
| 筛选后 VCF | `E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf` | ★ 待 liftover 的目标文件 |
| CN1 参考 fasta | `E:\genotype_data\CN1_v0.8.genome.fasta.gz` | ★ 已下载、gzip 校验通过、25 条长度已核对 |
| 筛选脚本 | `d:\ONT\filter_vcf.awk` | 按排除列表删样本列 |
| INFO 清理脚本 | `d:\ONT\clean_info.awk` | 删 AC/AN/SUPP/SUPP_VEC |
| 排除样本列表 | `d:\ONT\exclude_308.txt` | 308 个要剔除的样本 ID |

---

## 7. 下一步任务：把 SV 断点 liftover 到 GRCh38

这是**尚未开始**的核心任务。目标：把 `S957.1kGPont.merged.649clin_1027kgp.vcf` 的坐标从 CN1 转到 GRCh38。

### 7.1 待办步骤

1. 下载 **GRCh38 (hg38) fasta**（UCSC）：`https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`（GRCh38 chrom.sizes 已在 `d:\ONT` 下，或 `.../hg38.chrom.sizes`）。
2. 生成 CN1 fasta 的 `.fai`（已下 fasta，但 .fai 未下；用 `samtools faidx` 或脚本生成）。
3. **建 CN1→GRCh38 chain**（WSL 里）：
   - 快速近似：`minimap2 -cx asm5 --cs -t 8 GRCh38.fa CN1_v0.8.fasta | ...` 转 chain；
   - 或精确：`nucmer`/SyRI（慢但复杂区域更准）；
   - 备选：CRGD 有 `CHM13_to_CN1.paf`，可做 CN1→CHM13 再加 CHM13→GRCh38 的公开 chain。
4. **liftover VCF**：`picard LiftoverVcf` 或 `CrossMap vcf`，输入 CN1 chain + 目标 VCF + hg38 fasta（用于更新 REF），输出 GRCh38 版 VCF。注意要**重写 header 的 `##contig`/`##reference`**。
5. 从 GRCh38 VCF 生成 **breakpoint BED**（DEL/DUP/INV = POS+END 两个断点；INS = POS 单点；TRA/BND = (CHROM,POS)+(CHR2,END) 两端），并统计 **掉点率**（CN1 特有区域会掉点）。

### 7.2 关键注意事项（务必遵守）

- **liftover ≠ 重新 call**。它只是把坐标"翻译"过去，不重新验证变异在 GRCh38 下是否成立；REF/ALT 是 CN1 的（除非用目标 fasta 更新 REF）。
- **有损**：CN1 独有、GRCh38 没有的区域（着丝粒/rDNA/大片段）对应断点会掉点；跨复杂区域的 SV 坐标可能偏移。
- 对"断点富集找甲基化区域"这个目的，**liftover 坐标够用且是标准做法**，但方法里要交代 chain 来源、工具、掉点率。
- **SV 断点定义**：DEL/DUP/INV 用 POS+END；INS 用 POS（END==POS）；TRA（=BND）用 (CHROM,POS) 和 (CHR2,END) 两个断点。
- SV 的 `END`、`CHR2` 都要一起 lift，`SVLEN` 要更新。
- 版本注意：CRGD 只放 CN1 **v0.6**（`CN1_combine.v0.6`，与 v0.8 有微小差异：24/25 条 <2.5kb，chr21 差 1.4Mb rDNA），**不要用 v0.6 做 liftover**，要用已下的 **v0.8**。

---

## 8. 环境与工具

- OS：Windows 10 Pro，Git Bash（POSIX sh）。工作目录 `d:\ONT`，数据在 `E:\genotype_data`。
- Git Bash 里**没有** python/bcftools/conda/samtools；有 `gzip`、`curl`、`awk`、`grep`、`unzip`、`pdftotext`。
- **WSL Ubuntu-22.04 可用**（重工具装这里：minimap2、samtools、picard/CrossMap、nucmer）。
- 记忆目录：`C:\Users\Administrator\.claude\projects\d--ONT\memory\`（含 `ont-reference-genome-is-cn1.md`，已记录参考基因组结论）。

---

## 9. 可复现验证命令（重要）

**从 fasta 提取 25 条长度（已用，结果见 §5.2）：**
```bash
cd /e/genotype_data
gzip -dc CN1_v0.8.genome.fasta.gz | awk '
/^>/ { if (name != "") print name, len; name=substr($1,2); len=0; next }
{ len += length($0) }
END { if (name != "") print name, len }'
```

**从 VCF 头部提取 contig 长度：**
```bash
head -60 /e/genotype_data/S957.1kGPont.merged.649clin_1027kgp.vcf | grep '^##contig'
```

**参考基因组元数据：**
```bash
curl -s "https://ngdc.cncb.ac.cn/gwh/api/public/assembly/GWHCBHP00000000"
```

---

## 10. 一句话总结（给下一个 agent）

> 项目要复现"SV 断点富集找甲基化区域"论文。已有筛选好的 VCF（649 临床 + 1027 个 1kGP，`S957.1kGPont.merged.649clin_1027kgp.vcf`），参考基因组已确认是 **CN1_v0.8（汉族 T2T 单倍体参考，=Novogene 的 "T2T.v0.8"）**，fasta 已下载到 `E:\genotype_data\CN1_v0.8.genome.fasta.gz` 并逐碱基核对无误。**当前唯一未完成的核心任务：把这份 VCF 的 SV 断点 liftover 到 GRCh38（建 chain → LiftoverVcf → 断点 BED + 掉点统计）**。环境见 §8，注意事项见 §7.2。

---

## 11. Liftover 结果（2026-08-17 已完成）

### 11.1 数据流向澄清

3 份原始 VCF（全部基于 CN1_v0.8，25 contig）：

| 原始 VCF | 格式 | 样本数 | 说明 |
|---|---|---|---|
| S957.1kGPont.merged.vcf | VCFv4.1 | 1984 | SURVIVOR 合并，95,635 条 |
| S957.Sniffles2.merged.vcf | VCFv4.2 | 957 | Sniffles2 单工具，551,637 条 |
| S957.hc.merged.addmiss.vcf | VCFv4.2 | 957 | 高置信子集，119,798 条 |

**当前分析只用第 1 份** S957.1kGPont.merged.vcf，经两步处理：
1. filter_vcf.awk + exclude_308.txt：1984 -> 1676 样本（649 临床 + 1027 个 1kGP）。
2. clean_info.awk：删 INFO 的 SUPP/SUPP_VEC/AC/AN。

得到目标文件 S957.1kGPont.merged.649clin_1027kgp.vcf（1676 样本，95,635 条）。另外两份 VCF 未进入 liftover。

### 11.2 已完成步骤

1. WSL 装工具：Miniconda + env "liftover"（minimap2 2.31 / samtools 1.24 / bcftools 1.24 / CrossMap 0.7.3）。注意 WSL 无 sudo、无 pip、github 被墙，但 anaconda/pypi 可达，故走 conda。
2. 两侧 fasta 重命名（关键坑：两侧命名不一致，必须先统一）：
   - CN1：GWHCBHP00000001..25 -> chr1..chrM（按头部 OriSeqID）。-> liftover/CN1_chr.fasta
   - GRCh38：已有 d:/ONT/GRCh38.p14.fasta 是 GenBank NC_000001.11 命名（705 条），重命名并只留 25 条染色体 -> liftover/GRCh38_chr.fasta
3. 生成 .fai：samtools faidx 两个 fasta。
4. 建 chain：minimap2 -cx asm5 -c --secondary=no -t16 CN1_chr.fasta GRCh38_chr.fasta（579s，峰值 24.5 GB，9970 条对齐）-> 自写 paf2chain.py 转 chain（9970 条）。chain 覆盖 CN1 的 94.86%。
5. 提取断点：extract_breakpoints.awk，DEL/DUP/INV=POS+END 两断点，INS=POS 单点，TRA=(CHROM,POS)+(CHR2,END) 两端。-> CN1_breakpoints.bed（145,765 断点）。
6. liftover 断点：CrossMap bed（**不是 CrossMap vcf**，原因见 11.4）-> GRCh38_breakpoints.bed。

### 11.3 掉点率

断点级：145,765 -> lifted 96,506（66.21%），dropped 49,259（33.79%）。
SV 级（95,635 个 SV）：fully lifted 62,442（65.29%），partially 7,249（7.58%），fully dropped 25,944（27.13%）。

按 SVTYPE 断点掉点率：

| SVTYPE | 断点数 | 掉点数 | 掉点率 |
|---|---|---|---|
| DEL | 99,750 | 39,327 | 39.43% |
| INS | 45,505 | 9,878 | 21.71% |
| DUP | 374 | 39 | 10.43% |
| INV | 74 | 7 | 9.46% |
| TRA | 62 | 8 | 12.90% |

掉点主要来自 CN1 特有/异染色质区（chr9 着丝粒旁、acrocentric p 臂 rDNA 等），是 T2T->GRCh38 的生物学预期损失，不是技术错误。

### 11.4 关键教训（务必读）

- CrossMap 0.7.3 的 vcf 模式对 SV 是错的：它只 lift POS（第一个碱基），把 INFO 的 END 直接写成 lifted 后的 POS+1，且完全不更新 CHR2。因此对 DEL/DUP/INV/TRA 直接 CrossMap vcf 会得到错误的 END 和 CHR2。正确做法是先把断点抽成 BED，再 CrossMap bed 逐个断点 lift（本项目已这样处理）。
- chain 方向约定：CrossMap/UCSC 里 chain 的 tName = 源基因组（输入），qName = 目标基因组（输出）。要 CN1->GRCh38，用 minimap2 CN1.fa GRCh38.fa（tname=CN1, qname=GRCh38）。
- paftools.js 的 chain 命令在 minimap2 2.17/2.26/2.31 里都已被移除，故自写 paf2chain.py（用 -c 的 cg cigar tag 转 chain，已用合成数据 + CrossMap 验证正确）。
- 1,437 个断点因重复序列 multi-mapping 到多处 GRCh38 位点（在 GRCh38_breakpoints.bed 里有多行）。

### 11.5 产物文件

工作目录 E:/genotype_data/liftover：

| 文件 | 说明 |
|---|---|
| CN1_chr.fasta (+.fai) | CN1 重命名 fasta（25 条，chr 命名） |
| GRCh38_chr.fasta (+.fai) | GRCh38 重命名 fasta（25 条，chr 命名） |
| CN1_to_GRCh38.paf | minimap2 原始对齐（7 MB） |
| CN1_to_GRCh38.chain | UCSC chain（9970 条） |
| CN1_breakpoints.bed | CN1 断点（145,765 条，6 列：chrom/start/end/id/svtype/side） |
| GRCh38_breakpoints.bed | 星号 GRCh38 断点（lifted，100,003 行 = 96,506 唯一断点） |
| GRCh38_breakpoints.bed.unmap | 掉点断点（49,259 条） |

脚本 d:/ONT：extract_breakpoints.awk、paf2chain.py、analyze_liftover.py、run_alignment.sh、run_liftover.sh、rename_cn1.awk、rename_grch38.awk、setup_conda.sh。

### 11.6 下一步（供后续 agent）

1. 用 GRCh38_breakpoints.bed 与甲基化数据（GRCh38）做交集/富集（可先对每个断点取 ±窗口，如 ±100bp / ±1kb）。
2. 如需"完整 lifted VCF"（含正确的 END/CHR2/REF），需在断点 lift 结果基础上回填重建（当前只有断点 BED，未重建 VCF）。
3. 1,437 个 multi-mapping 断点需按分析需要去重或标记 ambiguous。
