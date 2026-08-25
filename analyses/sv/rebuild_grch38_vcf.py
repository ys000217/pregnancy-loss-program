import pandas as pd
import re

VCF  = r"E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf"
BED  = r"E:\genotype_data\liftover\GRCh38_breakpoints.unique.bed"
FAI  = r"E:\genotype_data\liftover\GRCh38_chr.fasta.fai"
OUT  = r"D:\ONT\S957.merged.GRCh38.vcf"
LOG  = r"D:\ONT\S957.merged.GRCh38.liftover.log"

# ---- 1. GRCh38 contig 长度 ----
contig_len = {}
with open(FAI) as f:
    for line in f:
        p = line.split("\t")
        contig_len[p[0]] = int(p[1])

# ---- 2. 已 liftover 断点: sv_idx -> side(L/R/S) -> (chrom, 1-based pos) ----
bed = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "svid", "svtype", "side"])
bed["sv_idx"] = bed["svid"].str.extract(r"SV(\d+)")[0].astype(int)
bed["s_letter"] = bed["svid"].str.extract(r":([SLR])$")[0]
bp = {}
for _, r in bed.iterrows():
    bp.setdefault(int(r.sv_idx), {})[r.s_letter] = (r.chrom, int(r.end))

# ---- 3. 读 VCF 头 ----
header_lines = []
chrom_line = None
with open(VCF) as f:
    for line in f:
        if line.startswith("#CHROM"):
            chrom_line = line
            break
        header_lines.append(line)

# ---- 4. 写新头(GRCh38) ----
out = open(OUT, "w", encoding="utf-8")
for line in header_lines:
    if line.startswith("##contig"):
        continue                     # 丢弃 CN1 contig
    if line.startswith("##reference"):
        continue                     # 丢弃旧 reference(若有)
    out.write(line)
for c in sorted(contig_len, key=lambda x: (x != "chrM", x != "chrY", x != "chrX",
                                           int(re.sub(r"\D", "", x) or 0))):
    out.write("##contig=<ID=%s,length=%d>\n" % (c, contig_len[c]))
out.write("##reference=GRCh38.p14\n")
out.write("##liftover=CN1_v0.8_to_GRCh38_CrossMap_bed;SOURCE=GRCh38_breakpoints.unique.bed\n")
out.write(chrom_line)

# ---- 5. 逐条重建数据 ----
log = open(LOG, "w", encoding="utf-8")
log.write("sv_idx\tsvtype\tstatus\tdetail\n")
n_full = n_partial = n_dropped = n_tra = 0
idx = 0
with open(VCF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        chrom, pos, rid = p[0], int(p[1]), p[2]
        info = p[7]
        svtype = end = chr2 = None
        for field in info.split(";"):
            if field.startswith("SVTYPE="):
                svtype = field[7:]
            elif field.startswith("END="):
                end = int(field[4:])
            elif field.startswith("CHR2="):
                chr2 = field[5:]
        b = bp.get(idx, {})
        if svtype in ("TRA", "BND"):
            # 易位: L 端 = CHROM/POS, R 端 = CHR2/END
            if "L" in b and "R" in b:
                Lc, Lp = b["L"]; Rc, Rp = b["R"]
                p[0] = Lc; p[1] = str(Lp)
                info = re.sub(r"END=\d+", "END=%d" % Rp, info)
                info = re.sub(r"CHR2=[^;]+", "CHR2=%s" % Rc, info)
                out.write("\t".join(p) + "\n")
                n_full += 1; n_tra += 1
                log.write("%d\t%s\tlifted\tTRA %s->%s\n" % (idx, svtype, chr2, Rc))
            else:
                n_dropped += 1
                log.write("%d\t%s\tdropped\tmissing breakpoint\n" % (idx, svtype))
        elif svtype == "INS":
            if "S" in b:
                Sc, Sp = b["S"]
                p[0] = Sc; p[1] = str(Sp)
                info = re.sub(r"END=\d+", "END=%d" % Sp, info)
                info = re.sub(r"CHR2=[^;]+", "CHR2=%s" % Sc, info)
                out.write("\t".join(p) + "\n")
                n_full += 1
                log.write("%d\t%s\tlifted\tINS\n" % (idx, svtype))
            else:
                n_dropped += 1
                log.write("%d\t%s\tdropped\tmissing breakpoint\n" % (idx, svtype))
        else:  # DEL / DUP / INV
            if "L" in b and "R" in b:
                Lc, Lp = b["L"]; Rc, Rp = b["R"]
                if Lc == Rc:
                    p[0] = Lc
                    p[1] = str(min(Lp, Rp))
                    info = re.sub(r"END=\d+", "END=%d" % max(Lp, Rp), info)
                    info = re.sub(r"CHR2=[^;]+", "CHR2=%s" % Lc, info)
                    out.write("\t".join(p) + "\n")
                    n_full += 1
                    log.write("%d\t%s\tlifted\n" % (idx, svtype))
                else:
                    n_partial += 1
                    log.write("%d\t%s\tpartial\tL/R lifted to different chrom (%s vs %s)\n"
                              % (idx, svtype, Lc, Rc))
            else:
                n_partial += 1
                log.write("%d\t%s\tpartial\tmissing one breakpoint\n" % (idx, svtype))
        idx += 1
out.close()
log.close()
print("总记录: %d" % idx)
print("fully lifted (写入 GRCh38 VCF): %d   (其中 TRA %d)" % (n_full, n_tra))
print("partial (仅一端 lifted 或两端异染色体): %d" % n_partial)
print("dropped (断点全缺): %d" % n_dropped)
print("输出:", OUT)
print("日志:", LOG)
