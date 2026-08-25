import time

VCF = r"E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf"
FID = r"D:\ONT\matrix_fid.txt"
OUT = r"D:\ONT\sv_carriers.tsv"

# 648 methylation FIDs (in matrix column order)
fids = [l.strip() for l in open(FID, encoding="utf-8") if l.strip()]
fid_set = set(fids)

# read header to locate sample columns
with open(VCF) as f:
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            break
samples = header[9:]
col_of_fid = {}
for i, s in enumerate(samples):
    if s in fid_set:
        col_of_fid[s] = i
print("fids=%d found_in_vcf=%d" % (len(fids), len(col_of_fid)), flush=True)

# target column indices in VCF data (aligned to fids order)
target_cols = [col_of_fid[f] for f in fids]

out = open(OUT, "w", encoding="utf-8")
out.write("sv_idx\tsvtype\tn_carriers\tcarriers\n")
idx = 0
t0 = time.time()
with open(VCF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split("\t")
        info = p[7]
        svtype = ""
        for field in info.split(";"):
            if field.startswith("SVTYPE="):
                svtype = field[7:]
                break
        if svtype == "TRA":
            idx += 1
            continue
        carriers = []
        for k, c in enumerate(target_cols):
            gt = p[9 + c]
            colon = gt.find(":")
            if colon > 0:
                gt = gt[:colon]
            if "1" in gt:
                carriers.append(fids[k])
        out.write("%d\t%s\t%d\t%s\n" % (idx, svtype, len(carriers), ",".join(carriers)))
        idx += 1
        if idx % 20000 == 0:
            print("processed %d records, %.1f s" % (idx, time.time() - t0), flush=True)
out.close()
print("DONE total_records=%d total_sec=%.1f" % (idx, time.time() - t0), flush=True)
