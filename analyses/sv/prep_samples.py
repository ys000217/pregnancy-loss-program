import pandas as pd

# 1. 649 临床样本 ID
clin = pd.read_csv(r"D:\ONT\clinical_649.tsv", sep="\t")
clin_ids = clin["Sample_ID"].astype(str).tolist()
print("clinical 样本数:", len(clin_ids), " 去重后:", len(set(clin_ids)))

# 2. VCF 头里的样本名
vcf = r"E:\genotype_data\S957.1kGPont.merged.649clin_1027kgp.vcf"
with open(vcf) as f:
    for line in f:
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            break
vcf_samples = header[9:]
print("VCF 样本数:", len(vcf_samples))

# 3. 核对
clin_set = set(clin_ids)
in_vcf = [s for s in clin_ids if s in set(vcf_samples)]
missing = [s for s in clin_ids if s not in set(vcf_samples)]
print("649 临床样本中在 VCF 里找到:", len(in_vcf), " 缺失:", len(missing))
if missing:
    print("缺失样本:", missing[:20])

# 写样本名单(用于 bcftools view -S)
with open(r"D:\ONT\clinical_649.samples.txt", "w") as w:
    for s in clin_ids:
        w.write(s + "\n")
print("已写 D:\ONT\clinical_649.samples.txt")
