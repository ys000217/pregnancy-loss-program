#!/bin/bash
set -u
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
IN=/mnt/d/ONT/clinical_649.vcf
CHAIN=/mnt/e/genotype_data/liftover/CN1_to_GRCh38.chain
FAI=/mnt/e/genotype_data/liftover/GRCh38_chr.fasta.fai
OUT=/mnt/d/ONT/clinical_649.GRCh38.correct.vcf

echo "=== A. bcftools 提取断点表 ==="
$TOOLS/bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/CHR2\n' $IN > /tmp/sv_table.tsv
echo "SV 记录数:"; wc -l < /tmp/sv_table.tsv
echo "ID 为 '.' 的记录数:"; awk -F'\t' '$1=="."' /tmp/sv_table.tsv | wc -l

echo "=== B. 转成断点 BED(CN1, 每断点一行) ==="
awk -F'\t' '{
  id=$1; chr=$2; pos=$3+0; end=$4+0; type=$5; chr2=$6
  if (type=="INS") print chr"\t"pos-1"\t"pos"\t"id"\tINS\tS"
  else if (type=="TRA"||type=="BND") {
    print chr"\t"pos-1"\t"pos"\t"id"\t"type"\tL"
    print chr2"\t"end-1"\t"end"\t"id"\t"type"\tR"
  } else {
    print chr"\t"pos-1"\t"pos"\t"id"\t"type"\tL"
    print chr"\t"end-1"\t"end"\t"id"\t"type"\tR"
  }
}' /tmp/sv_table.tsv > /tmp/cn1_bp.bed
echo "断点数:"; wc -l < /tmp/cn1_bp.bed

echo "=== C. CrossMap bed 抬升 CN1->GRCh38 ==="
$TOOLS/CrossMap bed $CHAIN /tmp/cn1_bp.bed /tmp/grch38_bp.bed 2>&1 | tail -2
echo "抬升后断点数:"; wc -l < /tmp/grch38_bp.bed
echo "unmap 断点数:"; wc -l < /tmp/grch38_bp.bed.unmap 2>/dev/null || echo 0

echo "=== D. awk 重建 VCF(把抬升断点回填 POS/END/CHR2) ==="
awk -F'\t' 'BEGIN{OFS="\t"}
  NR==FNR { id=$4; side=$6; if($2<0||$3<=0) next; c[id SUBSEP side]=$1; p[id SUBSEP side]=$3; next }
  /^#/ { next }
  {
    id=$3
    svtype=""; n=split($8,f,";")
    for(i=1;i<=n;i++) if(f[i]~/^SVTYPE=/) svtype=substr(f[i],8)
    if (svtype=="INS") {
      k=id SUBSEP "S"
      if (k in p) { $1=c[k]; $2=p[k]; sub(/END=[0-9]+/,"END="p[k],$8); sub(/CHR2=[^;]+/,"CHR2="c[k],$8); print }
    } else {
      kL=id SUBSEP "L"; kR=id SUBSEP "R"
      if (kL in p && kR in p) {
        if (svtype=="TRA"||svtype=="BND"||c[kL]==c[kR]) {
          $1=c[kL]; $2=(p[kL]<=p[kR]?p[kL]:p[kR])
          sub(/END=[0-9]+/,"END="(p[kL]<=p[kR]?p[kR]:p[kL]),$8)
          sub(/CHR2=[^;]+/,"CHR2="c[kR],$8)
          print
        }
      }
    }
  }' /tmp/grch38_bp.bed $IN > /tmp/body.vcf
echo "重建出的记录数:"; wc -l < /tmp/body.vcf

echo "=== E. 头部换成 GRCh38 contig ==="
{
  grep '^##' $IN | grep -v '^##contig'
  awk '{printf "##contig=<ID=%s,length=%d,assembly=GRCh38_chr.fasta>\n",$1,$2}' $FAI
  echo "##reference=GRCh38.p14"
  echo "##liftover=CN1_to_GRCh38_CrossMap_bed"
  grep '^#CHROM' $IN
} > /tmp/header.vcf
cat /tmp/header.vcf /tmp/body.vcf > $OUT
echo "=== 完成 ==="
echo "最终输出记录数:"; grep -vc '^#' $OUT
echo "输出文件:"; ls -la $OUT
