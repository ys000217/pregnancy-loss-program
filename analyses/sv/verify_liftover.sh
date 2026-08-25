#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
echo "=== 输出 VCF 头部 contig(前3条, 应=GRCh38长度) ==="
grep '^##contig' /mnt/d/ONT/clinical_649.GRCh38.vcf | head -3

echo ""
echo "=== 取输入中前2个 DEL 的 ID + POS/END(CN1) ==="
$TOOLS/bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
  /mnt/d/ONT/clinical_649.vcf 2>/dev/null | awk -F'\t' '$5=="DEL"' | head -2 | tee /tmp/del_input.txt

echo ""
echo "=== 输出中对应 ID 的 POS/END(GRCh38) ==="
while IFS=$'\t' read -r id chrom pos end svt; do
  echo "--- ID=$id ---"
  $TOOLS/bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
    /mnt/d/ONT/clinical_649.GRCh38.vcf 2>/dev/null | grep -F "$id" | head -1
done < /tmp/del_input.txt
