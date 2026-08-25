#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/liftover/bin
echo "=== 从输出(GRCh38)取前3个 DEL 的 ID + POS/END ==="
$TOOLS/bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
  /mnt/d/ONT/clinical_649.GRCh38.vcf 2>/dev/null | awk -F'\t' '$5=="DEL"' | head -3 | tee /tmp/del_out.txt

echo ""
echo "=== 回查输入(CN1)中相同 ID 的 POS/END ==="
while IFS=$'\t' read -r id chrom pos end svt; do
  echo "--- ID=$id (输出: $chrom $pos-$end) ---"
  $TOOLS/bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
    /mnt/d/ONT/clinical_649.vcf 2>/dev/null | grep -F "$id" | head -1 | \
    awk -F'\t' '{print "    输入(CN1): "$1" "$2"-"$3}'
done < /tmp/del_out.txt
