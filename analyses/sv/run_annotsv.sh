#!/bin/bash
export PATH=/home/administrator/miniconda3/envs/annotsv/bin:$PATH
echo "=== AnnotSV 开始 $(date) ==="
AnnotSV -SVinputFile /mnt/d/ONT/clinical_649.GRCh38.correct.vcf \
        -annotationsDir /mnt/d/ONT/AnnotSV_annotations \
        -genomeBuild GRCh38 \
        -outputFile /mnt/d/ONT/clinical_649.GRCh38.annotsv.tsv 2>&1 | tail -40
echo "=== AnnotSV 结束 $(date) ==="
echo "=== 输出文件 ==="
ls -la /mnt/d/ONT/clinical_649.GRCh38.annotsv.tsv 2>/dev/null || echo "输出文件未生成"
