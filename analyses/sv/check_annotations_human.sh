#!/bin/bash
DIR=/mnt/d/ONT/AnnotSV_annotations/Annotations_Human
echo "=== Annotations_Human 顶层内容 ==="
ls -la "$DIR"
echo ""
echo "=== 各子目录大小 ==="
du -sh "$DIR"/*/ 2>/dev/null | sort -rh | head -30
