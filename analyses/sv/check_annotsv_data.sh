#!/bin/bash
TOOLS=/home/administrator/miniconda3/envs/annotsv/bin
echo "=== AnnotSV 版本 ==="
$TOOLS/AnnotSV -version 2>&1 | head -3
echo "=== annotsv 环境里的 AnnotSV 相关目录 ==="
find /home/administrator/miniconda3/envs/annotsv -maxdepth 3 -iname '*annotsv*' 2>/dev/null
echo "=== share 目录 ==="
ls /home/administrator/miniconda3/envs/annotsv/share/ 2>/dev/null
echo "=== AnnotSV 数据目录(Annotations_Human) ==="
find /home/administrator/miniconda3/envs/annotsv -maxdepth 4 -iname 'Annotations_Human' 2>/dev/null
ls -d /home/administrator/miniconda3/envs/annotsv/share/AnnotSV/Annotations_Human 2>/dev/null && echo "有内置数据"
