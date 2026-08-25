#!/bin/bash
BIN=/home/administrator/miniconda3/envs/annotsv/bin
echo "=== bin 里的 tcl 相关 ==="
ls $BIN | grep -i tcl
echo "=== AnnotSV 脚本 shebang ==="
head -3 $BIN/AnnotSV
echo "=== AnnotSV 数据目录 ==="
ls /home/administrator/miniconda3/envs/annotsv/share/AnnotSV/ 2>/dev/null
echo "=== AnnotSV 配置里的数据路径 ==="
find /home/administrator/miniconda3/envs/annotsv/etc/AnnotSV -type f 2>/dev/null | head
cat /home/administrator/miniconda3/envs/annotsv/etc/AnnotSV/* 2>/dev/null | head -20
