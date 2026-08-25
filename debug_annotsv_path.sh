#!/bin/bash
BIN=/home/administrator/miniconda3/envs/annotsv/bin
SHARE=/home/administrator/miniconda3/envs/annotsv/share
echo "=== share 下的 tcl 目录 ==="
ls -d $SHARE/tcl* 2>/dev/null
echo "=== tcl8.6/AnnotSV 内容(前10) ==="
ls $SHARE/tcl8.6/AnnotSV/ 2>/dev/null | head -10
echo "=== AnnotSV 脚本里 tclDir 相关行 ==="
grep -n 'tclDir\|share/tcl\|tcl8' $BIN/AnnotSV | head -20
echo "=== AnnotSV 脚本 50-65 行 ==="
sed -n '50,65p' $BIN/AnnotSV
