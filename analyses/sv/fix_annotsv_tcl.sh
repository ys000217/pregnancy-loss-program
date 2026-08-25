#!/bin/bash
SHARE=/home/administrator/miniconda3/envs/annotsv/share
mkdir -p $SHARE/tcl
ln -sfn $SHARE/tcl8.6/AnnotSV $SHARE/tcl/AnnotSV
echo "=== 软链接已建, 验证 ==="
ls -ld $SHARE/tcl/AnnotSV
export PATH=/home/administrator/miniconda3/envs/annotsv/bin:$PATH
echo "=== AnnotSV -version ==="
AnnotSV -version 2>&1 | head -5
echo "=== 数据目录(share/AnnotSV) 现状 ==="
ls -la $SHARE/AnnotSV/
