#!/bin/bash
set -e
CONDA=/home/administrator/miniconda3/bin/conda
echo "=== 创建 annotsv 环境并安装 ==="
$CONDA create -y -n annotsv -c conda-forge -c bioconda annotsv 2>&1 | tail -30
echo "=== 验证 ==="
ls /home/administrator/miniconda3/envs/annotsv/bin/AnnotSV && echo "AnnotSV 安装成功"
