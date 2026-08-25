#!/bin/bash
echo "=== 网络检查 ==="
curl -sI --max-time 15 https://github.com/lgmgeo/AnnotSV 2>&1 | head -1
curl -sI --max-time 15 https://anaconda.org/bioconda/annotsv 2>&1 | head -1
echo "=== conda 版本 ==="
/home/administrator/miniconda3/bin/conda --version 2>&1
echo "=== perl 是否可用 ==="
command -v perl && perl -v 2>&1 | head -2
echo "=== 是否已装 annotsv ==="
ls /home/administrator/miniconda3/envs/*/bin/AnnotSV 2>/dev/null || echo "未安装"
echo "=== bioconda 已配置? ==="
/home/administrator/miniconda3/bin/conda config --show channels 2>/dev/null | head -5
