#!/bin/bash
CONDA=/home/administrator/miniconda3/bin/conda
echo "=== 安装 tcl(tclsh) ==="
$CONDA install -y -n annotsv -c conda-forge tcl 2>&1 | tail -5
echo "=== 验证 tclsh ==="
ls /home/administrator/miniconda3/envs/annotsv/bin/tclsh* 2>/dev/null
echo "=== AnnotSV -version ==="
/home/administrator/miniconda3/envs/annotsv/bin/AnnotSV -version 2>&1 | head -5
echo "=== AnnotSV -help(前30行) ==="
/home/administrator/miniconda3/envs/annotsv/bin/AnnotSV -help 2>&1 | head -30
