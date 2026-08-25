#!/bin/bash
echo "=== GitHub README (数据下载部分) ==="
curl -s --max-time 20 https://raw.githubusercontent.com/lgmgeo/AnnotSV/master/README.md 2>&1 | grep -iE 'download|data|makeAnnotSV|wget|curl|\.tar|\.gz' | head -30
echo ""
echo "=== repo 根目录文件 ==="
curl -s --max-time 20 https://api.github.com/repos/lgmgeo/AnnotSV/contents/ 2>&1 | grep '"name"' | head -30
