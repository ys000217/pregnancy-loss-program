#!/bin/bash
echo "=== README.md 全文(前60行) ==="
curl -s --max-time 20 https://raw.githubusercontent.com/lgmgeo/AnnotSV/master/README.md 2>&1 | head -60
echo ""
echo "=== bin 目录 ==="
curl -s --max-time 20 https://api.github.com/repos/lgmgeo/AnnotSV/contents/bin 2>&1 | grep '"name"'
