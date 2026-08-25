#!/bin/bash
echo "=== 通过 API 取 INSTALL_annotations.sh ==="
curl -s --max-time 20 "https://api.github.com/repos/lgmgeo/AnnotSV/contents/bin/INSTALL_annotations.sh" 2>&1 | \
  grep -oE '"download_url": "[^"]*"' | head -2
# 直接下载
DL=$(curl -s --max-time 20 "https://api.github.com/repos/lgmgeo/AnnotSV/contents/bin/INSTALL_annotations.sh" | grep -oE 'https://raw[^"]*' | head -1)
echo "raw url: $DL"
curl -s --max-time 20 "$DL" 2>&1 | head -100
