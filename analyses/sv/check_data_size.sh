#!/bin/bash
URL="https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz"
echo "=== 检查注释数据大小 ==="
curl -sI --max-time 30 "$URL" 2>&1 | grep -iE 'content-length|content-type|HTTP'
