#!/bin/bash
set -e
DIR=/mnt/d/ONT/AnnotSV_annotations
mkdir -p $DIR
cd $DIR
echo "=== 下载 Annotations_Human_3.5.tar.gz (5.3GB) $(date) ==="
curl -C - -L -o Annotations_Human_3.5.tar.gz \
  "https://cstb-icube.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz"
echo "=== 下载完成 $(date), 大小: ==="
ls -la Annotations_Human_3.5.tar.gz
echo "=== 解压 ==="
tar -xf Annotations_Human_3.5.tar.gz
rm -f Annotations_Human_3.5.tar.gz
echo "=== 解压完成, 目录结构 ==="
ls -la $DIR
