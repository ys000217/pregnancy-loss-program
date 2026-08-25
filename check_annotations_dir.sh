#!/bin/bash
DIR=/mnt/d/ONT/AnnotSV_annotations
echo "=== 目录是否存在 ==="
if [ -d "$DIR" ]; then
  echo "存在: $DIR"
  echo ""
  echo "=== 目录内容 ==="
  ls -la "$DIR"
  echo ""
  echo "=== 总大小 ==="
  du -sh "$DIR" 2>/dev/null
  echo ""
  echo "=== 是否有完整解压的 Annotations_Human ==="
  ls -d "$DIR"/Annotations_Human* 2>/dev/null || echo "无(未解压)"
else
  echo "目录不存在: $DIR"
fi
