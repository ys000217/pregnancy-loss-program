#!/bin/bash
echo "=== INSTALL_annotations.sh 内容(通过API base64解码) ==="
curl -s --max-time 25 "https://api.github.com/repos/lgmgeo/AnnotSV/contents/bin/INSTALL_annotations.sh" | \
  /home/administrator/miniconda3/bin/python3 -c "import sys,json,base64; d=json.load(sys.stdin); print(base64.b64decode(d['content']).decode())" 2>&1 | head -120
