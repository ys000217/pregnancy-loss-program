#!/bin/bash
URL="https://cstb-icube.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz"
echo "=== 最终 URL 的 Content-Length ==="
curl -sIL --max-time 30 "$URL" 2>&1 | grep -iE 'content-length|content-type|HTTP/' | tail -5
