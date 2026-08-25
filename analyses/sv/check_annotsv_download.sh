#!/bin/bash
echo "=== lbgi.fr AnnotSV 下载页 ==="
curl -s --max-time 20 https://lbgi.fr/AnnotSV/downloads/ 2>&1 | grep -oiE 'href="[^"]*"' | grep -iE 'tar|zip|data|annotsv' | head -20
echo "=== GitHub releases ==="
curl -s --max-time 20 https://api.github.com/repos/lgmgeo/AnnotSV/releases/latest 2>&1 | grep -E 'tag_name|browser_download_url' | head -20
