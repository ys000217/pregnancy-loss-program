#!/bin/bash
ps aux | grep tclsh | grep -v grep | awk '{print $3"% CPU", $4"% MEM", "cpu_time="$10}'
echo "--- 输出文件时间 ---"
ls -la --time-style=+%H:%M:%S /mnt/d/ONT/clinical_649.GRCh38.annotsv.tsv.tmp /mnt/d/ONT/clinical_649.GRCh38.annotsv.SV_RE_intersect.tmp 2>/dev/null
echo "当前时间:"
date +%H:%M:%S
