#!/bin/bash
SHARE=/home/administrator/miniconda3/envs/annotsv/share/tcl8.6/AnnotSV
echo "=== AnnotSV-config.tcl 里 annotationsDir 相关 ==="
grep -nE 'annotationsDir|Annotations_Human|Annotations_Users' $SHARE/AnnotSV-config.tcl | head -20
