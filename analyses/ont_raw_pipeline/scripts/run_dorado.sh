#!/bin/bash

GLIBC_LIB=/public/software/tools/glibc2.35/lib
GCC_LIB=/public/software/compiler/gcc/11.2.0/lib64
SYS_LIB=/usr/lib64
DORADO_BIN=/data/home/4250160/dorado-1.1.1-linux-x64/bin/dorado

LD_PRELOAD=$GCC_LIB/libstdc++.so.6 \
$GLIBC_LIB/ld-linux-x86-64.so.2 \
  --library-path $GLIBC_LIB:$SYS_LIB:$GCC_LIB \
  $DORADO_BIN "$@"

