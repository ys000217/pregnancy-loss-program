#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""兼容入口：转调 pathogenic_element_enrichment.py（新富集逻辑）。

旧逻辑用「显著 vs 全部非显著」做 Fisher，会把海量非显著 pair 塞进背景。
新逻辑：
  - 富集 p 值只在显著 pair 内、相对 Roadmap 基因组覆盖检验；
  - 另随机抽样非显著 pair 仅作分布对照。
"""
from pathogenic_element_enrichment import main

if __name__ == "__main__":
    main()
