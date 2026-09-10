#!/usr/bin/env Rscript
# Export symbol -> canonical CDS length (bp) from TxDb knownGene.
suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicFeatures)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else
  "D:/ONT/analyses/burden/tables/loose_dual_filter/gene_cds_length_bp.tsv"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx_lens <- transcriptLengths(txdb, with.cds_len = TRUE)
tx_lens <- tx_lens[!is.na(tx_lens$gene_id) & tx_lens$cds_len > 0, ]
tx_lens <- tx_lens[order(tx_lens$cds_len, decreasing = TRUE), ]
tx_lens <- tx_lens[!duplicated(tx_lens$gene_id), ]

eg <- as.character(tx_lens$gene_id)
sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = eg,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
sym <- unname(sym[eg])

df <- data.frame(
  gene = sym,
  entrez = eg,
  cds_bp = as.integer(tx_lens$cds_len),
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$gene) & nzchar(df$gene), ]
# if duplicate symbols, keep longest cds
df <- df[order(df$cds_bp, decreasing = TRUE), ]
df <- df[!duplicated(df$gene), ]
df <- df[order(df$gene), ]

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(df), " genes -> ", out)
