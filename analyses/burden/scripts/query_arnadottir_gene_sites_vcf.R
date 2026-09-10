#!/usr/bin/env Rscript
# Tabix-query 648 VCF for Arnadottir gene-level ANNOVAR sites.

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
})

sites <- read.delim(
  "D:/ONT/analyses/burden/tables/external_arnadottir2025/arnadottir_gene_annovar_sites.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
pheno <- read.delim(
  "D:/ONT/analyses/burden/tables/suspect_abn/phenotype.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
vcf <- "D:/ONT/figure2/WGS_ONT_Intersection_648samples.vcf.gz"
tbx <- TabixFile(vcf)
hdr <- headerTabix(tbx)

con <- gzfile(vcf, "rt")
vcf_samples <- NULL
repeat {
  line <- readLines(con, n = 1L)
  if (!length(line)) break
  if (startsWith(line, "#CHROM")) {
    vcf_samples <- strsplit(line, "\t", fixed = TRUE)[[1]][-(1:9)]
    break
  }
}
close(con)

pheno$Condition <- tolower(pheno$Condition)
id_map <- setNames(pheno$Condition, pheno$VCF_Sample_ID)
id_map_l <- setNames(pheno$Condition, tolower(pheno$VCF_Sample_ID))

gt_has_alt <- function(gt) {
  gt <- sub(":.*", "", gt)
  if (gt %in% c(".", "./.", ".|.")) return(FALSE)
  parts <- unlist(strsplit(gt, "[/|]", perl = TRUE))
  any(parts != "0" & parts != ".")
}

gt_is_hom_alt <- function(gt) {
  gt <- sub(":.*", "", gt)
  parts <- unlist(strsplit(gt, "[/|]", perl = TRUE))
  length(parts) >= 2L && all(parts != "." & parts != "0")
}

allele_ok <- function(vcf_ref, vcf_alt, site_ref, site_alt) {
  vcf_ref <- toupper(vcf_ref)
  site_ref <- toupper(site_ref)
  site_alt <- toupper(site_alt)
  alts <- strsplit(toupper(vcf_alt), ",", fixed = TRUE)[[1]]
  if (vcf_ref == site_ref && site_alt %in% alts) return(TRUE)
  FALSE
}

n <- nrow(sites)
sites$in_vcf <- FALSE
sites$vcf_filter <- NA_character_
sites$n_carriers <- 0L
sites$n_hom <- 0L
sites$n_abnormal <- 0L
sites$n_normal <- 0L
sites$n_control <- 0L
sites$carrier_ids <- ""

carrier_long <- list()
k_long <- 0L

for (i in seq_len(n)) {
  chrom <- as.character(sites$chrom[i])
  pos <- as.integer(sites$pos[i])
  if (!(chrom %in% hdr$seqnames) && sub("^chr", "", chrom) %in% hdr$seqnames) {
    chrom <- sub("^chr", "", chrom)
  }
  win <- GRanges(chrom, IRanges(max(1L, pos - 2L), pos + 2L))
  hits <- tryCatch(scanTabix(tbx, param = win)[[1]], error = function(e) character())
  best <- NULL
  if (length(hits)) {
    for (line in hits) {
      f <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (as.integer(f[2]) == pos && allele_ok(f[4], f[5], sites$ref[i], sites$alt[i])) {
        best <- f
        break
      }
    }
  }
  if (is.null(best)) next
  sites$in_vcf[i] <- TRUE
  sites$vcf_filter[i] <- best[7]
  gts <- best[-(1:9)]
  ids <- character()
  n_hom <- 0L
  by <- c(abnormal = 0L, normal = 0L, control = 0L)
  for (j in seq_along(gts)) {
    if (!gt_has_alt(gts[j])) next
    sid <- vcf_samples[j]
    ids <- c(ids, sid)
    if (gt_is_hom_alt(gts[j])) n_hom <- n_hom + 1L
    cond <- id_map[[sid]]
    if (is.null(cond) || is.na(cond)) cond <- id_map_l[[tolower(sid)]]
    if (!is.null(cond) && !is.na(cond) && cond %in% names(by)) {
      by[[cond]] <- by[[cond]] + 1L
    }
    k_long <- k_long + 1L
    carrier_long[[k_long]] <- data.frame(
      canon_gene = sites$canon_gene[i],
      chrom = sites$chrom[i],
      pos = pos,
      ref = sites$ref[i],
      alt = sites$alt[i],
      consequence = sites$consequence[i],
      sample = sid,
      condition = ifelse(is.null(cond) || is.na(cond), "", cond),
      gt = sub(":.*", "", gts[j]),
      stringsAsFactors = FALSE
    )
  }
  sites$n_carriers[i] <- length(ids)
  sites$n_hom[i] <- n_hom
  sites$n_abnormal[i] <- by[["abnormal"]]
  sites$n_normal[i] <- by[["normal"]]
  sites$n_control[i] <- by[["control"]]
  sites$carrier_ids[i] <- paste(ids, collapse = ",")
  if (i %% 200L == 0L) message("queried ", i, "/", n)
}

out <- "D:/ONT/analyses/burden/tables/external_arnadottir2025"
write.table(sites, file.path(out, "arnadottir_gene_sites_in_648vcf.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (k_long) {
  write.table(do.call(rbind, carrier_long),
              file.path(out, "arnadottir_gene_carriers_long.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(
    data.frame(canon_gene = character(), chrom = character(), pos = integer(),
               ref = character(), alt = character(), consequence = character(),
               sample = character(), condition = character(), gt = character()),
    file.path(out, "arnadottir_gene_carriers_long.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
message("sites=", n, " in_vcf=", sum(sites$in_vcf), " with_carrier=", sum(sites$n_carriers > 0))
