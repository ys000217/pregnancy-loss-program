#!/usr/bin/env Rscript
# Query 648-sample hg38 VCF for Arnadottir pathogenic SSVs.

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
})

sites <- read.delim(
  "D:/ONT/analyses/burden/tables/external_arnadottir2025/causal_ssv_grch38.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
pheno <- read.delim(
  "D:/ONT/figure2/sample_phenotype_648.tsv",
  stringsAsFactors = FALSE, check.names = FALSE
)
vcf <- "D:/ONT/figure2/WGS_ONT_Intersection_648samples.vcf.gz"
tbx <- TabixFile(vcf)

hdr <- headerTabix(tbx)
message("tabix seqnames sample: ", paste(head(hdr$seqnames, 5), collapse = ","))

# VCF sample IDs
con <- gzfile(vcf, "rt")
vcf_samples <- NULL
repeat {
  line <- readLines(con, n = 1L)
  if (length(line) == 0L) break
  if (startsWith(line, "#CHROM")) {
    vcf_samples <- strsplit(line, "\t", fixed = TRUE)[[1]][-(1:9)]
    break
  }
}
close(con)
message("VCF samples: ", length(vcf_samples))

pheno$Condition <- tolower(pheno$Condition)
id_map <- setNames(pheno$Condition, pheno$VCF_Sample_ID)
# tolerate case differences
id_map_l <- setNames(pheno$Condition, tolower(pheno$VCF_Sample_ID))

gt_allele_n <- function(gt) {
  gt <- sub(":.*", "", gt)
  if (gt %in% c(".", "./.", ".|.")) return(NA_integer_)
  parts <- unlist(strsplit(gt, "[/|]", perl = TRUE))
  if (any(parts == ".")) {
    # still count known alts
    parts <- parts[parts != "."]
    if (!length(parts)) return(NA_integer_)
  }
  sum(as.integer(parts) > 0L)
}

norm_seq <- function(x) toupper(as.character(x))

allele_match <- function(vcf_ref, vcf_alt, site_ref, site_alt) {
  vcf_ref <- norm_seq(vcf_ref)
  site_ref <- norm_seq(site_ref)
  site_alt <- norm_seq(site_alt)
  alts <- strsplit(norm_seq(vcf_alt), ",", fixed = TRUE)[[1]]
  # exact
  if (vcf_ref == site_ref && site_alt %in% alts && nzchar(site_ref) && nzchar(site_alt)) {
    return(which(alts == site_alt)[1])
  }
  # deletion SPDI: ref=deleted, alt empty
  if (nzchar(site_ref) && !nzchar(site_alt)) {
    for (i in seq_along(alts)) {
      a <- alts[i]
      if (nchar(vcf_ref) > nchar(a) && grepl(site_ref, vcf_ref, fixed = TRUE) &&
          (vcf_ref == paste0(a, site_ref) || vcf_ref == paste0(site_ref, a) ||
           grepl(paste0(a, site_ref), vcf_ref, fixed = TRUE))) {
        return(i)
      }
    }
  }
  # insertion SPDI: ref empty, alt=inserted
  if (!nzchar(site_ref) && nzchar(site_alt)) {
    for (i in seq_along(alts)) {
      a <- alts[i]
      if (nchar(a) > nchar(vcf_ref) && grepl(site_alt, a, fixed = TRUE)) {
        return(i)
      }
    }
  }
  NA_integer_
}

out_rows <- list()
carrier_rows <- list()

for (i in seq_len(nrow(sites))) {
  s <- sites[i, ]
  rec <- data.frame(
    gene = s$gene,
    allele_n = s$allele_n,
    hgvs_c = s$hgvs_c,
    genotype_class = s$genotype_class,
    karyotype = s$karyotype,
    euploid_likely_causal = s$euploid_likely_causal,
    chrom = s$chrom,
    pos = s$pos,
    ref = s$ref,
    alt = s$alt,
    mapping_error = s$error,
    queried = FALSE,
    vcf_chrom = NA_character_,
    vcf_pos = NA_integer_,
    vcf_ref = NA_character_,
    vcf_alt = NA_character_,
    n_alt_carriers = 0L,
    n_hom_alt = 0L,
    n_abnormal = 0L,
    n_normal = 0L,
    n_control = 0L,
    carrier_ids = "",
    stringsAsFactors = FALSE
  )
  if (!nzchar(as.character(s$chrom)) || is.na(s$pos) || s$pos == "") {
    out_rows[[i]] <- rec
    next
  }
  pos <- as.integer(s$pos)
  chrom <- as.character(s$chrom)
  if (!(chrom %in% hdr$seqnames) && sub("^chr", "", chrom) %in% hdr$seqnames) {
    chrom <- sub("^chr", "", chrom)
  }
  win <- GRanges(chrom, IRanges(max(1L, pos - 25L), pos + 25L))
  rec$queried <- TRUE
  hits <- tryCatch(scanTabix(tbx, param = win)[[1]], error = function(e) character())
  best <- NULL
  best_ai <- NA_integer_
  if (length(hits)) {
    for (line in hits) {
      f <- strsplit(line, "\t", fixed = TRUE)[[1]]
      ai <- allele_match(f[4], f[5], s$ref, s$alt)
      if (!is.na(ai)) {
        best <- f
        best_ai <- ai
        break
      }
    }
  }
  if (is.null(best)) {
    out_rows[[i]] <- rec
    next
  }
  rec$vcf_chrom <- best[1]
  rec$vcf_pos <- as.integer(best[2])
  rec$vcf_ref <- best[4]
  rec$vcf_alt <- best[5]
  gts <- best[-(1:9)]
  n_alt <- 0L
  n_hom <- 0L
  ids <- character()
  by_g <- c(abnormal = 0L, normal = 0L, control = 0L)
  for (j in seq_along(gts)) {
    n <- gt_allele_n(gts[j])
    if (is.na(n) || n <= 0L) next
    # if multi-allelic, crude: any alt counted; exact alt index not decoded from GT unless == best_ai
    gt <- sub(":.*", "", gts[j])
    parts <- unlist(strsplit(gt, "[/|]", perl = TRUE))
    parts <- suppressWarnings(as.integer(parts[parts != "."]))
    if (!length(parts) || !any(parts == best_ai)) {
      # biallelic common case: alt index 1
      if (!(length(strsplit(best[5], ",", fixed = TRUE)[[1]]) == 1L && any(parts > 0L))) {
        next
      }
    }
    n_alt <- n_alt + 1L
    if (length(parts) >= 2L && all(parts == best_ai)) n_hom <- n_hom + 1L
    sid <- vcf_samples[j]
    ids <- c(ids, sid)
    cond <- id_map[[sid]]
    if (is.null(cond) || is.na(cond)) cond <- id_map_l[[tolower(sid)]]
    if (!is.null(cond) && !is.na(cond) && cond %in% names(by_g)) {
      by_g[[cond]] <- by_g[[cond]] + 1L
    }
    carrier_rows[[length(carrier_rows) + 1L]] <- data.frame(
      gene = s$gene, hgvs_c = s$hgvs_c, sample = sid, condition = ifelse(is.null(cond), NA, cond),
      gt = gt, stringsAsFactors = FALSE
    )
  }
  rec$n_alt_carriers <- n_alt
  rec$n_hom_alt <- n_hom
  rec$n_abnormal <- by_g[["abnormal"]]
  rec$n_normal <- by_g[["normal"]]
  rec$n_control <- by_g[["control"]]
  rec$carrier_ids <- paste(ids, collapse = ",")
  out_rows[[i]] <- rec
  message(s$gene, " ", s$hgvs_c, " carriers=", n_alt)
}

res <- do.call(rbind, out_rows)
out_dir <- "D:/ONT/analyses/burden/tables/external_arnadottir2025"
write.table(res, file.path(out_dir, "causal_ssv_in_648vcf.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (length(carrier_rows)) {
  write.table(do.call(rbind, carrier_rows), file.path(out_dir, "causal_ssv_648_carriers.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(
    data.frame(gene = character(), hgvs_c = character(), sample = character(),
               condition = character(), gt = character()),
    file.path(out_dir, "causal_ssv_648_carriers.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
message("done rows=", nrow(res), " with_carriers=", sum(res$n_alt_carriers > 0))
