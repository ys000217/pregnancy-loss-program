#!/usr/bin/env Rscript
# Overlap Arnadottir all DNMs (no pathogenicity filter) vs 648 VCF alleles.

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(GenomicRanges)
})

dnm <- fread(
  "D:/ONT/analyses/burden/tables/external_arnadottir2025/dnms.tsv",
  sep = "\t"
)
dnm[, Pos := as.integer(Pos)]
u <- unique(dnm[, .(Chr, Pos, Ref, Alt)])
setkey(u, Chr, Pos)
message("unique DNM sites=", nrow(u), " rows=", nrow(dnm))

vcf <- "D:/ONT/figure2/WGS_ONT_Intersection_648samples.vcf.gz"
tbx <- TabixFile(vcf)
sn <- headerTabix(tbx)$seqnames

parse_hits <- function(lines, pos, ref, alt) {
  pos_hit <- FALSE
  allele_hit <- FALSE
  vcf_alleles <- character()
  if (!length(lines)) {
    return(list(pos_hit = FALSE, allele_hit = FALSE, vcf_alleles = ""))
  }
  ref <- toupper(ref)
  alt <- toupper(alt)
  for (line in lines) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    vp <- as.integer(f[2])
    if (vp != pos) next
    pos_hit <- TRUE
    vr <- toupper(f[4])
    va <- strsplit(toupper(f[5]), ",", fixed = TRUE)[[1]]
    vcf_alleles <- c(vcf_alleles, paste0(vr, ">", paste(va, collapse = ",")))
    if (vr == ref && alt %in% va) allele_hit <- TRUE
    # indel: same event if deleted/inserted seq matches ignoring pad base
    if (!allele_hit && nchar(ref) != nchar(alt)) {
      if (vr == ref && alt %in% va) allele_hit <- TRUE
    }
  }
  list(
    pos_hit = pos_hit,
    allele_hit = allele_hit,
    vcf_alleles = paste(unique(vcf_alleles), collapse = ";")
  )
}

chroms <- unique(u$Chr)
out <- vector("list", length(chroms))
n_pos <- 0L
n_all <- 0L
k <- 0L

for (chr in chroms) {
  k <- k + 1L
  qchr <- chr
  if (!(qchr %in% sn) && sub("^chr", "", qchr) %in% sn) {
    qchr <- sub("^chr", "", qchr)
  }
  sub <- u[Chr == chr]
  gr <- GRanges(qchr, IRanges(sub$Pos, sub$Pos))
  hits <- tryCatch(scanTabix(tbx, param = gr), error = function(e) {
    message("ERR ", chr, " ", conditionMessage(e))
    vector("list", nrow(sub))
  })
  if (length(hits) != nrow(sub)) {
    message("length mismatch ", chr, " ", length(hits), " vs ", nrow(sub))
  }
  pos_hit <- logical(nrow(sub))
  allele_hit <- logical(nrow(sub))
  vcf_al <- character(nrow(sub))
  for (i in seq_len(nrow(sub))) {
    h <- if (i <= length(hits)) hits[[i]] else character()
    r <- parse_hits(h, sub$Pos[i], sub$Ref[i], sub$Alt[i])
    pos_hit[i] <- r$pos_hit
    allele_hit[i] <- r$allele_hit
    vcf_al[i] <- r$vcf_alleles
  }
  sub[, `:=`(pos_overlap = pos_hit, allele_overlap = allele_hit, vcf_alleles = vcf_al)]
  out[[k]] <- sub
  n_pos <- n_pos + sum(pos_hit)
  n_all <- n_all + sum(allele_hit)
  message(chr, " sites=", nrow(sub), " pos_overlap=", sum(pos_hit), " allele_overlap=", sum(allele_hit))
}

res <- rbindlist(out)
od <- "D:/ONT/analyses/burden/tables/external_arnadottir2025"
fwrite(res[allele_overlap == TRUE], file.path(od, "dnm_allele_overlap_648vcf.tsv"), sep = "\t")
fwrite(res[pos_overlap == TRUE & allele_overlap == FALSE],
       file.path(od, "dnm_pos_only_overlap_648vcf.tsv"), sep = "\t")

summ <- data.table(
  unique_dnm_sites = nrow(res),
  snv = sum(nchar(res$Ref) == 1L & nchar(res$Alt) == 1L),
  indel = sum(!(nchar(res$Ref) == 1L & nchar(res$Alt) == 1L)),
  pos_overlap = sum(res$pos_overlap),
  allele_overlap = sum(res$allele_overlap),
  pos_only = sum(res$pos_overlap & !res$allele_overlap)
)
fwrite(summ, file.path(od, "dnm_overlap_648vcf_summary.tsv"), sep = "\t")
print(summ)
message("done")
