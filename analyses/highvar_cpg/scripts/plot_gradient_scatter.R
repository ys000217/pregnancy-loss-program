#!/usr/bin/env Rscript
# Sample scatter plots (PCA + Pearson PCoA) at CpG-count gradients.
# No k-cut. Points colored by clinical Class3.

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "gradient_scatter")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MIN_SAMPLE_COVERAGE <- 0.90
FULL_TOP_N <- c(5000L, 10000L, 20000L, 50000L)
NC_TOP_N   <- c(2000L, 5000L, 10000L)

class_cols <- c(control = "#4DAF4A", normal_case = "#377EB8", abnormal_case = "#E41A1C")
class_pch  <- c(control = 16, normal_case = 17, abnormal_case = 15)

pearson_dist_samples <- function(mat) {
  cm <- cor(mat, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

impute_cpg_mean <- function(mat) {
  out <- mat
  na <- is.na(out)
  rs <- rowMeans(out, na.rm = TRUE)
  out[na] <- rs[row(out)[na]]
  out
}

coords_pca <- function(mat) {
  X <- t(impute_cpg_mean(mat))
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  data.table(
    sample = rownames(pc$x),
    X = pc$x[, 1],
    Y = pc$x[, 2],
    xlab = sprintf("PC1 (%.1f%%)", ve[1]),
    ylab = sprintf("PC2 (%.1f%%)", ve[2])
  )
}

coords_pcoa <- function(mat) {
  d <- pearson_dist_samples(mat)
  fit <- cmdscale(d, k = 2, eig = TRUE)
  eig <- fit$eig
  eig[eig < 0] <- 0
  ve <- 100 * eig[1:2] / sum(eig)
  xy <- fit$points
  data.table(
    sample = rownames(xy),
    X = xy[, 1],
    Y = xy[, 2],
    xlab = sprintf("PCoA1 (%.1f%%)", ve[1]),
    ylab = sprintf("PCoA2 (%.1f%%)", ve[2])
  )
}

plot_scatter <- function(dt, class3, main) {
  cls <- as.character(class3[dt$sample])
  cols <- class_cols[cls]
  pch  <- class_pch[cls]
  xlab <- dt$xlab[1]
  ylab <- dt$ylab[1]
  plot(dt$X, dt$Y, col = adjustcolor(cols, 0.75), pch = pch, cex = 0.85,
       xlab = xlab, ylab = ylab, main = main, las = 1)
  present <- unique(cls)
  present <- present[present %in% names(class_cols)]
  legend("topleft", legend = present, col = class_cols[present],
         pch = class_pch[present], bty = "n", pt.cex = 1.1, cex = 0.85)
}

save_one <- function(dt, class3, stem, panel_title) {
  fwrite(dt[, .(sample, Class3 = as.character(class3[sample]), X, Y, xlab, ylab)],
         file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
  png(file.path(OUT, paste0(stem, ".png")), width = 1600, height = 1400, res = 180)
  par(mar = c(4.2, 4.2, 3.2, 1.2))
  plot_scatter(dt, class3, panel_title)
  dev.off()
}

embed_gradient <- function(mat, class3, top_ns, var_ids, prefix, panel_name) {
  pca_list <- list()
  pcoa_list <- list()
  labels <- character()
  for (n in top_ns) {
    cpgs <- head(var_ids, min(n, length(var_ids)))
    message(sprintf("%s top%d: %d CpGs, %d samples", prefix, n, length(cpgs), ncol(mat)))
    sub <- mat[cpgs, , drop = FALSE]
    pca <- coords_pca(sub)
    pcoa <- coords_pcoa(sub)
    lab <- sprintf("%s | top %d", panel_name, length(cpgs))
    save_one(pca, class3, sprintf("%s_pca_top%d", prefix, n), paste(lab, "PCA"))
    save_one(pcoa, class3, sprintf("%s_pcoa_top%d", prefix, n), paste(lab, "PCoA (1-Pearson)"))
    pca_list[[length(pca_list) + 1L]] <- pca
    pcoa_list[[length(pcoa_list) + 1L]] <- pcoa
    labels <- c(labels, lab)
    rm(sub)
    gc()
  }
  list(pca = pca_list, pcoa = pcoa_list, labels = labels)
}

draw_grid <- function(coord_list, labels, class3, file_stub, method) {
  n <- length(coord_list)
  ncol <- if (n <= 3L) n else 2L
  nrow <- ceiling(n / ncol)
  pdf(file.path(OUT, paste0(file_stub, ".pdf")), width = 5.2 * ncol, height = 5 * nrow)
  par(mfrow = c(nrow, ncol), mar = c(4.2, 4.2, 3.2, 1.2))
  for (i in seq_len(n)) {
    plot_scatter(coord_list[[i]], class3, paste(labels[i], method))
  }
  dev.off()
  png(file.path(OUT, paste0(file_stub, ".png")),
      width = 900 * ncol, height = 860 * nrow, res = 140)
  par(mfrow = c(nrow, ncol), mar = c(4.2, 4.2, 3.2, 1.2))
  for (i in seq_len(n)) {
    plot_scatter(coord_list[[i]], class3, paste(labels[i], method))
  }
  dev.off()
}

# ---- Panel 1: full cohort ----
message("=== Full cohort ===")
cache <- readRDS(file.path(BASE, "clustering_output", "val_mat_samples.rds"))
val_mat <- cache$mat
sample_class3 <- cache$sample_class3
sample_na_rate <- colMeans(is.na(val_mat))
bad_samples <- names(sample_na_rate)[sample_na_rate > 0.20]
if (length(bad_samples)) {
  val_mat <- val_mat[, !(colnames(val_mat) %in% bad_samples), drop = FALSE]
  sample_class3 <- sample_class3[colnames(val_mat)]
}

var_dt <- fread(file.path(BASE, "CpG_variance.tsv"), header = FALSE,
                col.names = c("CpG_ID", "count", "variance"))
setorder(var_dt, -variance)
var_dt <- var_dt[CpG_ID %in% rownames(val_mat)]
min_non_na <- ceiling(ncol(val_mat) * MIN_SAMPLE_COVERAGE)
coverage_ok <- rowSums(!is.na(val_mat)) >= min_non_na
var_ids_ok <- var_dt[CpG_ID %in% names(coverage_ok)[coverage_ok], CpG_ID]

full <- embed_gradient(val_mat, sample_class3, FULL_TOP_N, var_ids_ok,
                       "full", "Full cohort")
draw_grid(full$pca, full$labels, sample_class3, "full_pca_gradient", "PCA")
draw_grid(full$pcoa, full$labels, sample_class3, "full_pcoa_gradient", "PCoA")

rm(val_mat, cache)
gc()

# ---- Panel 2: NC ----
message("=== NC (abnormal removed) ===")
nc_rds <- file.path(BASE, "clustering_output", "gradient_heatmaps", "nc_mat_cache.rds")
if (!file.exists(nc_rds)) stop("Missing NC cache: ", nc_rds)
nc <- readRDS(nc_rds)
nc_mat <- nc$mat
nc_cls <- nc$class3

var_nc <- fread(file.path(BASE, "CpG_variance_NC.tsv"), header = FALSE,
                col.names = c("CpG_ID", "count", "variance"))
setorder(var_nc, -variance)
var_nc_ids <- var_nc[CpG_ID %in% rownames(nc_mat), CpG_ID]
if (length(var_nc_ids) < nrow(nc_mat)) {
  var_nc_ids <- c(var_nc_ids, setdiff(rownames(nc_mat), var_nc_ids))
}

nc_emb <- embed_gradient(nc_mat, nc_cls, NC_TOP_N, var_nc_ids, "NC", "NC (no abnormal)")
draw_grid(nc_emb$pca, nc_emb$labels, nc_cls, "NC_pca_gradient", "PCA")
draw_grid(nc_emb$pcoa, nc_emb$labels, nc_cls, "NC_pcoa_gradient", "PCoA")

writeLines(c(
  "Sample scatter plots at CpG-count gradients (no k).",
  "",
  "PCA: samples x CpGs after per-CpG mean fill of NA (visualization only).",
  "PCoA: classical MDS of 1 - Pearson (same distance as Ward clustering).",
  "Color/shape = Class3, not unsupervised cluster IDs.",
  "",
  "full_pca_gradient / full_pcoa_gradient: all samples, cohort variance 5k/10k/20k/50k.",
  "NC_pca_gradient / NC_pcoa_gradient: abnormal removed, NC variance 2k/5k/10k.",
  "Per-panel coordinates: *_pca_topN.tsv and *_pcoa_topN.tsv."
), file.path(OUT, "README.txt"))

message("=== Done ===")
message("Output: ", OUT)
