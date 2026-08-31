#!/usr/bin/env Rscript
# k-scan (global top10000) + normal_case sub-clustering (within-cohort variance)
# No imputation; pairwise Pearson distance; Ward.D2

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output")
SUB  <- file.path(OUT, "kscan_and_normal_sub")
dir.create(SUB, recursive = TRUE, showWarnings = FALSE)

PRIMARY_TOP_N <- 10000L
MIN_SAMPLE_COVERAGE <- 0.90
K_RANGE <- 2:10

message("=== Load data ===")
clin <- fread("D:/ONT/clinical_649.tsv")
clin[, Sample_ID := as.character(Sample_ID)]
clin[, Class3 := fifelse(tolower(Group1) == "control", "control",
                 fifelse(tolower(Group4) == "abnormal", "abnormal_case", "normal_case"))]

cache <- readRDS(file.path(OUT, "val_mat_samples.rds"))
val_mat <- cache$mat
sample_class3 <- cache$sample_class3

sample_na_rate <- colMeans(is.na(val_mat))
bad <- names(sample_na_rate)[sample_na_rate > 0.20]
if (length(bad)) {
  val_mat <- val_mat[, !(colnames(val_mat) %in% bad), drop = FALSE]
  sample_class3 <- sample_class3[colnames(val_mat)]
}

clin_key <- clin[, .(key = tolower(Sample_ID), Group2, Group1, Group4, Class3)]
setkey(clin_key, key)
meta <- data.table(
  sample = colnames(val_mat),
  Class3 = sample_class3,
  Group2 = clin_key[J(tolower(colnames(val_mat))), Group2]
)

var_dt <- fread(file.path(BASE, "CpG_variance.tsv"), header = FALSE,
                col.names = c("CpG_ID", "count", "variance"))
setorder(var_dt, -variance)
var_dt <- var_dt[CpG_ID %in% rownames(val_mat)]

min_non_na <- ceiling(ncol(val_mat) * MIN_SAMPLE_COVERAGE)
var_ids_ok <- var_dt[CpG_ID %in% rownames(val_mat)[rowSums(!is.na(val_mat)) >= min_non_na], CpG_ID]

pearson_dist_samples <- function(mat) {
  cm <- cor(mat, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

fowlkes_mallows <- function(lab1, lab2) {
  # External validation only: cluster partition vs known reference labels.
  lab1 <- as.integer(factor(lab1))
  lab2 <- as.integer(factor(lab2))
  tab <- table(lab1, lab2)
  tk <- sum(choose(tab, 2))
  pk <- sum(choose(rowSums(tab), 2))
  qk <- sum(choose(colSums(tab), 2))
  if (pk == 0 || qk == 0) return(NA_real_)
  sqrt(tk / pk * tk / qk)
}

row_var <- function(m) {
  apply(m, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::var(x)
  })
}

pick_top_var <- function(mat, ids_ordered, n, min_cov_frac) {
  min_n <- ceiling(ncol(mat) * min_cov_frac)
  ok <- rownames(mat)[rowSums(!is.na(mat)) >= min_n]
  head(ids_ordered[ids_ordered %in% ok], n)
}

cluster_at_k <- function(mat, cpgs, k) {
  sub <- mat[cpgs, , drop = FALSE]
  d <- pearson_dist_samples(sub)
  n <- attr(d, "Size")
  hc <- hclust(d, method = "ward.D2")
  labs <- cutree(hc, k = k)
  avg_sil <- NA_real_
  if (k >= 2L && k < n) {
    sil <- tryCatch(silhouette(labs, d)[, "sil_width"], error = function(e) NULL)
    if (!is.null(sil)) avg_sil <- mean(sil)
  }
  list(hc = hc, dist = d, labels = labs, avg_silhouette = avg_sil, mat = sub)
}

abnormal_enrichment <- function(labs, class3) {
  tab <- table(cluster = labs, Class3 = class3)
  abn_col <- which(colnames(tab) == "abnormal_case")
  if (length(abn_col) == 0) return(NA_real_)
  abn_per_cluster <- tab[, abn_col]
  max(abn_per_cluster) / sum(abn_per_cluster)
}

# ---------------------------------------------------------------------------
# Part 1: Global k-scan (top10000, all Class3)
# ---------------------------------------------------------------------------
message("=== Part 1: Global k-scan (top", PRIMARY_TOP_N, ") ===")
global_cpgs <- pick_top_var(val_mat, var_dt$CpG_ID, PRIMARY_TOP_N, MIN_SAMPLE_COVERAGE)
global_base <- cluster_at_k(val_mat, global_cpgs, k = 3L)

kscan_global <- rbindlist(lapply(K_RANGE, function(k) {
  r <- cluster_at_k(val_mat, global_cpgs, k)
  abn_in_best <- abnormal_enrichment(r$labels, meta$Class3)
  data.table(
    k = k,
    avg_silhouette = r$avg_silhouette,
    FM_vs_Class3 = fowlkes_mallows(r$labels, meta$Class3),
    abnormal_in_top_cluster_frac = abn_in_best,
    n_clusters = length(unique(r$labels))
  )
}))
fwrite(kscan_global, file.path(SUB, "kscan_global_silhouette.tsv"), sep = "\t")
fwrite(
  kscan_global[, .(k, FM_vs_Class3, avg_silhouette)],
  file.path(SUB, "fm_external_global_vs_class3.tsv"),
  sep = "\t"
)

# Cross-tabs per k
labs_k <- lapply(K_RANGE, function(k) cluster_at_k(val_mat, global_cpgs, k)$labels)
for (i in seq_along(K_RANGE)) {
  k <- K_RANGE[i]
  labs <- labs_k[[i]]
  xt <- as.data.table(as.data.frame.matrix(table(labs, meta$Class3)), keep.rownames = "cluster")
  fwrite(xt, file.path(SUB, sprintf("crosstab_global_k%d.tsv", k)), sep = "\t")
}

best_k_global <- kscan_global[which.max(avg_silhouette), k]
message("Best k by silhouette (global): ", best_k_global)

pdf(file.path(SUB, "kscan_global_silhouette.pdf"), width = 7, height = 5)
plot(kscan_global$k, kscan_global$avg_silhouette, type = "b", pch = 19,
     xlab = "Number of clusters (k)", ylab = "Average silhouette width",
     main = sprintf("Global top%d sample clustering", PRIMARY_TOP_N))
abline(v = 3, lty = 2, col = "grey50")
legend("topright", legend = c("max silhouette", "k=3 (clinical)"),
       lty = c(1, 2), pch = c(1, NA), col = c("black", "grey50"), bty = "n")
dev.off()

# ---------------------------------------------------------------------------
# Part 2: normal_case sub-clustering (within-cohort variance)
# ---------------------------------------------------------------------------
message("=== Part 2: normal_case sub-clustering ===")
normal_samp <- meta[Class3 == "normal_case", sample]
mat_normal <- val_mat[, normal_samp, drop = FALSE]
meta_normal <- meta[Class3 == "normal_case"]

message(sprintf("normal_case samples: %d", length(normal_samp)))
message("Computing variance within normal_case only...")
var_normal <- row_var(mat_normal)
var_normal_dt <- data.table(CpG_ID = rownames(mat_normal), variance = var_normal)
var_normal_dt <- var_normal_dt[!is.na(variance)]
setorder(var_normal_dt, -variance)

min_n_normal <- ceiling(ncol(mat_normal) * MIN_SAMPLE_COVERAGE)
cov_ok <- rowSums(!is.na(mat_normal)) >= min_n_normal
var_ids_normal <- var_normal_dt[CpG_ID %in% names(cov_ok)[cov_ok], CpG_ID]

for (top_n in c(5000L, 10000L)) {
  cpgs_n <- head(var_ids_normal, min(top_n, length(var_ids_normal)))
  message(sprintf("  normal-only top%d: %d CpGs", top_n, length(cpgs_n)))

  kscan_n <- rbindlist(lapply(K_RANGE, function(k) {
    r <- cluster_at_k(mat_normal, cpgs_n, k)
    data.table(
      k = k,
      avg_silhouette = r$avg_silhouette,
      FM_vs_Group2 = fowlkes_mallows(r$labels, meta_normal$Group2)
    )
  }))
  fwrite(kscan_n, file.path(SUB, sprintf("kscan_normal_top%d_silhouette.tsv", top_n)), sep = "\t")
  fwrite(
    kscan_n[, .(top_n = top_n, k, FM_vs_Group2, avg_silhouette)],
    file.path(SUB, sprintf("fm_external_normal_top%d_vs_group2.tsv", top_n)),
    sep = "\t"
  )

  best_k <- kscan_n[which.max(avg_silhouette), k]
  cl <- cluster_at_k(mat_normal, cpgs_n, best_k)

  # Group2 cross-tab (SPL / RPL / control mislabels shouldn't appear)
  cl_df <- data.table(
    sample = colnames(cl$mat),
    cluster = cl$labels,
    Group2 = meta_normal$Group2,
    Class3 = meta_normal$Class3
  )
  fwrite(cl_df, file.path(SUB, sprintf("normal_top%d_cluster_assignments_k%d.tsv", top_n, best_k)), sep = "\t")

  xt_g2 <- as.data.table(as.data.frame.matrix(table(cl$labels, meta_normal$Group2)), keep.rownames = "cluster")
  fwrite(xt_g2, file.path(SUB, sprintf("crosstab_normal_top%d_k%d_vs_Group2.tsv", top_n, best_k)), sep = "\t")

  # Heatmap for best k
  ht_opt$message <- FALSE
  g2_cols <- c(SPL = "#984EA3", RPL = "#FF7F00", control = "#999999")
  g2_use <- intersect(names(g2_cols), unique(meta_normal$Group2))
  ha <- HeatmapAnnotation(
    Group2 = meta_normal$Group2[match(colnames(cl$mat), meta_normal$sample)],
    col = list(Group2 = g2_cols[g2_use])
  )
  mat_hm <- cl$mat
  mat_hm[mat_hm < 0] <- 0
  mat_hm[mat_hm > 1] <- 1
  vord <- var_normal_dt[match(rownames(mat_hm), CpG_ID), variance]
  pdf(file.path(SUB, sprintf("heatmap_normal_top%d_k%d.pdf", top_n, best_k)), width = 12, height = 9)
  draw(Heatmap(
    mat_hm,
    name = "beta",
    col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE,
    row_order = order(vord, decreasing = TRUE),
    cluster_columns = cl$hc,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha,
    column_title = sprintf(
      "normal_case only (n=%d) | top%d CpGs (normal variance) | k=%d (max silhouette)",
      ncol(mat_hm), top_n, best_k
    ),
    row_title = "CpG (variance order)"
  ))
  dev.off()

  pdf(file.path(SUB, sprintf("kscan_normal_top%d_silhouette.pdf", top_n)), width = 7, height = 5)
  plot(kscan_n$k, kscan_n$avg_silhouette, type = "b", pch = 19,
       xlab = "k", ylab = "Average silhouette",
       main = sprintf("normal_case sub-clustering (top%d)", top_n))
  abline(v = best_k, lty = 2, col = "red")
  dev.off()
}

writeLines(c(
  "k-scan and normal_case sub-clustering",
  paste0("- Global: top ", PRIMARY_TOP_N, " CpGs (cohort variance), k=", paste(K_RANGE, collapse = ",")),
  paste0("- normal_case: recompute variance in n=", length(normal_samp), " samples only"),
  "- No imputation; Ward.D2 + pairwise Pearson distance",
  paste0("- Best global k by silhouette: ", best_k_global),
  "- Fowlkes-Mallows (external): cluster labels vs Class3 (global) or Group2 (normal subcluster)",
  "- Silhouette: internal cluster quality (not FM)"
), file.path(SUB, "README.txt"))

message("=== Done ===")
message("Output: ", SUB)
