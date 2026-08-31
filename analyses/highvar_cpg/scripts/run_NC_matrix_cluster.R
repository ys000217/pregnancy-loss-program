#!/usr/bin/env Rscript
# 步骤 2 敏感性测试（已完成）：去掉 abnormal 后用 NC 专用高变矩阵聚类。
# 30 例 normal_case 小簇仍在 → abnormal-like，不是临床 abnormal。
# Input: CpG_matrix_NC.tsv + NC_matrix_sample_annotation.tsv
# No imputation; pairwise Pearson distance; Ward.D2

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "NC_matrix_analysis")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MAT_FILE  <- file.path(BASE, "CpG_matrix_NC.tsv")
ANN_FILE  <- file.path(BASE, "NC_matrix_sample_annotation.tsv")
VAR_FILE  <- file.path(BASE, "CpG_variance_NC.tsv")
TOP_FILE  <- file.path(BASE, "CpG_topNC_10000.list")
K_RANGE   <- 2:8

message("=== Load annotation ===")
ann <- fread(ANN_FILE)
ann[, sample_id := as.character(sample_id)]
cls <- factor(ann$Class3, levels = c("control", "normal_case"))
names(cls) <- ann$sample_id
print(table(cls))

message("=== Load matrix (space-separated body) ===")
hdr <- strsplit(readLines(MAT_FILE, n = 1L), "\t", fixed = TRUE)[[1]]
mat_samples <- hdr[-1]
stopifnot(identical(mat_samples, ann$sample_id))

lines <- readLines(MAT_FILE)[-1L]
n_cpg <- length(lines)
n_samp <- length(mat_samples)
mat <- matrix(NA_real_, nrow = n_cpg, ncol = n_samp)
cpg_ids <- character(n_cpg)
for (i in seq_along(lines)) {
  sp <- strsplit(lines[i], " ", fixed = TRUE)[[1]]
  cpg_ids[i] <- sp[1L]
  x <- sp[-1L]
  if (length(x) != n_samp) {
    # fallback tab
    sp <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    cpg_ids[i] <- sp[1L]
    x <- sp[-1L]
  }
  x[x == "NA"] <- NA
  mat[i, ] <- suppressWarnings(as.numeric(x))
}
rm(lines)
gc()
rownames(mat) <- cpg_ids
colnames(mat) <- mat_samples
cat(sprintf("Matrix: %d CpGs x %d samples\n", nrow(mat), ncol(mat)))
cat(sprintf("Overall NA rate: %.4f\n", mean(is.na(mat))))

# variance order for row display (from NC variance file)
message("=== Load NC variance ranks for heatmap row order ===")
top_ids <- fread(TOP_FILE, header = FALSE)$V1
var_dt <- fread(VAR_FILE, header = FALSE, col.names = c("CpG_ID", "count", "variance"))
setorder(var_dt, -variance)
# keep only matrix CpGs
var_lookup <- var_dt[CpG_ID %in% cpg_ids]
vord <- match(rownames(mat), var_lookup$CpG_ID)
# if match fails somehow, use existing order
row_order_var <- order(var_lookup[match(rownames(mat), CpG_ID), variance],
                       decreasing = TRUE, na.last = TRUE)

pearson_dist_samples <- function(m) {
  cm <- cor(m, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

fowlkes_mallows <- function(lab1, lab2) {
  # External validation only: cluster partition vs known Class3 labels.
  lab1 <- as.integer(factor(lab1))
  lab2 <- as.integer(factor(lab2))
  tab <- table(lab1, lab2)
  tk <- sum(choose(tab, 2))
  pk <- sum(choose(rowSums(tab), 2))
  qk <- sum(choose(colSums(tab), 2))
  if (pk == 0 || qk == 0) return(NA_real_)
  sqrt(tk / pk * tk / qk)
}

message("=== Sample distance + hierarchical clustering ===")
d <- pearson_dist_samples(mat)
hc <- hclust(d, method = "ward.D2")
saveRDS(list(hc = hc, dist = d), file.path(OUT, "hclust_ward.rds"))

# k-scan
message("=== k-scan silhouette ===")
kscan <- rbindlist(lapply(K_RANGE, function(k) {
  labs <- cutree(hc, k = k)
  sil <- mean(silhouette(labs, d)[, "sil_width"])
  xt <- table(cluster = labs, Class3 = cls[colnames(mat)])
  data.table(
    k = k,
    avg_silhouette = sil,
    FM_vs_Class3 = fowlkes_mallows(labs, cls[colnames(mat)]),
    crosstab = paste(capture.output(print(xt)), collapse = " | ")
  )
}))
fwrite(kscan[, .(k, avg_silhouette, FM_vs_Class3)], file.path(OUT, "kscan_silhouette.tsv"), sep = "\t")
fwrite(
  kscan[, .(k, FM_vs_Class3, avg_silhouette)],
  file.path(OUT, "fm_external_vs_class3.tsv"),
  sep = "\t"
)
print(kscan[, .(k, avg_silhouette, FM_vs_Class3)])

best_k <- kscan[which.max(avg_silhouette), k]
message("Best k by silhouette: ", best_k)

pdf(file.path(OUT, "kscan_silhouette.pdf"), width = 7, height = 5)
plot(kscan$k, kscan$avg_silhouette, type = "b", pch = 19,
     xlab = "k", ylab = "Average silhouette",
     main = "NC high-var CpGs: sample clustering")
abline(v = 2, lty = 2, col = "grey50")
abline(v = best_k, lty = 2, col = "red")
legend("topright", legend = c("k=2 (binary)", paste0("best k=", best_k)),
       lty = 2, col = c("grey50", "red"), bty = "n")
dev.off()

# Primary: k=2 (binary control vs normal_case question)
message("=== Primary k=2 vs Class3 ===")
labs2 <- cutree(hc, k = 2L)
assign2 <- data.table(
  sample = colnames(mat),
  Class3 = as.character(cls[colnames(mat)]),
  cluster_k2 = labs2
)
fwrite(assign2, file.path(OUT, "cluster_assignments_k2.tsv"), sep = "\t")

xt2 <- as.data.table(as.data.frame.matrix(table(labs2, cls)), keep.rownames = "cluster")
fwrite(xt2, file.path(OUT, "crosstab_k2.tsv"), sep = "\t")
cat("\nCross-tab k=2:\n")
print(table(labs2, cls))

pur2 <- assign2[, .(
  n = .N,
  dominant = names(sort(table(Class3), decreasing = TRUE))[1],
  purity = max(table(Class3)) / .N,
  n_control = sum(Class3 == "control"),
  n_normal = sum(Class3 == "normal_case")
), by = cluster_k2]
fwrite(pur2, file.path(OUT, "cluster_purity_k2.tsv"), sep = "\t")
print(pur2)

# mapped accuracy (majority vote)
mapped <- pur2[, .(cluster_k2, pred = dominant)]
pred <- mapped[match(labs2, mapped$cluster_k2), pred]
acc <- mean(pred == as.character(cls))
baseline <- max(table(cls)) / length(cls)
cat(sprintf("Mapped accuracy k=2: %.3f (baseline majority=%.3f)\n", acc, baseline))

# Also report best_k crosstab
labs_b <- cutree(hc, k = best_k)
assign_b <- data.table(sample = colnames(mat), Class3 = as.character(cls),
                       cluster = labs_b)
fwrite(assign_b, file.path(OUT, sprintf("cluster_assignments_k%d.tsv", best_k)), sep = "\t")
xtb <- as.data.table(as.data.frame.matrix(table(labs_b, cls)), keep.rownames = "cluster")
fwrite(xtb, file.path(OUT, sprintf("crosstab_k%d.tsv", best_k)), sep = "\t")
cat(sprintf("\nCross-tab best k=%d:\n", best_k))
print(table(labs_b, cls))

# Heatmap k=2 annotation
message("=== Heatmap ===")
ht_opt$message <- FALSE
class_cols <- c(control = "#4DAF4A", normal_case = "#377EB8")
ha <- HeatmapAnnotation(
  Class3 = as.character(cls[colnames(mat)]),
  col = list(Class3 = class_cols)
)
mat_hm <- mat
mat_hm[mat_hm < 0] <- 0
mat_hm[mat_hm > 1] <- 1

pdf(file.path(OUT, "heatmap_NC_top10000_class3.pdf"), width = 12, height = 9)
draw(Heatmap(
  mat_hm,
  name = "beta",
  col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = FALSE,
  row_order = row_order_var,
  cluster_columns = hc,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  column_title = sprintf(
    "control vs normal_case (n=%d) | NC-specific top10000 var CpGs | no imputation",
    ncol(mat_hm)
  ),
  row_title = "CpG (NC variance order)"
))
dev.off()

# PCA (zero-fill only for visualization)
message("=== PCA ===")
X <- t(mat)
X[is.na(X)] <- 0
pc <- prcomp(X, center = TRUE, scale. = FALSE)
pc_df <- data.table(
  sample = rownames(X),
  Class3 = as.character(cls),
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2]
)
fwrite(pc_df, file.path(OUT, "pca_coordinates.tsv"), sep = "\t")
ve <- summary(pc)$importance[2, 1:2] * 100

pdf(file.path(OUT, "pca_NC_top10000.pdf"), width = 7, height = 6)
cols <- class_cols[pc_df$Class3]
plot(pc_df$PC1, pc_df$PC2, col = cols, pch = 19, cex = 0.65,
     xlab = sprintf("PC1 (%.1f%%)", ve[1]),
     ylab = sprintf("PC2 (%.1f%%)", ve[2]),
     main = "NC high-var CpGs PCA (abnormal excluded)")
legend("topright", legend = names(class_cols), col = class_cols, pch = 19, bty = "n")
dev.off()

# Summary
summary_dt <- data.table(
  n_samples = ncol(mat),
  n_control = sum(cls == "control"),
  n_normal_case = sum(cls == "normal_case"),
  n_cpgs = nrow(mat),
  overall_na_rate = mean(is.na(mat)),
  best_k_silhouette = best_k,
  max_silhouette = kscan[which.max(avg_silhouette), avg_silhouette],
  silhouette_k2 = kscan[k == 2, avg_silhouette],
  FM_vs_Class3_k2 = kscan[k == 2, FM_vs_Class3],
  mapped_accuracy_k2 = acc,
  baseline_majority = baseline
)
fwrite(summary_dt, file.path(OUT, "analysis_summary.tsv"), sep = "\t")

writeLines(c(
  "NC-specific high-var CpG clustering (abnormal excluded at variance stage)",
  paste0("- Matrix: CpG_matrix_NC.tsv (", nrow(mat), " x ", ncol(mat), ")"),
  "- Distance: 1 - pairwise Pearson; method: Ward.D2",
  "- No imputation",
  paste0("- Best k by silhouette: ", best_k, " (sil=", round(summary_dt$max_silhouette, 3), ")"),
  paste0("- Fowlkes-Mallows (external): cluster labels vs Class3; k=2 FM=", round(summary_dt$FM_vs_Class3_k2, 3)),
  "- Silhouette: internal cluster quality (not FM)",
  paste0("- k=2 mapped accuracy: ", round(acc, 3), " vs majority baseline ", round(baseline, 3))
), file.path(OUT, "README.txt"))

message("=== Done ===")
message("Output: ", OUT)
print(summary_dt)
