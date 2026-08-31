#!/usr/bin/env Rscript
# High-variance CpGs computed ONLY in normal_case + control (abnormal excluded)
# Then cluster / evaluate separation between normal_case vs control

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "normal_control_hvar")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_N <- c(5000L, 10000L)
MIN_SAMPLE_COVERAGE <- 0.90
PRIMARY_TOP_N <- 10000L

message("=== Load ===")
clin <- fread("D:/ONT/clinical_649.tsv")
clin[, Sample_ID := as.character(Sample_ID)]
clin[, Class3 := fifelse(tolower(Group1) == "control", "control",
                 fifelse(tolower(Group4) == "abnormal", "abnormal_case", "normal_case"))]

cache <- readRDS(file.path(BASE, "clustering_output", "val_mat_samples.rds"))
val_mat <- cache$mat
sample_class3 <- cache$sample_class3

sample_na_rate <- colMeans(is.na(val_mat))
bad <- names(sample_na_rate)[sample_na_rate > 0.20]
if (length(bad)) {
  val_mat <- val_mat[, !(colnames(val_mat) %in% bad), drop = FALSE]
  sample_class3 <- sample_class3[colnames(val_mat)]
}

# Keep normal_case + control only
keep <- sample_class3 %in% c("normal_case", "control")
mat <- val_mat[, keep, drop = FALSE]
cls <- factor(sample_class3[keep], levels = c("control", "normal_case"))

cat("Samples (abnormal excluded):\n")
print(table(cls))

message("=== Variance within normal_case + control only ===")
row_var_nc <- apply(mat, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::var(x)
})
var_nc <- data.table(CpG_ID = rownames(mat), variance = row_var_nc)
var_nc <- var_nc[!is.na(variance)]
setorder(var_nc, -variance)
fwrite(var_nc, file.path(OUT, "CpG_variance_normal_control.tsv"), sep = "\t")

min_non_na <- ceiling(ncol(mat) * MIN_SAMPLE_COVERAGE)
cov_ok <- rowSums(!is.na(mat)) >= min_non_na
var_ids_ok <- var_nc[CpG_ID %in% names(cov_ok)[cov_ok], CpG_ID]
cat(sprintf("CpGs with >=%.0f%% coverage in NC subset: %d\n",
            MIN_SAMPLE_COVERAGE * 100, length(var_ids_ok)))

pearson_dist_samples <- function(m) {
  cm <- cor(m, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

pick_top <- function(n) head(var_ids_ok, min(n, length(var_ids_ok)))

cluster_nc <- function(cpgs) {
  sub <- mat[cpgs, , drop = FALSE]
  d <- pearson_dist_samples(sub)
  hc <- hclust(d, method = "ward.D2")
  list(mat = sub, hc = hc, dist = d)
}

eval_k2 <- function(labs, truth) {
  tab <- table(cluster = labs, truth = truth)
  # best cluster-truth alignment purity
  pur <- apply(tab, 1, max) / rowSums(tab)
  data.table(
    cluster = rownames(tab),
    n = rowSums(tab),
    dominant = colnames(tab)[apply(tab, 1, which.max)],
    purity = pur
  )
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

results <- list()
for (n in TOP_N) {
  cpgs <- pick_top(n)
  cl <- cluster_nc(cpgs)
  labs2 <- cutree(cl$hc, k = 2L)
  sil <- mean(silhouette(labs2, cl$dist)[, "sil_width"])

  pur <- eval_k2(labs2, cls)
  fwrite(pur, file.path(OUT, sprintf("cluster_purity_top%d_k2.tsv", n)), sep = "\t")

  # cross-tab
  xt <- as.data.table(as.data.frame.matrix(table(labs2, cls)), keep.rownames = "cluster")
  fwrite(xt, file.path(OUT, sprintf("crosstab_top%d_k2.tsv", n)), sep = "\t")

  # assignments
  assign <- data.table(sample = colnames(cl$mat), Class3 = cls, cluster_k2 = labs2)
  fwrite(assign, file.path(OUT, sprintf("cluster_assignments_top%d_k2.tsv", n)), sep = "\t")

  # simple accuracy if cluster labels mapped to dominant class
  mapped <- pur[, .(cluster = as.integer(cluster), pred = dominant)]
  pred <- mapped[match(labs2, mapped$cluster), pred]
  acc <- mean(pred == cls)

  results[[paste0("top", n)]] <- c(cl, list(
    top_n = n, cpgs = cpgs, labels_k2 = labs2,
    avg_silhouette_k2 = sil, accuracy_k2 = acc,
    FM_vs_Class3 = fowlkes_mallows(labs2, cls)
  ))
  message(sprintf("top%d: k=2 silhouette=%.3f, mapped accuracy=%.3f, FM_vs_Class3=%.3f",
                  n, sil, acc, results[[paste0("top", n)]]$FM_vs_Class3))
}

summary_dt <- rbindlist(lapply(results, function(r) {
  data.table(top_n = r$top_n, n_cpgs = length(r$cpgs),
             silhouette_k2 = r$avg_silhouette_k2,
             FM_vs_Class3 = r$FM_vs_Class3,
             accuracy_k2 = r$accuracy_k2)
}))
fwrite(summary_dt, file.path(OUT, "separation_summary.tsv"), sep = "\t")
fwrite(
  summary_dt[, .(top_n, k = 2L, FM_vs_Class3, silhouette_k2, accuracy_k2)],
  file.path(OUT, "fm_external_vs_class3.tsv"),
  sep = "\t"
)
print(summary_dt)

# Primary heatmap (top10000)
primary <- results[[paste0("top", PRIMARY_TOP_N)]]
ht_opt$message <- FALSE
class_cols <- c(control = "#4DAF4A", normal_case = "#377EB8")
ha <- HeatmapAnnotation(
  Class3 = cls[match(colnames(primary$mat), names(cls))],
  col = list(Class3 = class_cols)
)
mat_hm <- primary$mat
mat_hm[mat_hm < 0] <- 0
mat_hm[mat_hm > 1] <- 1
vord <- var_nc[match(rownames(mat_hm), CpG_ID), variance]

pdf(file.path(OUT, "heatmap_top10000_normal_vs_control.pdf"), width = 12, height = 9)
draw(Heatmap(
  mat_hm,
  name = "beta",
  col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = FALSE,
  row_order = order(vord, decreasing = TRUE),
  cluster_columns = primary$hc,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  column_title = sprintf(
    "control vs normal_case (n=%d) | top%d NC-variance CpGs | abnormal excluded",
    ncol(mat_hm), PRIMARY_TOP_N
  ),
  row_title = "CpG (NC variance order)"
))
dev.off()

# PCA for visual check (complete samples x top cpgs with pairwise... use zero only for PCA optional)
# Use samples with low NA on selected cpgs
message("=== PCA (top10000 NC-var CpGs) ===")
sub <- primary$mat
pc_mat <- t(sub)  # samples x cpgs
pc_mat[is.na(pc_mat)] <- 0  # for PCA visualization only
pc <- prcomp(pc_mat, center = TRUE, scale. = FALSE)
pc_df <- data.table(
  sample = rownames(pc_mat),
  Class3 = cls,
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2]
)
fwrite(pc_df, file.path(OUT, "pca_coordinates.tsv"), sep = "\t")

pdf(file.path(OUT, "pca_top10000_normal_vs_control.pdf"), width = 7, height = 6)
cols <- class_cols[as.character(pc_df$Class3)]
plot(pc_df$PC1, pc_df$PC2, col = cols, pch = 19, cex = 0.7,
     xlab = "PC1", ylab = "PC2",
     main = "normal_case vs control (NC high-var CpGs, abnormal excluded)")
legend("topright", legend = names(class_cols), col = class_cols, pch = 19, bty = "n")
dev.off()

writeLines(c(
  "normal_case + control only (abnormal_case excluded)",
  paste0("- Variance recomputed in n=", ncol(mat), " samples"),
  paste0("- Coverage filter: >=", MIN_SAMPLE_COVERAGE * 100, "% non-NA per CpG"),
  "- No imputation for clustering; PCA uses zero-fill for missing values only",
  paste0("- Primary top-N: ", PRIMARY_TOP_N),
  "- Fowlkes-Mallows (external): cluster k=2 vs Class3 (normal_case vs control)",
  "- Silhouette: internal cluster quality"
), file.path(OUT, "README.txt"))

message("Done: ", OUT)
