#!/usr/bin/env Rscript
# 步骤 1（已完成）：全体样本高变 CpG 聚类。
# 结论：abnormal 富集枝可分，并带上 30 例 later 标为 abnormal-like 的 normal_case。
# Approach A: control + normal_case + abnormal_case together
# No imputation: pairwise complete Pearson correlation for distances

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_N_TARGETS <- c(10000L, 20000L, 50000L)
PRIMARY_TOP_N <- 10000L
CLUSTER_K <- 3L
MIN_SAMPLE_COVERAGE <- 0.90

message("=== Step 1: Clinical Class3 ===")
clin <- fread("D:/ONT/clinical_649.tsv")
clin[, Sample_ID := as.character(Sample_ID)]
clin[, Class3 := fifelse(tolower(Group1) == "control", "control",
                 fifelse(tolower(Group4) == "abnormal", "abnormal_case", "normal_case"))]
fwrite(clin[, .(Sample_ID, Group1, Group2, Group4, Class3)],
       file.path(OUT, "sample_class3.tsv"), sep = "\t")

message("=== Step 2: Load methylation matrix ===")
mat_file <- file.path(BASE, "CpG_top100k_matrix.tsv")
hdr <- strsplit(readLines(mat_file, n = 1), "\t", fixed = TRUE)[[1]]
mat_samples <- hdr[-1]

clin_map <- clin[, .(key = tolower(Sample_ID), Class3)]
setkey(clin_map, key)
map_one <- function(s) {
  hit <- clin_map[list(tolower(s))]
  if (nrow(hit)) as.character(hit$Class3[1]) else NA_character_
}
sample_class3 <- vapply(mat_samples, map_one, character(1))
names(sample_class3) <- mat_samples

matched <- !is.na(sample_class3)
cat(sprintf("Matrix samples: %d | Matched to clinical: %d\n",
            length(mat_samples), sum(matched)))

keep_samples <- mat_samples[matched]
val_samples  <- keep_samples

mat_rds <- file.path(OUT, "val_mat_samples.rds")
if (file.exists(mat_rds)) {
  message("Loading cached sample matrix...")
  cache <- readRDS(mat_rds)
  val_mat <- cache$mat
  sample_class3 <- cache$sample_class3
} else {
  message("Reading matrix body...")
  lines <- readLines(mat_file)[-1L]
  n_samp <- length(mat_samples)
  val_mat <- matrix(NA_real_, nrow = length(lines), ncol = n_samp)
  cpg_ids <- character(length(lines))
  for (i in seq_along(lines)) {
    line <- lines[i]
    sp <- strsplit(line, " ", fixed = TRUE)[[1]]
    cpg_ids[i] <- sp[1L]
    x <- sp[-1L]
    x[x == "NA"] <- NA
    val_mat[i, ] <- suppressWarnings(as.numeric(x))
  }
  rm(lines)
  gc()
  colnames(val_mat) <- mat_samples
  rownames(val_mat) <- cpg_ids
  val_mat <- val_mat[, keep_samples, drop = FALSE]
  sample_class3 <- sample_class3[keep_samples]
  saveRDS(list(mat = val_mat, sample_class3 = sample_class3), mat_rds)
}

sample_na_rate <- colMeans(is.na(val_mat))
fwrite(data.table(sample = names(sample_na_rate), na_rate = sample_na_rate,
                  Class3 = sample_class3),
       file.path(OUT, "sample_na_rate.tsv"), sep = "\t")

bad_samples <- names(sample_na_rate)[sample_na_rate > 0.20]
if (length(bad_samples)) {
  message(sprintf("Dropping %d samples with >20%% NA.", length(bad_samples)))
  val_mat <- val_mat[, !(colnames(val_mat) %in% bad_samples), drop = FALSE]
  sample_class3 <- sample_class3[colnames(val_mat)]
}

cat("Class3 counts:\n")
print(table(sample_class3))

message("=== Step 3: Variance order + coverage filter ===")
var_dt <- fread(file.path(BASE, "CpG_variance.tsv"), header = FALSE,
                col.names = c("CpG_ID", "count", "variance"))
setorder(var_dt, -variance)
var_dt <- var_dt[CpG_ID %in% rownames(val_mat)]

min_non_na <- ceiling(ncol(val_mat) * MIN_SAMPLE_COVERAGE)
non_na_count <- rowSums(!is.na(val_mat))
coverage_ok <- non_na_count >= min_non_na
var_ids_ok <- var_dt[CpG_ID %in% names(coverage_ok)[coverage_ok], CpG_ID]
cat(sprintf("CpGs with >=%.0f%% sample coverage: %d / %d\n",
            MIN_SAMPLE_COVERAGE * 100, length(var_ids_ok), nrow(val_mat)))

pick_top_n <- function(n) {
  head(var_ids_ok, min(n, length(var_ids_ok)))
}

pearson_dist_samples <- function(mat) {
  cm <- cor(mat, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

pearson_dist_cpgs <- function(mat) {
  cm <- cor(t(mat), use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

fowlkes_mallows <- function(lab1, lab2) {
  # External cluster validation: compare two partitions of the same samples.
  # Use with (cluster_labels, known_truth_labels), not two unsupervised runs.
  lab1 <- as.integer(factor(lab1))
  lab2 <- as.integer(factor(lab2))
  tab <- table(lab1, lab2)
  tk <- sum(choose(tab, 2))
  pk <- sum(choose(rowSums(tab), 2))
  qk <- sum(choose(colSums(tab), 2))
  if (pk == 0 || qk == 0) return(NA_real_)
  sqrt(tk / pk * tk / qk)
}

if (.Platform$OS.type == "windows") {
  tryCatch(memory.limit(size = 32000), error = function(e) NULL)
}

cluster_samples <- function(cpgs, mat) {
  sub <- mat[cpgs, , drop = FALSE]
  hc <- hclust(pearson_dist_samples(sub), method = "ward.D2")
  list(mat = sub, hc_sample = hc, cluster_k3 = cutree(hc, k = CLUSTER_K))
}

message("=== Step 4: Multi-gradient sample clustering ===")
results <- list()
for (n in TOP_N_TARGETS) {
  cpgs <- pick_top_n(n)
  message(sprintf("  top%d: using %d CpGs...", n, length(cpgs)))
  results[[paste0("top", n)]] <- c(list(top_n = n, cpgs = cpgs), cluster_samples(cpgs, val_mat))
}

usage <- rbindlist(lapply(results, function(r) {
  data.table(target_topN = r$top_n, cpgs_used = length(r$cpgs),
             min_non_na_required = min_non_na)
}))
fwrite(usage, file.path(OUT, "complete_cpg_counts.tsv"), sep = "\t")
print(usage)

message("=== Step 5: External FM — cluster k vs Class3 (known labels) ===")
fm_external <- rbindlist(lapply(results, function(r) {
  labs <- r$cluster_k3
  truth <- sample_class3[colnames(r$mat)]
  data.table(
    top_n = r$top_n,
    k = CLUSTER_K,
    n_cpgs = length(r$cpgs),
    FM_vs_Class3 = fowlkes_mallows(labs, truth)
  )
}))
fwrite(fm_external, file.path(OUT, "fm_external_vs_class3.tsv"), sep = "\t")
print(fm_external)

message(sprintf("=== Step 6: Cluster vs Class3 (top%d) ===", PRIMARY_TOP_N))
primary <- results[[paste0("top", PRIMARY_TOP_N)]]

cluster_df <- data.table(
  sample = colnames(primary$mat),
  Class3 = sample_class3,
  cluster_k3 = primary$cluster_k3
)
fwrite(cluster_df, file.path(OUT, "cluster_vs_class3.tsv"), sep = "\t")

xt <- as.data.frame.matrix(table(primary$cluster_k3, cluster_df$Class3))
xt$cluster <- rownames(xt)
fwrite(as.data.table(xt), file.path(OUT, "crosstab_cluster_class3.tsv"), sep = "\t")
cat("\nCross-tab (cluster k=3 vs Class3):\n")
print(table(primary$cluster_k3, cluster_df$Class3))

pur <- cluster_df[, .(
  n = .N,
  dominant = names(sort(table(Class3), decreasing = TRUE))[1],
  purity = max(table(Class3)) / .N
), by = cluster_k3]
fwrite(pur, file.path(OUT, "cluster_purity.tsv"), sep = "\t")
print(pur)

message(sprintf("=== Step 7: Heatmap (top%d) ===", PRIMARY_TOP_N))
ht_opt$message <- FALSE
ht_opt$raster_temp_image_max_width <- 12000
ht_opt$raster_temp_image_max_height <- 12000
class_cols <- c(control = "#4DAF4A", normal_case = "#377EB8", abnormal_case = "#E41A1C")
ha <- HeatmapAnnotation(
  Class3 = cluster_df$Class3[match(colnames(primary$mat), cluster_df$sample)],
  col = list(Class3 = class_cols)
)
mat_hm <- primary$mat
mat_hm[mat_hm < 0] <- 0
mat_hm[mat_hm > 1] <- 1
var_lookup <- var_dt[, .(CpG_ID, variance)]
vord <- var_lookup[match(rownames(mat_hm), CpG_ID), variance]
row_order <- order(vord, decreasing = TRUE)

hm_args <- list(
  mat_hm,
  name = "beta",
  col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = FALSE,
  row_order = row_order,
  cluster_columns = primary$hc_sample,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  column_title = sprintf(
    "n=%d samples | top%d CpGs (>=%.0f%% coverage, pairwise cor, no imputation)",
    ncol(mat_hm), PRIMARY_TOP_N, MIN_SAMPLE_COVERAGE * 100
  ),
  row_title = "CpG sites (variance order)"
)

pdf(file.path(OUT, "heatmap_top10000_class3.pdf"), width = 14, height = 10)
do.call(Heatmap, hm_args)
dev.off()

message("Heatmap saved: heatmap_top10000_class3.pdf")

saveRDS(results, file.path(OUT, "all_cluster_results.rds"))

writeLines(c(
  "Method summary",
  paste0("- Primary analysis: top ", PRIMARY_TOP_N, " variable CpGs"),
  paste0("- Samples: all Class3 groups (Approach A), n=", ncol(val_mat)),
  paste0("- CpG filter: variance-ranked, >=", MIN_SAMPLE_COVERAGE * 100, "% non-NA per CpG"),
  "- Missing: no imputation; Pearson cor uses pairwise.complete.obs",
  "- Clustering: Ward.D2 on 1 - Pearson(sample correlation)",
  paste0("- Fowlkes-Mallows (external): cluster k=", CLUSTER_K, " vs known Class3 labels"),
  "- Silhouette / cross-tab / purity: see cluster_purity.tsv"
), file.path(OUT, "method_summary.txt"))

message("=== Done ===")
message("Output: ", OUT)
