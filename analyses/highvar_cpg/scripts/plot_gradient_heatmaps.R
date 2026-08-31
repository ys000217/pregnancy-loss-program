#!/usr/bin/env Rscript
# Draw sample-clustering heatmaps at multiple CpG-count gradients.
# k is NOT used: Ward dendrogram is shown intact; only clinical Class3 is annotated.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "gradient_heatmaps")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MIN_SAMPLE_COVERAGE <- 0.90
FULL_TOP_N <- c(5000L, 10000L, 20000L, 50000L)
NC_TOP_N   <- c(2000L, 5000L, 10000L)

class_cols <- c(control = "#4DAF4A", normal_case = "#377EB8", abnormal_case = "#E41A1C")

pearson_dist_samples <- function(mat) {
  cm <- cor(mat, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

read_space_matrix <- function(path) {
  hdr <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]]
  mat_samples <- hdr[-1]
  lines <- readLines(path)[-1L]
  n_samp <- length(mat_samples)
  val_mat <- matrix(NA_real_, nrow = length(lines), ncol = n_samp)
  cpg_ids <- character(length(lines))
  for (i in seq_along(lines)) {
    sp <- strsplit(lines[i], " ", fixed = TRUE)[[1]]
    if (length(sp) != n_samp + 1L) {
      sp <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    }
    cpg_ids[i] <- sp[1L]
    x <- sp[-1L]
    x[x == "NA"] <- NA
    val_mat[i, ] <- suppressWarnings(as.numeric(x))
  }
  rm(lines)
  gc()
  rownames(val_mat) <- cpg_ids
  colnames(val_mat) <- mat_samples
  val_mat
}

draw_one <- function(mat, hc, class3, title, pdf_path, png_path, var_order_ids) {
  ht_opt$message <- FALSE
  ht_opt$raster_temp_image_max_width <- 12000
  ht_opt$raster_temp_image_max_height <- 12000
  mat_hm <- mat
  mat_hm[mat_hm < 0] <- 0
  mat_hm[mat_hm > 1] <- 1
  row_order <- match(var_order_ids, rownames(mat_hm))
  row_order <- row_order[!is.na(row_order)]
  if (length(row_order) != nrow(mat_hm)) {
    row_order <- seq_len(nrow(mat_hm))
  }
  keep_lv <- intersect(names(class_cols), unique(as.character(class3)))
  ha <- HeatmapAnnotation(
    Class3 = as.character(class3),
    col = list(Class3 = class_cols[keep_lv]),
    annotation_name_side = "left"
  )
  ht <- Heatmap(
    mat_hm,
    name = "beta",
    col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE,
    row_order = row_order,
    cluster_columns = hc,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha,
    column_title = title,
    row_title = "CpG (variance rank)",
    use_raster = TRUE,
    raster_quality = 5
  )
  pdf(pdf_path, width = 14, height = 10)
  draw(ht)
  dev.off()
  png(png_path, width = 2800, height = 2000, res = 180)
  draw(ht)
  dev.off()
  message("Wrote: ", pdf_path)
}

# ---------------------------------------------------------------------------
# Panel 1: all samples, cohort-wide high-var CpGs
# ---------------------------------------------------------------------------
message("=== Panel 1: full cohort ===")
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

saved <- file.path(BASE, "clustering_output", "all_cluster_results.rds")
old <- if (file.exists(saved)) readRDS(saved) else list()

for (n in FULL_TOP_N) {
  key <- paste0("top", n)
  cpgs <- head(var_ids_ok, min(n, length(var_ids_ok)))
  if (!is.null(old[[key]]$hc_sample) &&
      identical(old[[key]]$cpgs, cpgs) &&
      identical(colnames(old[[key]]$mat), colnames(val_mat))) {
    message(sprintf("Reuse cached hclust for full top%d", n))
    hc <- old[[key]]$hc_sample
    sub <- old[[key]]$mat
  } else {
    message(sprintf("Clustering full cohort, top%d (%d CpGs)...", n, length(cpgs)))
    sub <- val_mat[cpgs, , drop = FALSE]
    hc <- hclust(pearson_dist_samples(sub), method = "ward.D2")
  }
  draw_one(
    mat = sub,
    hc = hc,
    class3 = sample_class3[colnames(sub)],
    title = sprintf(
      "Full cohort (n=%d) | top %d cohort-var CpGs | Ward + Pearson | no k-cut",
      ncol(sub), length(cpgs)
    ),
    pdf_path = file.path(OUT, sprintf("full_top%d.pdf", n)),
    png_path = file.path(OUT, sprintf("full_top%d.png", n)),
    var_order_ids = cpgs
  )
  rm(sub, hc)
  gc()
}

# ---------------------------------------------------------------------------
# Panel 2: abnormal removed; NC-specific high-var CpGs (matrix is top 10000)
# ---------------------------------------------------------------------------
message("=== Panel 2: NC (abnormal excluded) ===")
nc_rds <- file.path(OUT, "nc_mat_cache.rds")
if (file.exists(nc_rds)) {
  nc <- readRDS(nc_rds)
  nc_mat <- nc$mat
  nc_cls <- nc$class3
} else {
  ann <- fread(file.path(BASE, "NC_matrix_sample_annotation.tsv"))
  ann[, sample_id := as.character(sample_id)]
  nc_mat <- read_space_matrix(file.path(BASE, "CpG_matrix_NC.tsv"))
  stopifnot(identical(colnames(nc_mat), ann$sample_id))
  nc_cls <- ann$Class3
  names(nc_cls) <- ann$sample_id
  saveRDS(list(mat = nc_mat, class3 = nc_cls), nc_rds)
}

var_nc <- fread(file.path(BASE, "CpG_variance_NC.tsv"), header = FALSE,
                col.names = c("CpG_ID", "count", "variance"))
setorder(var_nc, -variance)
var_nc_ids <- var_nc[CpG_ID %in% rownames(nc_mat), CpG_ID]
if (length(var_nc_ids) < nrow(nc_mat)) {
  extra <- setdiff(rownames(nc_mat), var_nc_ids)
  var_nc_ids <- c(var_nc_ids, extra)
}

for (n in NC_TOP_N) {
  cpgs <- head(var_nc_ids, min(n, length(var_nc_ids)))
  message(sprintf("Clustering NC, top%d (%d CpGs)...", n, length(cpgs)))
  sub <- nc_mat[cpgs, , drop = FALSE]
  hc <- hclust(pearson_dist_samples(sub), method = "ward.D2")
  draw_one(
    mat = sub,
    hc = hc,
    class3 = nc_cls[colnames(sub)],
    title = sprintf(
      "NC only, abnormal removed (n=%d) | top %d NC-var CpGs | Ward + Pearson | no k-cut",
      ncol(sub), length(cpgs)
    ),
    pdf_path = file.path(OUT, sprintf("NC_top%d.pdf", n)),
    png_path = file.path(OUT, sprintf("NC_top%d.png", n)),
    var_order_ids = cpgs
  )
  rm(sub, hc)
  gc()
}

writeLines(c(
  "Gradient heatmaps (no k-cut)",
  "",
  "k is not used. Sample order is the Ward dendrogram on 1-Pearson.",
  "Color bar is clinical Class3 only (not unsupervised cluster IDs).",
  "Rows: CpGs in decreasing variance rank; rows are not clustered.",
  "",
  "Panel 1 full_top*.pdf: all Class3 groups; variance from full cohort.",
  paste0("  gradients: ", paste(FULL_TOP_N, collapse = ", ")),
  "Panel 2 NC_top*.pdf: control + normal_case; variance from NC-only job 09.",
  paste0("  gradients: ", paste(NC_TOP_N, collapse = ", "),
         " (NC downloaded matrix has 10000 sites)"),
  "",
  "Missing values: no imputation; pairwise complete Pearson."
), file.path(OUT, "README.txt"))

message("=== Done ===")
message("Output: ", OUT)