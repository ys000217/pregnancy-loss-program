#!/usr/bin/env Rscript
# PPT figures Fig1–Fig3 + matching Ward heatmaps (Nature style).
# Fig4 removed per user request.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "ppt_figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MM <- 1 / 25.4
COL_ABN  <- "#D55E00"
COL_CTRL <- "#009E73"
COL_NORM <- "#0072B2"
COL_SUS  <- "#E69F00"

like30 <- fread("D:/ONT/analyses/highvar_cpg/metadata/abnormal_like_normal_case_30.txt",
                header = FALSE)$V1
like_key <- unique(tolower(like30))

nature_theme <- function(base_size = 7) {
  theme_bw(base_size = base_size, base_family = "sans") %+replace%
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.2, "mm"),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(3, "mm"),
      plot.tag = element_text(face = "bold", size = 8),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
}

save_gg <- function(stem, plot, w_mm, h_mm) {
  pdf_path <- file.path(OUT, paste0(stem, ".pdf"))
  png_path <- file.path(OUT, paste0(stem, ".png"))
  tryCatch(
    ggsave(pdf_path, plot, width = w_mm * MM, height = h_mm * MM,
           device = cairo_pdf, useDingbats = FALSE),
    error = function(e) ggsave(pdf_path, plot, width = w_mm * MM, height = h_mm * MM,
                               useDingbats = FALSE)
  )
  ggsave(png_path, plot, width = w_mm * MM, height = h_mm * MM, dpi = 600)
  message("Wrote ", stem)
}

read_space_matrix <- function(path) {
  hdr <- strsplit(readLines(path, n = 1L), "\t", fixed = TRUE)[[1]]
  mat_samples <- hdr[-1]
  lines <- readLines(path)[-1L]
  n_samp <- length(mat_samples)
  val <- matrix(NA_real_, nrow = length(lines), ncol = n_samp)
  ids <- character(length(lines))
  for (i in seq_along(lines)) {
    sp <- strsplit(lines[i], " ", fixed = TRUE)[[1]]
    if (length(sp) != n_samp + 1L) sp <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    ids[i] <- sp[1L]
    x <- sp[-1L]
    x[x == "NA"] <- NA
    val[i, ] <- suppressWarnings(as.numeric(x))
  }
  rownames(val) <- ids
  colnames(val) <- mat_samples
  val
}

pearson_dist <- function(mat) {
  cm <- cor(mat, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

assign_group <- function(samples, class3) {
  g <- as.character(class3[samples])
  sus <- tolower(samples) %in% like_key & g == "normal_case"
  g[sus] <- "suspected_abnormal"
  factor(g, levels = c("control", "normal_case", "suspected_abnormal", "abnormal_case"))
}

group_cols <- c(
  control = COL_CTRL,
  normal_case = COL_NORM,
  suspected_abnormal = COL_SUS,
  abnormal_case = COL_ABN
)

save_heatmap <- function(stem, mat, groups, title, pdf_w = 7, pdf_h = 5) {
  ht_opt$message <- FALSE
  ht_opt$raster_temp_image_max_width <- 12000
  ht_opt$raster_temp_image_max_height <- 12000
  mat_hm <- mat
  mat_hm[mat_hm < 0] <- 0
  mat_hm[mat_hm > 1] <- 1
  hc <- hclust(pearson_dist(mat_hm), method = "ward.D2")
  lv <- intersect(names(group_cols), levels(groups))
  ha <- HeatmapAnnotation(
    Group = as.character(groups),
    col = list(Group = group_cols[lv]),
    annotation_name_side = "left",
    simple_anno_size = unit(3, "mm")
  )
  ht <- Heatmap(
    mat_hm,
    name = "beta",
    col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE,
    cluster_columns = hc,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha,
    column_title = title,
    row_title = "CpG (variance order)",
    use_raster = TRUE,
    raster_quality = 5
  )
  pdf(file.path(OUT, paste0(stem, ".pdf")), width = pdf_w, height = pdf_h)
  draw(ht)
  dev.off()
  png(file.path(OUT, paste0(stem, ".png")), width = pdf_w * 300, height = pdf_h * 300, res = 300)
  draw(ht)
  dev.off()
  message("Wrote ", stem, "_heatmap")
}

pca_df_from_mat <- function(mat, class3) {
  X <- t(mat)
  na <- is.na(X)
  rs <- rowMeans(X, na.rm = TRUE)
  X[na] <- rs[row(X)[na]]
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  data.table(
    sample = rownames(pc$x),
    Class3 = as.character(class3[rownames(pc$x)]),
    PC1 = pc$x[, 1], PC2 = pc$x[, 2],
    xlab = sprintf("PC1 (%.1f%%)", ve[1]),
    ylab = sprintf("PC2 (%.1f%%)", ve[2])
  )
}

# ---------- Fig1 ----------
message("Fig1")
cache <- readRDS(file.path(BASE, "clustering_output", "val_mat_samples.rds"))
val_mat <- cache$mat
sample_class3 <- cache$sample_class3
na_rate <- colMeans(is.na(val_mat))
bad <- names(na_rate)[na_rate > 0.20]
if (length(bad)) {
  val_mat <- val_mat[, !(colnames(val_mat) %in% bad), drop = FALSE]
  sample_class3 <- sample_class3[colnames(val_mat)]
}
var_dt <- fread(file.path(BASE, "CpG_variance.tsv"), header = FALSE,
                col.names = c("CpG_ID", "count", "variance"))
setorder(var_dt, -variance)
min_non_na <- ceiling(ncol(val_mat) * 0.90)
ok <- rowSums(!is.na(val_mat)) >= min_non_na
cpgs <- head(var_dt[CpG_ID %in% rownames(val_mat[ok, , drop = FALSE]), CpG_ID], 10000L)
sub1 <- val_mat[cpgs, , drop = FALSE]

df1 <- pca_df_from_mat(sub1, sample_class3)
df1[, group := assign_group(sample, sample_class3)]

p1 <- ggplot(df1, aes(PC1, PC2, colour = group, shape = group)) +
  geom_point(data = df1[group %in% c("control", "normal_case")],
             size = 1.0, alpha = 0.65, stroke = 0.15) +
  geom_point(data = df1[group %in% c("suspected_abnormal", "abnormal_case")],
             size = 1.4, alpha = 0.95, stroke = 0.2) +
  scale_colour_manual(values = group_cols,
    labels = c(control = "Control", normal_case = "Normal case",
               suspected_abnormal = "Suspected abnormal (n=30)",
               abnormal_case = "abnormal"), name = NULL) +
  scale_shape_manual(values = c(control = 16, normal_case = 17,
               suspected_abnormal = 15, abnormal_case = 15),
    labels = c(control = "Control", normal_case = "Normal case",
               suspected_abnormal = "Suspected abnormal (n=30)",
               abnormal_case = "abnormal"), name = NULL) +
  labs(x = df1$xlab[1], y = df1$ylab[1]) +
  nature_theme(7) + theme(legend.position = "bottom", legend.box = "vertical")
save_gg("Fig1", p1, 120, 100)
fwrite(df1[, .(sample, Class3, group, PC1, PC2)], file.path(OUT, "Fig1_source.csv"))

g1 <- assign_group(colnames(sub1), sample_class3)
save_heatmap("Fig1", sub1, g1,
  sprintf("Full cohort n=%d | top10k cohort-var CpGs | Ward, no k-cut", ncol(sub1)))

rm(val_mat, cache)
gc()

# ---------- Fig2 ----------
message("Fig2")
nc <- readRDS(file.path(BASE, "clustering_output", "gradient_heatmaps", "nc_mat_cache.rds"))
nc_mat <- nc$mat
nc_cls <- nc$class3
df2a <- pca_df_from_mat(nc_mat, nc_cls)
df2a[, group := assign_group(sample, nc_cls)]
df2a <- df2a[group != "abnormal_case"]

p2a <- ggplot(df2a, aes(PC1, PC2, colour = group, shape = group)) +
  geom_point(data = df2a[group != "suspected_abnormal"], size = 1.0, alpha = 0.65) +
  geom_point(data = df2a[group == "suspected_abnormal"], size = 1.5, alpha = 0.95) +
  scale_colour_manual(values = group_cols[c("control", "normal_case", "suspected_abnormal")],
    labels = c("Control", "Normal case", "Suspected abnormal"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15),
    labels = c("Control", "Normal case", "Suspected abnormal"), name = NULL) +
  labs(x = df2a$xlab[1], y = df2a$ylab[1]) +
  nature_theme(7) + theme(legend.position = "bottom")

clin <- fread(file.path(BASE, "clinical_649.tsv"))
clin[, sid := as.character(Sample_ID)]
clin[, Class3 := fifelse(tolower(Group1) == "control", "control",
                 fifelse(tolower(Group4) == "abnormal", "abnormal_case", "normal_case"))]
clin[, gw_week := as.integer(Group3)]
clin[, is_like := tolower(sid) %in% like_key]
# Group label for composition within each gestational week
clin[, group := fifelse(is_like, "Suspected abnormal",
                fifelse(Class3 == "control", "Control",
                fifelse(Class3 == "abnormal_case", "abnormal",
                        "Other normal case")))]
# Denominator = all samples in that gestational week
gw_n <- clin[, .(week_n = .N), by = gw_week]
gw_plot <- clin[, .N, by = .(gw_week, group)]
gw_plot <- merge(gw_plot, gw_n, by = "gw_week")
gw_plot[, pct := 100 * N / week_n]
gw_plot[, group := factor(group, levels = c(
  "Control", "Other normal case", "Suspected abnormal", "abnormal"
))]

p2b <- ggplot(gw_plot, aes(factor(gw_week), pct, fill = group)) +
  geom_col(position = "stack", width = 0.75) +
  scale_fill_manual(values = c(
    Control = COL_CTRL,
    `Other normal case` = COL_NORM,
    `Suspected abnormal` = COL_SUS,
    abnormal = COL_ABN
  ), name = NULL) +
  labs(x = "Gestational week", y = "Percent of week") +
  nature_theme(7) + theme(legend.position = "bottom")

fig2 <- (p2a | p2b) + plot_annotation(tag_levels = "a")
save_gg("Fig2", fig2, 180, 85)
fwrite(df2a[, .(sample, Class3, group, PC1, PC2)], file.path(OUT, "Fig2a_source.csv"))
fwrite(gw_plot, file.path(OUT, "Fig2b_source.csv"))

g2 <- assign_group(colnames(nc_mat), nc_cls)
save_heatmap("Fig2", nc_mat, g2,
  sprintf("NC n=%d | clinical abnormal removed | NC top10k | Ward", ncol(nc_mat)),
  pdf_w = 7, pdf_h = 5)

rm(nc_mat, nc)
gc()

# ---------- Fig3 ----------
message("Fig3")
WEEKS <- c("W8", "W9", "W10")
wk_list <- lapply(WEEKS, function(w) {
  wdir <- file.path(BASE, "week8910_hvar", "noabn", w)
  ann <- fread(file.path(wdir, "sample_annotation.tsv"))
  mat <- read_space_matrix(file.path(wdir, "CpG_matrix.tsv"))
  list(week = w, ann = ann, mat = mat)
})

# PCA panel (reuse or recompute)
wk_pca <- rbindlist(lapply(wk_list, function(x) {
  cls <- setNames(x$ann$Class3, x$ann$sample_id)
  df <- pca_df_from_mat(x$mat, cls)
  df[, week := x$week]
  df[, n := ncol(x$mat)]
  df
}))
wk_pca[, panel := factor(
  paste0(c(W8 = "a", W9 = "b", W10 = "c")[week], "  ", week, " (n=", n, ")"),
  levels = paste0(c(W8 = "a", W9 = "b", W10 = "c")[WEEKS], "  ", WEEKS,
                  " (n=", sapply(wk_list, function(x) ncol(x$mat)), ")")
)]
p3 <- ggplot(wk_pca, aes(PC1, PC2, colour = Class3, shape = Class3)) +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_colour_manual(values = c(control = COL_CTRL, normal_case = COL_NORM),
                      labels = c("Control", "Normal case"), name = NULL) +
  scale_shape_manual(values = c(16, 17), labels = c("Control", "Normal case"), name = NULL) +
  facet_wrap(~ panel, nrow = 1, scales = "free") +
  labs(x = "PC1", y = "PC2",
       subtitle = "Clinical abnormal removed; suspected n=30 not in weeks 8–10") +
  nature_theme(7) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 6))
save_gg("Fig3", p3, 180, 75)
fwrite(wk_pca[, .(sample, week, Class3, PC1, PC2)], file.path(OUT, "Fig3_source.csv"))

# Fig3 heatmap: 3 weeks stacked in one PDF
ht_opt$message <- FALSE
pdf(file.path(OUT, "Fig3_heatmap.pdf"), width = 8, height = 12)
for (i in seq_along(wk_list)) {
  x <- wk_list[[i]]
  cls <- setNames(x$ann$Class3, x$ann$sample_id)
  mat <- x$mat
  mat[mat < 0] <- 0
  mat[mat > 1] <- 1
  hc <- hclust(pearson_dist(mat), method = "ward.D2")
  ha <- HeatmapAnnotation(
    Class3 = as.character(cls[colnames(mat)]),
    col = list(Class3 = c(control = COL_CTRL, normal_case = COL_NORM)),
    annotation_name_side = "left"
  )
  draw(Heatmap(
    mat, name = "beta",
    col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE, cluster_columns = hc,
    show_row_names = FALSE, show_column_names = FALSE,
    top_annotation = ha,
    column_title = sprintf("%s n=%d | noabn week-var top10k | Ward",
                           x$week, ncol(mat)),
    row_title = "CpG",
    use_raster = TRUE
  ), newpage = i < length(wk_list))
}
dev.off()
message("Wrote Fig3_heatmap")

# Also single-page PNG for W8 (main week) for PPT slide if needed
x8 <- wk_list[[1]]
cls8 <- setNames(x8$ann$Class3, x8$ann$sample_id)
save_heatmap("Fig3_W8", x8$mat, factor(cls8[colnames(x8$mat)], levels = c("control", "normal_case")),
  sprintf("Week 8 noabn n=%d | clinical abnormal removed", ncol(x8$mat)),
  pdf_w = 7, pdf_h = 5)

# Remove old Fig4 outputs
old4 <- list.files(OUT, pattern = "^Fig4", full.names = TRUE)
if (length(old4)) file.remove(old4)

writeLines(c(
  "PPT slides (Fig1–Fig3 only; Fig4 deleted)",
  "",
  "Fig1     PCA + Fig1_heatmap",
  "         Full cohort; all samples; cohort-var top10k.",
  "",
  "Fig2     PCA + GW + Fig2_heatmap",
  "         Clinical abnormal REMOVED; NC-var top10k.",
  "         Suspected 30 still present (sensitivity test).",
  "",
  "Fig3     PCA + Fig3_heatmap (W8/W9/W10 pages) + optional Fig3_W8_heatmap",
  "         Clinical abnormal REMOVED only (noabn).",
  "         Suspected 30 NOT explicitly removed — they are NOT in weeks 8–10.",
  "         So Fig3 ≈ mainstream control + remaining normal at each week.",
  "",
  "Heatmaps: Ward.D2 on 1-Pearson; dendrogram only (no k-cut colouring).",
  "",
  paste0("Folder: ", OUT)
), file.path(OUT, "PPT_slide_notes.txt"))

message("=== Done ===")
message(OUT)
