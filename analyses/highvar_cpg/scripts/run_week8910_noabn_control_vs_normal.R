#!/usr/bin/env Rscript
# Week-stratified (8/9/10) analysis on noabn high-var CpGs:
# Can control vs normal_case separate within week? Do clusters recur across weeks?
# Figures follow Nature Portfolio artwork conventions (.cursor/skills/nature).

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(ggplot2)
  library(patchwork)
})

BASE <- "D:/ONT/筛选高变CpG"
WEEK_ROOT <- file.path(BASE, "week8910_hvar", "noabn")
OUT <- file.path(BASE, "clustering_output", "week8910_noabn_control_vs_normal")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

WEEKS <- c("W8", "W9", "W10")
K_RANGE <- 2:6

# Okabe–Ito (Nature-recommended)
COL_CTRL <- "#009E73"
COL_NORM <- "#0072B2"
COL_CLUST <- c("#000000", "#E69F00", "#56B4E9", "#CC79A7", "#D55E00", "#F0E442")

MM <- 1 / 25.4
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
      plot.title = element_text(size = base_size, face = "plain", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0),
      plot.margin = margin(2, 2, 2, 2, "mm")
    )
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
    if (length(sp) != n_samp + 1L) {
      sp <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    }
    ids[i] <- sp[1L]
    x <- sp[-1L]
    x[x == "NA"] <- NA
    val[i, ] <- suppressWarnings(as.numeric(x))
  }
  rownames(val) <- ids
  colnames(val) <- mat_samples
  val
}

impute_cpg_mean <- function(mat) {
  out <- mat
  na <- is.na(out)
  rs <- rowMeans(out, na.rm = TRUE)
  out[na] <- rs[row(out)[na]]
  out
}

pearson_dist <- function(mat) {
  cm <- cor(mat, use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  diag(cm) <- 1
  as.dist(1 - cm)
}

fowlkes_mallows <- function(lab1, lab2) {
  lab1 <- as.integer(factor(lab1))
  lab2 <- as.integer(factor(lab2))
  tab <- table(lab1, lab2)
  tk <- sum(choose(tab, 2))
  pk <- sum(choose(rowSums(tab), 2))
  qk <- sum(choose(colSums(tab), 2))
  if (pk == 0 || qk == 0) return(NA_real_)
  sqrt(tk / pk * tk / qk)
}

map_accuracy <- function(cluster, truth) {
  # majority-label accuracy (descriptive; baseline = majority class)
  tab <- table(cluster, truth)
  mapped <- sum(apply(tab, 1, max))
  mapped / sum(tab)
}

majority_baseline <- function(truth) {
  max(table(truth)) / length(truth)
}

analyze_week <- function(week) {
  message("=== ", week, " ===")
  wdir <- file.path(WEEK_ROOT, week)
  ann <- fread(file.path(wdir, "sample_annotation.tsv"))
  ann[, sample_id := as.character(sample_id)]
  mat <- read_space_matrix(file.path(wdir, "CpG_matrix.tsv"))
  stopifnot(identical(colnames(mat), ann$sample_id))
  cls <- factor(ann$Class3, levels = c("control", "normal_case"))
  names(cls) <- ann$sample_id

  d <- pearson_dist(mat)
  hc <- hclust(d, method = "ward.D2")

  sil_rows <- list()
  for (k in K_RANGE) {
    labs <- cutree(hc, k = k)
    sil <- silhouette(labs, d)
    sil_rows[[length(sil_rows) + 1L]] <- data.table(
      week = week, k = k,
      avg_silhouette = mean(sil[, "sil_width"]),
      FM_vs_Class3 = fowlkes_mallows(labs, cls),
      mapped_accuracy = map_accuracy(labs, cls),
      majority_baseline = majority_baseline(cls)
    )
  }
  sil_dt <- rbindlist(sil_rows)
  best_k <- sil_dt$k[which.max(sil_dt$avg_silhouette)]
  labs_best <- cutree(hc, k = best_k)
  labs_k2 <- cutree(hc, k = 2L)

  # PCA for scatter (mean-fill NA only for viz)
  X <- t(impute_cpg_mean(mat))
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  pca_df <- data.table(
    sample = rownames(pc$x),
    Class3 = as.character(cls[rownames(pc$x)]),
    cluster_best = as.integer(labs_best[rownames(pc$x)]),
    cluster_k2 = as.integer(labs_k2[rownames(pc$x)]),
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    week = week,
    pc1_ve = ve[1],
    pc2_ve = ve[2]
  )

  # PCoA of clustering distance
  fit <- cmdscale(d, k = 2, eig = TRUE)
  eig <- fit$eig
  eig[eig < 0] <- 0
  pve <- 100 * eig[1:2] / sum(eig)
  pcoa_df <- data.table(
    sample = rownames(fit$points),
    Class3 = as.character(cls[rownames(fit$points)]),
    cluster_best = as.integer(labs_best[rownames(fit$points)]),
    cluster_k2 = as.integer(labs_k2[rownames(fit$points)]),
    PCoA1 = fit$points[, 1],
    PCoA2 = fit$points[, 2],
    week = week,
    pcoa1_ve = pve[1],
    pcoa2_ve = pve[2]
  )

  xt_best <- as.data.table(as.data.frame.matrix(table(labs_best, cls)), keep.rownames = "cluster")
  xt_best[, week := week]
  xt_k2 <- as.data.table(as.data.frame.matrix(table(labs_k2, cls)), keep.rownames = "cluster")
  xt_k2[, week := week]

  assign_df <- data.table(
    sample = colnames(mat),
    Class3 = as.character(cls),
    week = week,
    cluster_best = as.integer(labs_best),
    cluster_k2 = as.integer(labs_k2),
    best_k = best_k
  )

  list(
    week = week, n = ncol(mat), best_k = best_k,
    sil = sil_dt, pca = pca_df, pcoa = pcoa_df,
    xt_best = xt_best, xt_k2 = xt_k2, assign = assign_df,
    top_cpgs = rownames(mat)
  )
}

results <- lapply(WEEKS, analyze_week)
names(results) <- WEEKS

sil_all <- rbindlist(lapply(results, `[[`, "sil"))
pca_all <- rbindlist(lapply(results, `[[`, "pca"))
pcoa_all <- rbindlist(lapply(results, `[[`, "pcoa"))
assign_all <- rbindlist(lapply(results, `[[`, "assign"))
xt_k2_all <- rbindlist(lapply(results, `[[`, "xt_k2"))
xt_best_all <- rbindlist(lapply(results, `[[`, "xt_best"))

fwrite(sil_all, file.path(OUT, "kscan_silhouette.tsv"), sep = "\t")
fwrite(pca_all, file.path(OUT, "pca_coordinates.tsv"), sep = "\t")
fwrite(pcoa_all, file.path(OUT, "pcoa_coordinates.tsv"), sep = "\t")
fwrite(assign_all, file.path(OUT, "cluster_assignments.tsv"), sep = "\t")
fwrite(xt_k2_all, file.path(OUT, "crosstab_k2.tsv"), sep = "\t")
fwrite(xt_best_all, file.path(OUT, "crosstab_bestk.tsv"), sep = "\t")

# Top-CpG overlap across weeks (Jaccard)
top_lists <- lapply(results, `[[`, "top_cpgs")
ov <- CJ(week_a = WEEKS, week_b = WEEKS)
ov[, `:=`(
  n_a = lengths(top_lists)[week_a],
  n_b = lengths(top_lists)[week_b],
  n_overlap = mapply(function(a, b) length(intersect(top_lists[[a]], top_lists[[b]])), week_a, week_b)
)]
ov[, jaccard := n_overlap / (n_a + n_b - n_overlap)]
fwrite(ov, file.path(OUT, "top10000_overlap_across_weeks.tsv"), sep = "\t")

# Within each week: enrichment of clusters for Class3 (Fisher descriptive table)
enrich_rows <- list()
for (w in WEEKS) {
  a <- results[[w]]$assign
  for (cl in sort(unique(a$cluster_k2))) {
    in_cl <- a$cluster_k2 == cl
    tab <- table(in_cl, a$Class3 == "normal_case")
    # ensure 2x2
    if (all(dim(tab) == c(2, 2))) {
      ft <- fisher.test(tab)
      enrich_rows[[length(enrich_rows) + 1L]] <- data.table(
        week = w, cluster_k2 = cl,
        n = sum(in_cl),
        n_control = sum(in_cl & a$Class3 == "control"),
        n_normal = sum(in_cl & a$Class3 == "normal_case"),
        pct_normal = 100 * mean(a$Class3[in_cl] == "normal_case"),
        fisher_OR = unname(ft$estimate),
        fisher_p = ft$p.value
      )
    }
  }
}
enrich_dt <- rbindlist(enrich_rows)
fwrite(enrich_dt, file.path(OUT, "cluster_k2_class_enrichment.tsv"), sep = "\t")

summary_dt <- sil_all[k == 2, .(
  week, n = sapply(week, function(w) results[[w]]$n),
  best_k_by_sil = sapply(week, function(w) results[[w]]$best_k),
  sil_k2 = avg_silhouette,
  FM_vs_Class3, mapped_accuracy, majority_baseline,
  delta_vs_majority = mapped_accuracy - majority_baseline
)]
summary_dt[, sil_best := {
  mapply(function(w, bk) sil_all[week == w & k == bk, avg_silhouette],
         week, best_k_by_sil)
}]
fwrite(summary_dt, file.path(OUT, "separation_summary.tsv"), sep = "\t")

# -------------------- Figures (Nature) --------------------
pca_all[, week := factor(week, levels = WEEKS)]
pcoa_all[, week := factor(week, levels = WEEKS)]
sil_all[, week := factor(week, levels = WEEKS)]

# Panel labels as facet strip (a–c weeks); Nature lowercase bold via strip
labeller_week <- function(w) {
  idx <- match(as.character(w), WEEKS)
  paste0(letters[idx], "  ", w, " (n=", sapply(as.character(w), function(x) results[[x]]$n), ")")
}
pca_all[, panel := labeller_week(week)]
pcoa_all[, panel := labeller_week(week)]
pca_all[, panel := factor(panel, levels = unique(labeller_week(WEEKS)))]
pcoa_all[, panel := factor(panel, levels = unique(labeller_week(WEEKS)))]

# axis labels with % var (use week-specific; facet shares generic labels via annotate)
# ggplot facet: put VE in subtitle via free scales and custom - use mean note in axis
make_pca_plot <- function(df, colour_by = c("Class3", "cluster_k2")) {
  colour_by <- match.arg(colour_by)
  if (colour_by == "Class3") {
    p <- ggplot(df, aes(PC1, PC2, colour = Class3, shape = Class3)) +
      geom_point(size = 1.1, alpha = 0.85, stroke = 0.2) +
      scale_colour_manual(values = c(control = COL_CTRL, normal_case = COL_NORM),
                          name = NULL) +
      scale_shape_manual(values = c(control = 16, normal_case = 17), name = NULL)
  } else {
    df[, cl := factor(cluster_k2)]
    p <- ggplot(df, aes(PC1, PC2, colour = cl, shape = Class3)) +
      geom_point(size = 1.1, alpha = 0.85, stroke = 0.2) +
      scale_colour_manual(values = COL_CLUST[seq_len(max(df$cluster_k2))], name = "Cluster") +
      scale_shape_manual(values = c(control = 16, normal_case = 17), name = NULL)
  }
  p +
    facet_wrap(~ panel, nrow = 1, scales = "free") +
    labs(x = "PC1", y = "PC2") +
    nature_theme(7) +
    theme(legend.position = "bottom", legend.box = "horizontal")
}

p_pca_class <- make_pca_plot(pca_all, "Class3")
p_pca_cl <- make_pca_plot(pca_all, "cluster_k2")

# Silhouette vs k
p_sil <- ggplot(sil_all, aes(k, avg_silhouette, colour = week, group = week)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = c(W8 = "#000000", W9 = "#E69F00", W10 = "#56B4E9"),
                      name = NULL) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "k", y = "Mean silhouette") +
  nature_theme(7) +
  theme(legend.position = "right")

# Crosstab bars for k=2: fraction normal in each cluster
bar_dt <- enrich_dt[, .(week, cluster_k2, pct_normal, n)]
bar_dt[, week := factor(week, levels = WEEKS)]
p_bar <- ggplot(bar_dt, aes(factor(cluster_k2), pct_normal, fill = week)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c(W8 = "#000000", W9 = "#E69F00", W10 = "#56B4E9"),
                    name = NULL) +
  labs(x = "Cluster (k = 2)", y = "normal_case (%)") +
  nature_theme(7) +
  theme(legend.position = "right")

# Jaccard heatmap of top CpGs
ov_tri <- ov[week_a <= week_b]
p_jac <- ggplot(ov, aes(week_a, week_b, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.2, colour = "black") +
  scale_fill_viridis_c(limits = c(0, 1), name = "Jaccard", option = "C") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  nature_theme(7) +
  theme(legend.position = "right")

fig3 <- (p_sil | p_bar | p_jac) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8, family = "sans"))

save_fig <- function(path_pdf, path_png, plot, w_mm, h_mm) {
  tryCatch({
    ggsave(path_pdf, plot, width = w_mm * MM, height = h_mm * MM,
           device = cairo_pdf, useDingbats = FALSE)
  }, error = function(e) {
    ggsave(path_pdf, plot, width = w_mm * MM, height = h_mm * MM, useDingbats = FALSE)
  })
  ggsave(path_png, plot, width = w_mm * MM, height = h_mm * MM, dpi = 600)
}

save_fig(file.path(OUT, "Fig_week8910_noabn_PCA_Class3.pdf"),
         file.path(OUT, "Fig_week8910_noabn_PCA_Class3.png"),
         p_pca_class, 180, 70)
save_fig(file.path(OUT, "Fig_week8910_noabn_PCA_cluster.pdf"),
         file.path(OUT, "Fig_week8910_noabn_PCA_cluster.png"),
         p_pca_cl, 180, 70)
save_fig(file.path(OUT, "Fig_week8910_noabn_stability.pdf"),
         file.path(OUT, "Fig_week8910_noabn_stability.png"),
         fig3, 180, 70)

# Source data
fwrite(pca_all, file.path(OUT, "source_Fig_PCA_Class3.csv"))
fwrite(sil_all, file.path(OUT, "source_Fig_stability_silhouette.csv"))
fwrite(bar_dt, file.path(OUT, "source_Fig_stability_cluster_pct.csv"))
fwrite(ov, file.path(OUT, "source_Fig_stability_jaccard.csv"))

writeLines(c(
  "Week-stratified noabn analysis (control vs normal_case)",
  "",
  "Narrative: prior full-cohort + NC analyses identified ~30 normal_case as",
  "suspected abnormal (abnormal-like). Those samples are NOT in weeks 8–10.",
  "Here we ask whether, within each of weeks 8/9/10 after dropping clinical",
  "abnormal, control and remaining normal_case separate on week-specific",
  "high-var CpGs, and whether any clusters are stable across weeks.",
  "",
  "Method: Ward.D2 on 1-Pearson; PCA for display (per-CpG mean fill of NA);",
  "no k forced — report silhouette best k and always show k=2 for Class3.",
  "Figures: Nature print width 180 mm, Okabe–Ito colours, vector PDF + 600 dpi PNG.",
  "",
  paste0("Output: ", OUT)
), file.path(OUT, "README.txt"))

cat("\n=== Separation summary (k=2) ===\n")
print(summary_dt)
cat("\n=== Crosstab k=2 ===\n")
print(xt_k2_all)
cat("\n=== Top10000 Jaccard across weeks ===\n")
print(ov)
cat("\nDone: ", OUT, "\n", sep = "")
