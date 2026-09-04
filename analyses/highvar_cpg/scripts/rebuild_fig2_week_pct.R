#!/usr/bin/env Rscript
# Rebuild Fig2 only: GW bars use denominator = samples in that gestational week.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "ppt_figures")
MM <- 1 / 25.4
COL_ABN  <- "#D55E00"
COL_CTRL <- "#009E73"
COL_NORM <- "#0072B2"
COL_SUS  <- "#E69F00"

like_key <- unique(tolower(
  fread("D:/ONT/analyses/highvar_cpg/metadata/abnormal_like_normal_case_30.txt",
        header = FALSE)$V1))

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

assign_group <- function(samples, class3) {
  g <- as.character(class3[samples])
  sus <- tolower(samples) %in% like_key & g == "normal_case"
  g[sus] <- "suspected_abnormal"
  factor(g, levels = c("control", "normal_case", "suspected_abnormal", "abnormal_case"))
}

group_cols <- c(
  control = COL_CTRL, normal_case = COL_NORM,
  suspected_abnormal = COL_SUS, abnormal_case = COL_ABN
)

message("Fig2a PCA")
nc <- readRDS(file.path(BASE, "clustering_output", "gradient_heatmaps", "nc_mat_cache.rds"))
df2a <- pca_df_from_mat(nc$mat, nc$class3)
df2a[, group := assign_group(sample, nc$class3)]
df2a <- df2a[group != "abnormal_case"]

p2a <- ggplot(df2a, aes(PC1, PC2, colour = group, shape = group)) +
  geom_point(data = df2a[group != "suspected_abnormal"], size = 1.0, alpha = 0.65) +
  geom_point(data = df2a[group == "suspected_abnormal"], size = 1.5, alpha = 0.95) +
  scale_colour_manual(
    values = group_cols[c("control", "normal_case", "suspected_abnormal")],
    labels = c("Control", "Normal case", "Suspected abnormal"), name = NULL) +
  scale_shape_manual(
    values = c(16, 17, 15),
    labels = c("Control", "Normal case", "Suspected abnormal"), name = NULL) +
  labs(x = df2a$xlab[1], y = df2a$ylab[1]) +
  nature_theme(7) + theme(legend.position = "bottom")

message("Fig2b week composition")
clin <- fread(file.path(BASE, "clinical_649.tsv"))
clin[, sid := as.character(Sample_ID)]
clin[, Class3 := fifelse(tolower(Group1) == "control", "control",
                 fifelse(tolower(Group4) == "abnormal", "abnormal_case", "normal_case"))]
clin[, gw_week := as.integer(Group3)]
clin[, is_like := tolower(sid) %in% like_key]
clin[, group := fifelse(is_like, "Suspected abnormal",
                fifelse(Class3 == "control", "Control",
                fifelse(Class3 == "abnormal_case", "abnormal",
                        "Other normal case")))]
gw_n <- clin[, .(week_n = .N), by = gw_week]
gw_plot <- clin[, .N, by = .(gw_week, group)]
gw_plot <- merge(gw_plot, gw_n, by = "gw_week")
gw_plot[, pct := 100 * N / week_n]
gw_plot[, group := factor(group, levels = c(
  "Control", "Other normal case", "Suspected abnormal", "abnormal"
))]
setorder(gw_plot, gw_week, group)
print(gw_plot[gw_week == 11])
message("Week11 pct sum = ", gw_plot[gw_week == 11, sum(pct)])

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
ggsave(file.path(OUT, "Fig2.pdf"), fig2, width = 180 * MM, height = 85 * MM, useDingbats = FALSE)
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 180 * MM, height = 85 * MM, dpi = 600)
fwrite(df2a[, .(sample, Class3, group, PC1, PC2)], file.path(OUT, "Fig2a_source.csv"))
fwrite(gw_plot, file.path(OUT, "Fig2b_source.csv"))
message("Done: ", OUT)
