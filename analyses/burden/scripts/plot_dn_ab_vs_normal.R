#!/usr/bin/env Rscript
# Plots for DN-background abnormal vs normal analysis (raw-column strategy).
# Usage: Rscript plot_dn_ab_vs_normal.R [table_dir] [plot_dir] [subtitle]
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(ggvenn)
})

root <- "D:/ONT/analyses/burden"
args <- commandArgs(trailingOnly = TRUE)
tab <- if (length(args) >= 1) args[[1]] else file.path(root, "tables/loose_dual_filter")
fig <- if (length(args) >= 2) args[[2]] else file.path(root, "plots/dn_ab_vs_normal")
sub_extra <- if (length(args) >= 3) args[[3]] else "Loose dual-filter; abnormal=77 (47+30 CpG suspects) vs normal"
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

sample_load <- read_tsv(file.path(tab, "dn_sample_load_ab_normal.tsv"), show_col_types = FALSE)
genes <- read_tsv(file.path(tab, "dn_gene_burden_ab_vs_normal.tsv"), show_col_types = FALSE)
sites <- read_tsv(file.path(tab, "dn_site_patterns_ab_vs_normal.tsv"), show_col_types = FALSE)
sens <- read_tsv(file.path(tab, "dn_sensitivity_site_subsets.tsv"), show_col_types = FALSE)

cond_cols <- c(abnormal = "#D55E00", normal = "#0072B2")

# 1) Sample DN site load
p1 <- ggplot(sample_load, aes(Condition, n_dn_sites, fill = Condition)) +
  geom_boxplot(outlier.size = 0.6, width = 0.55, alpha = 0.85) +
  scale_fill_manual(values = cond_cols) +
  labs(
    title = "DN-background SNV load per sample",
    subtitle = sub_extra,
    x = NULL, y = "DN candidate sites per sample"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(fig, "fig_sample_dn_site_load.pdf"), p1, width = 4.2, height = 4.0)
ggsave(file.path(fig, "fig_sample_dn_site_load.png"), p1, width = 4.2, height = 4.0, dpi = 300)

p1b <- ggplot(sample_load, aes(Condition, n_dn_genes, fill = Condition)) +
  geom_boxplot(outlier.size = 0.6, width = 0.55, alpha = 0.85) +
  scale_fill_manual(values = cond_cols) +
  labs(
    title = "DN-background gene hits per sample",
    x = NULL, y = "Genes with >=1 DN candidate"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(fig, "fig_sample_dn_gene_load.pdf"), p1b, width = 4.2, height = 4.0)
ggsave(file.path(fig, "fig_sample_dn_gene_load.png"), p1b, width = 4.2, height = 4.0, dpi = 300)

# 2) Gene volcano (raw Fisher); color by pattern
genes <- genes %>%
  mutate(
    p_ab_num = suppressWarnings(as.numeric(p_ab_greater)),
    p_norm_num = suppressWarnings(as.numeric(p_norm_greater)),
    p_two_num = suppressWarnings(as.numeric(p_fisher_two_sided)),
    p_plot = dplyr::case_when(
      rate_abnormal >= rate_normal & !is.na(p_ab_num) ~ p_ab_num,
      rate_normal > rate_abnormal & !is.na(p_norm_num) ~ p_norm_num,
      !is.na(p_two_num) ~ p_two_num,
      TRUE ~ 1
    ),
    neglog10p = pmin(-log10(pmax(p_plot, 1e-300)), 30),
    log2or = log2(pmax(or_ab_vs_normal, 1e-6))
  )

p2 <- ggplot(genes, aes(log2or, neglog10p, color = pattern)) +
  geom_point(alpha = 0.55, size = 1.1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  labs(
    title = "Gene burden: abnormal vs normal",
    subtitle = "Full-union FDR (two-sided + two one-sided); color = pattern label",
    x = "log2(OR abnormal vs normal)", y = expression(-log[10](P))
  ) +
  theme_bw(base_size = 11)
ggsave(file.path(fig, "fig_gene_volcano.pdf"), p2, width = 6.2, height = 4.5)
ggsave(file.path(fig, "fig_gene_volcano.png"), p2, width = 6.2, height = 4.5, dpi = 300)

# Top genes by full-union one-sided FDR
top_ab <- genes %>%
  filter(pattern == "abnormal_enriched") %>%
  mutate(fdr_use = suppressWarnings(as.numeric(fdr_ab_greater))) %>%
  arrange(fdr_use, p_ab_greater) %>%
  slice_head(n = 20)
top_norm <- genes %>%
  filter(pattern == "normal_enriched") %>%
  mutate(fdr_use = suppressWarnings(as.numeric(fdr_norm_greater))) %>%
  arrange(fdr_use, p_norm_greater) %>%
  slice_head(n = 20)
top <- bind_rows(top_ab, top_norm)
if (nrow(top) > 0) {
  top <- top %>%
    mutate(gene = factor(gene, levels = rev(unique(gene))))
  p3 <- ggplot(top, aes(or_ab_vs_normal, gene, fill = pattern)) +
    geom_col(width = 0.7) +
    scale_x_continuous(trans = "log10") +
    labs(
      title = "Top FDR gene enrichments",
      x = "OR (abnormal vs normal, log10)", y = NULL
    ) +
    theme_bw(base_size = 10)
  ggsave(file.path(fig, "fig_gene_top_or.pdf"), p3, width = 6.5, height = 5.5)
  ggsave(file.path(fig, "fig_gene_top_or.png"), p3, width = 6.5, height = 5.5, dpi = 300)
}

# Exploratory: top nominal p (FDR may be empty with n_ab=47)
top_nom <- genes %>%
  mutate(p_nom = suppressWarnings(as.numeric(p_fisher_two_sided))) %>%
  filter(!is.na(p_nom), p_nom < 0.05) %>%
  arrange(p_nom) %>%
  slice_head(n = 30) %>%
  mutate(
    direction = ifelse(rate_abnormal > rate_normal, "ab_higher", "norm_higher"),
    gene = factor(gene, levels = rev(gene))
  )
if (nrow(top_nom) > 0) {
  p3b <- ggplot(top_nom, aes(or_ab_vs_normal, gene, fill = direction)) +
    geom_col(width = 0.7) +
    scale_x_continuous(trans = "log10") +
    labs(
      title = "Top nominal gene burden (p<0.05, exploratory)",
      subtitle = "No gene passes FDR<0.05; shown for hypothesis generation only",
      x = "OR (abnormal vs normal, log10)", y = NULL
    ) +
    theme_bw(base_size = 10)
  ggsave(file.path(fig, "fig_gene_top_nominal_or.pdf"), p3b, width = 6.5, height = 6.5)
  ggsave(file.path(fig, "fig_gene_top_nominal_or.png"), p3b, width = 6.5, height = 6.5, dpi = 300)
}

# Length-adjusted: rate difference per 100kb among genes with cds
glen <- genes %>%
  filter(!is.na(rate_ab_per_100kb), !is.na(rate_norm_per_100kb)) %>%
  mutate(delta = rate_ab_per_100kb - rate_norm_per_100kb) %>%
  arrange(desc(abs(delta))) %>%
  slice_head(n = 30)
if (nrow(glen) > 0) {
  glen <- glen %>% mutate(gene = factor(gene, levels = rev(gene)))
  p4 <- ggplot(glen, aes(delta, gene, fill = pattern)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    labs(
      title = "Length-adjusted carrier-rate difference",
      subtitle = "rate_ab - rate_norm per 100 kb CDS (top |delta|)",
      x = "Delta carrier rate / 100kb", y = NULL
    ) +
    theme_bw(base_size = 10)
  ggsave(file.path(fig, "fig_gene_length_adj_delta.pdf"), p4, width = 6.5, height = 6.0)
  ggsave(file.path(fig, "fig_gene_length_adj_delta.png"), p4, width = 6.5, height = 6.0, dpi = 300)
}

# 3) Site pattern counts
sp <- sites %>% count(pattern, name = "n")
p5 <- ggplot(sp, aes(reorder(pattern, n), n, fill = pattern)) +
  geom_col() +
  coord_flip() +
  labs(title = "Site pattern counts (DN union)", x = NULL, y = "Sites") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(fig, "fig_site_pattern_counts.pdf"), p5, width = 5.5, height = 3.8)
ggsave(file.path(fig, "fig_site_pattern_counts.png"), p5, width = 5.5, height = 3.8, dpi = 300)

# Site-pattern Venn: sites with >=1 carrier in abnormal vs normal
id_ab <- paste(sites$chrom, sites$pos, sites$ref, sites$alt, sep = ":")[sites$n_abnormal > 0]
id_nm <- paste(sites$chrom, sites$pos, sites$ref, sites$alt, sep = ":")[sites$n_normal > 0]
p_venn <- ggvenn(
  list(abnormal = id_ab, normal = id_nm),
  fill_color = c("#D55E00", "#0072B2"),
  stroke_size = 0.4,
  set_name_size = 4,
  text_size = 4
) +
  labs(
    title = "Site pattern Venn (carrier presence)",
    subtitle = "Sites with >=1 abnormal carrier vs >=1 normal carrier"
  ) +
  theme_void(base_size = 11)
ggsave(file.path(fig, "fig_site_pattern_venn.pdf"), p_venn, width = 5.2, height = 4.6)
ggsave(file.path(fig, "fig_site_pattern_venn.png"), p_venn, width = 5.2, height = 4.6, dpi = 300)

if ("n_control" %in% names(sites) && any(sites$n_control > 0, na.rm = TRUE)) {
  id_ct <- paste(sites$chrom, sites$pos, sites$ref, sites$alt, sep = ":")[sites$n_control > 0]
  p_venn3 <- ggvenn(
    list(abnormal = id_ab, normal = id_nm, control = id_ct),
    fill_color = c("#D55E00", "#0072B2", "#009E73"),
    stroke_size = 0.35,
    set_name_size = 3.5,
    text_size = 3.5
  ) +
    labs(
      title = "Site pattern Venn (including control)",
      subtitle = "Sites with >=1 carrier in each group"
    ) +
    theme_void(base_size = 11)
  ggsave(file.path(fig, "fig_site_pattern_venn_three.pdf"), p_venn3, width = 5.6, height = 5.0)
  ggsave(file.path(fig, "fig_site_pattern_venn_three.png"), p_venn3, width = 5.6, height = 5.0, dpi = 300)
}

# 4) Sensitivity stacked patterns
sens_long <- sens %>%
  tidyr::pivot_longer(
    cols = c(abnormal_enriched, normal_enriched, shared, abnormal_only, normal_only),
    names_to = "pattern", values_to = "n"
  )
p6 <- ggplot(sens_long, aes(subset, n, fill = pattern)) +
  geom_col(position = "stack") +
  labs(
    title = "Sensitivity: site patterns by raw-column subset",
    x = NULL, y = "Sites"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig, "fig_sensitivity_patterns.pdf"), p6, width = 7.0, height = 4.2)
ggsave(file.path(fig, "fig_sensitivity_patterns.png"), p6, width = 7.0, height = 4.2, dpi = 300)

message("Plots written to ", fig)
