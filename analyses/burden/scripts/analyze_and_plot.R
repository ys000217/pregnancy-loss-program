# SNV / SV burden: statistics, mut/Mb, and plots

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
})

# Outputs live in analyses/burden; raw inputs stay under figure2 (BURDEN_ROOT).
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  module <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
                          winslash = "/", mustWork = TRUE)
} else {
  module <- Sys.getenv("BURDEN_MODULE", unset = "D:/ONT/analyses/burden")
}

tab_dir <- Sys.getenv("BURDEN_OUT", unset = file.path(module, "tables"))
fig_dir <- Sys.getenv("BURDEN_PLOT", unset = file.path(module, "plots"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

condition_colors <- c(
  abnormal = "#E41A1C",
  normal = "#377EB8",
  control = "#4DAF4A"
)

compute_canonical_exon_mb <- function() {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  tx_by_gene <- transcriptsBy(txdb, by = "gene")
  canonical_tx_ids <- vapply(
    tx_by_gene,
    function(x) x$tx_id[which.max(width(x))],
    FUN.VALUE = integer(1)
  )
  all_exons <- exonsBy(txdb, by = "tx")
  canonical_exons <- all_exons[as.character(canonical_tx_ids)]
  gr <- reduce(unlist(canonical_exons, use.names = FALSE))
  sum(width(gr)) / 1e6
}

normalize_gw <- function(gw) {
  gw <- as.character(gw)
  dplyr::case_when(
    grepl("^g8$|^8$|8\\+", gw, ignore.case = TRUE) ~ "g8",
    grepl("^g9$|^9$|9\\+", gw, ignore.case = TRUE) ~ "g9",
    grepl("^g10$|^10$|10\\+", gw, ignore.case = TRUE) ~ "g10",
    TRUE ~ NA_character_
  )
}

compare_two_groups <- function(x, g, g1, g2) {
  x1 <- x[g == g1]
  x2 <- x[g == g2]
  list(
    group1 = g1,
    group2 = g2,
    n1 = length(x1),
    n2 = length(x2),
    median1 = median(x1),
    median2 = median(x2),
    mean1 = mean(x1),
    mean2 = mean(x2),
    carrier_rate1 = mean(x1 > 0),
    carrier_rate2 = mean(x2 > 0),
    p_wilcox = wilcox.test(x1, x2, exact = FALSE)$p.value
  )
}

run_comparisons <- function(df, metric, label) {
  x <- df[[metric]]
  comps <- list(
    compare_two_groups(x, df$CaseControl, "case", "control"),
    compare_two_groups(x, df$GroupCompare, "abnormal", "control"),
    compare_two_groups(x, df$GroupCompare, "abnormal", "normal"),
    compare_two_groups(x, df$ConditionGroup, "abnormal", "non_abnormal")
  )
  bind_rows(lapply(comps, as.data.frame)) |>
    mutate(metric = label, .before = 1)
}

plot_burden_box <- function(df, y_col, y_label, out_file, title) {
  p <- ggplot(df, aes(x = Condition, y = .data[[y_col]], fill = Condition)) +
    geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.08, alpha = 0.25, size = 0.7) +
    scale_fill_manual(values = condition_colors) +
    labs(title = title, x = NULL, y = y_label) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave(out_file, p, width = 7, height = 5, dpi = 160)
}

plot_group_compare <- function(df, y_col, y_label, out_file, title) {
  p <- ggplot(df, aes(x = GroupCompare, y = .data[[y_col]], fill = GroupCompare)) +
    geom_boxplot(outlier.alpha = 0.35) +
    geom_jitter(width = 0.12, alpha = 0.2, size = 0.7) +
    labs(title = title, x = NULL, y = y_label) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  ggsave(out_file, p, width = 8, height = 5, dpi = 160)
}

plot_carrier_bar <- function(df, metric, out_file, title) {
  summary_df <- df |>
    group_by(Condition) |>
    summarise(
      n = n(),
      carrier_n = sum(.data[[metric]] > 0),
      carrier_rate = carrier_n / n,
      .groups = "drop"
    )
  p <- ggplot(summary_df, aes(x = Condition, y = carrier_rate, fill = Condition)) +
    geom_col(width = 0.65, alpha = 0.9) +
    geom_text(aes(label = sprintf("%d/%d", carrier_n, n)), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = title, x = NULL, y = "Carrier rate") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave(out_file, p, width = 6.5, height = 5, dpi = 160)
}

message("Computing canonical exon Mb (knownGene)...")
coding_mb <- compute_canonical_exon_mb()
message(sprintf("Eligible coding exon size: %.3f Mb", coding_mb))

burden <- read_tsv(file.path(tab_dir, "sample_burden.tsv"), show_col_types = FALSE) |>
  mutate(
    CaseControl = if_else(Condition %in% c("abnormal", "normal"), "case", "control"),
    GroupCompare = Condition,
    ConditionGroup = if_else(Condition == "abnormal", "abnormal", "non_abnormal"),
    GW_bin = normalize_gw(Gestational_Week),
    SNV_mut_per_Mb = SNV_nonsyn_count / coding_mb
  )

write_tsv(
  tibble(coding_mb_canonical_exons = coding_mb, source = "TxDb.Hsapiens.UCSC.hg38.knownGene"),
  file.path(tab_dir, "coding_mb_canonical_exons.tsv")
)

stats_all <- bind_rows(
  run_comparisons(burden, "SNV_nonsyn_count", "SNV_nonsyn_count"),
  run_comparisons(burden, "SNV_mut_per_Mb", "SNV_mut_per_Mb"),
  run_comparisons(burden, "SV_total", "SV_total"),
  run_comparisons(burden, "SV_plp_count", "SV_plp_count"),
  run_comparisons(burden, "SV_plp_rare_count", "SV_plp_rare_count"),
  run_comparisons(burden, "SV_plp_strict_rare_count", "SV_plp_strict_rare_count")
)
write_tsv(stats_all, file.path(tab_dir, "group_comparison_stats.tsv"))

# Sensitivity: g8 / g9 / g10 only (case vs control)
burden_gw <- burden |> filter(!is.na(GW_bin))
stats_gw <- bind_rows(lapply(sort(unique(burden_gw$GW_bin)), function(gw) {
  df <- burden_gw |>
    filter(GW_bin == gw) |>
    mutate(GroupCompare = CaseControl)
  bind_rows(
    compare_two_groups(df$SNV_mut_per_Mb, df$GroupCompare, "case", "control") |>
      as.data.frame() |>
      mutate(metric = "SNV_mut_per_Mb", GW_bin = gw),
    compare_two_groups(df$SV_total, df$GroupCompare, "case", "control") |>
      as.data.frame() |>
      mutate(metric = "SV_total", GW_bin = gw),
    compare_two_groups(df$SV_plp_count, df$GroupCompare, "case", "control") |>
      as.data.frame() |>
      mutate(metric = "SV_plp_count", GW_bin = gw),
    compare_two_groups(df$SV_plp_rare_count, df$GroupCompare, "case", "control") |>
      as.data.frame() |>
      mutate(metric = "SV_plp_rare_count", GW_bin = gw)
  )
}))
write_tsv(stats_gw, file.path(tab_dir, "group_comparison_stats_by_gw.tsv"))

# plots: primary
plot_burden_box(
  burden,
  "SNV_mut_per_Mb",
  sprintf("Nonsynonymous SNV rate (mut/Mb)\nexon denominator = %.2f Mb", coding_mb),
  file.path(fig_dir, "01_SNV_mutMb_by_condition.png"),
  "SNV burden by condition (648 samples)"
)

plot_burden_box(
  burden,
  "SV_total",
  "PASS SV count per sample",
  file.path(fig_dir, "02_SV_total_by_condition.png"),
  "Total SV burden by condition"
)

plot_burden_box(
  burden,
  "SV_plp_count",
  "ACMG P/LP SV count per sample",
  file.path(fig_dir, "03_SV_plp_by_condition.png"),
  "Pathogenic SV burden by condition"
)

case_control <- burden |>
  mutate(GroupCompare = CaseControl)
plot_group_compare(
  case_control,
  "SNV_mut_per_Mb",
  "Nonsynonymous SNV rate (mut/Mb)",
  file.path(fig_dir, "04_SNV_mutMb_case_vs_control.png"),
  "SNV burden: case vs control"
)
plot_group_compare(
  case_control,
  "SV_total",
  "PASS SV count",
  file.path(fig_dir, "05_SV_total_case_vs_control.png"),
  "SV burden: case vs control"
)

plot_carrier_bar(
  burden,
  "SV_plp_count",
  file.path(fig_dir, "06_SV_plp_carrier_rate.png"),
  "Samples carrying >=1 P/LP SV (all)"
)

plot_burden_box(
  burden,
  "SV_plp_rare_count",
  "Rare P/LP SV count per sample\n(no gnomAD AF or pop AF < 1%)",
  file.path(fig_dir, "09_SV_plp_rare_by_condition.png"),
  "Rare P/LP SV burden by condition"
)

plot_group_compare(
  case_control,
  "SV_plp_rare_count",
  "Rare P/LP SV count",
  file.path(fig_dir, "10_SV_plp_rare_case_vs_control.png"),
  "Rare P/LP SV burden: case vs control"
)

plot_carrier_bar(
  burden,
  "SV_plp_rare_count",
  file.path(fig_dir, "11_SV_plp_rare_carrier_rate.png"),
  "Samples carrying >=1 rare P/LP SV"
)

# SV type stacked summary
sv_type_cols <- grep("^SV_", names(burden), value = TRUE)
sv_type_cols <- setdiff(
  sv_type_cols,
  c(
    "SV_total", "SV_plp_count", "SV_plp_binary",
    "SV_plp_rare_count", "SV_plp_rare_binary",
    "SV_plp_strict_rare_count", "SV_plp_strict_rare_binary"
  )
)
if (length(sv_type_cols) > 0) {
  sv_long <- burden |>
    select(Sample_ID, Condition, all_of(sv_type_cols)) |>
    pivot_longer(all_of(sv_type_cols), names_to = "SV_type", values_to = "count") |>
    mutate(SV_type = sub("^SV_", "", SV_type))
  sv_type_summary <- sv_long |>
    group_by(Condition, SV_type) |>
    summarise(mean_count = mean(count), .groups = "drop")
  p_types <- ggplot(sv_type_summary, aes(x = Condition, y = mean_count, fill = SV_type)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Mean SV count by type and condition", x = NULL, y = "Mean SV count") +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_dir, "07_SV_type_composition.png"), p_types, width = 8, height = 5, dpi = 160)
}

# P/LP enrichment top hits
plp <- read_tsv(file.path(tab_dir, "sv_plp_enrichment_case_vs_control.tsv"), show_col_types = FALSE)
plp_top <- plp |>
  arrange(p_value) |>
  slice_head(n = min(20, nrow(plp))) |>
  mutate(label = paste0(coord_key, " (", gene, ")"))
if (nrow(plp_top) > 0) {
  p_enrich <- ggplot(plp_top, aes(x = reorder(label, -p_value), y = -log10(p_value), fill = factor(acmg_class))) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(
      title = "Top P/LP SV enrichment (case vs control)",
      x = NULL,
      y = expression(-log[10](p)),
      fill = "ACMG class"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(fig_dir, "08_SV_plp_enrichment_top20.png"), p_enrich, width = 10, height = 6, dpi = 160)
  write_tsv(plp_top, file.path(tab_dir, "sv_plp_enrichment_top20.tsv"))
}

# Rare P/LP enrichment
plp_rare_path <- file.path(tab_dir, "sv_plp_rare_enrichment_case_vs_control.tsv")
if (file.exists(plp_rare_path)) {
  plp_rare <- read_tsv(plp_rare_path, show_col_types = FALSE)
  plp_rare_top <- plp_rare |>
    arrange(p_value) |>
    slice_head(n = min(20, nrow(plp_rare))) |>
    mutate(label = paste0(coord_key, " (", gene, ")"))
  if (nrow(plp_rare_top) > 0) {
    p_rare <- ggplot(plp_rare_top, aes(x = reorder(label, -p_value), y = -log10(p_value), fill = factor(acmg_class))) +
      geom_col(width = 0.75) +
      coord_flip() +
      labs(
        title = "Rare P/LP SV enrichment (case vs control)",
        subtitle = "Population-rare: gnomAD-SV AF missing or < 1%",
        x = NULL,
        y = expression(-log[10](p)),
        fill = "ACMG class"
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(fig_dir, "12_SV_plp_rare_enrichment_top20.png"), p_rare, width = 10, height = 6, dpi = 160)
    write_tsv(plp_rare_top, file.path(tab_dir, "sv_plp_rare_enrichment_top20.tsv"))
  }
}

# Genome-wide PASS SV locus enrichment (abnormal / normal / control)
plot_locus_enrichment <- function(path, p_col, title, out_png) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }
  df <- read_tsv(path, show_col_types = FALSE)
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }
  df_top <- df |>
    arrange(.data[[p_col]]) |>
    slice_head(n = min(20, nrow(df))) |>
    mutate(label = paste0(coord_key, " (", sv_id, ")"))
  p <- ggplot(df_top, aes(x = reorder(label, -.data[[p_col]]), y = -log10(.data[[p_col]]), fill = svtype)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(title = title, x = NULL, y = expression(-log[10](p)), fill = "SV type") +
    theme_bw(base_size = 11)
  ggsave(out_png, p, width = 10, height = 6, dpi = 160)
}

plot_locus_enrichment(
  file.path(tab_dir, "sv_locus_enrichment_abnormal_specific.tsv"),
  "p_abnormal_vs_control",
  "Abnormal-specific PASS SV loci (also enriched vs normal)",
  file.path(fig_dir, "13_SV_locus_enrichment_abnormal_specific_top20.png")
)
plot_locus_enrichment(
  file.path(tab_dir, "sv_locus_enrichment_top50_abnormal_vs_normal.tsv"),
  "p_abnormal_vs_normal",
  "Top PASS SV loci: abnormal vs normal",
  file.path(fig_dir, "14_SV_locus_enrichment_ab_vs_normal_top20.png")
)
plot_locus_enrichment(
  file.path(tab_dir, "sv_locus_enrichment_top50_normal_vs_control.tsv"),
  "p_normal_vs_control",
  "Top PASS SV loci: normal vs control",
  file.path(fig_dir, "15_SV_locus_enrichment_normal_vs_control_top20.png")
)

message("Analysis and plots complete.")
message(sprintf("Tables: %s", tab_dir))
message(sprintf("Plots: %s", fig_dir))
