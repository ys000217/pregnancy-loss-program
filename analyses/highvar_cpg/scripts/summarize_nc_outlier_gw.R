#!/usr/bin/env Rscript
# Gestational week distribution of the NC k=2 outlier normal_case samples (~30).

suppressPackageStartupMessages(library(data.table))

BASE <- "D:/ONT/筛选高变CpG"
OUT  <- file.path(BASE, "clustering_output", "NC_outlier_gestational_week")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

assign <- fread(file.path(BASE, "clustering_output", "NC_matrix_analysis",
                          "cluster_assignments_k2.tsv"))
assign[, sample := as.character(sample)]
clin <- fread("D:/ONT/clinical_649.tsv")
clin[, Sample_ID := as.character(Sample_ID)]

parse_gw <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "na")] <- NA_character_
  vapply(x, function(s) {
    if (is.na(s)) return(NA_real_)
    if (grepl("+", s, fixed = TRUE)) {
      sp <- strsplit(s, "+", fixed = TRUE)[[1]]
      as.numeric(sp[1]) + as.numeric(sp[2]) / 7
    } else {
      as.numeric(s)
    }
  }, numeric(1))
}

clin[, gw := parse_gw(get("Gestational_Week(w)"))]
clin[, gw_week := as.integer(Group3)]
clin[, key := tolower(Sample_ID)]

m <- merge(
  assign[, key := tolower(sample)],
  clin[, .(key, Sample_ID, Group1, Group2, gw, gw_week,
           raw_gw = get("Gestational_Week(w)"), Group4)],
  by = "key", all.x = TRUE, sort = FALSE
)
m[is.na(Sample_ID), Sample_ID := sample]

out30 <- m[Class3 == "normal_case" & cluster == 2]
rest  <- m[Class3 == "normal_case" & cluster == 1]
ctrl  <- m[Class3 == "control"]

summarize_gw <- function(dt, label) {
  x <- dt$gw
  data.table(
    group = label,
    n = nrow(dt),
    n_missing_gw = sum(is.na(x)),
    min = min(x, na.rm = TRUE),
    q25 = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q75 = as.numeric(quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
}

summ <- rbind(
  summarize_gw(out30, "NC_cluster2_normal_case"),
  summarize_gw(rest, "NC_cluster1_normal_case"),
  summarize_gw(ctrl, "NC_control")
)
fwrite(summ, file.path(OUT, "gw_summary.tsv"), sep = "\t")

wt <- wilcox.test(out30$gw, rest$gw)
ks <- ks.test(out30$gw, rest$gw)
tests <- data.table(
  comparison = "cluster2_normal vs cluster1_normal",
  wilcox_W = unname(wt$statistic),
  wilcox_p = wt$p.value,
  ks_D = unname(ks$statistic),
  ks_p = ks$p.value
)
fwrite(tests, file.path(OUT, "gw_tests.tsv"), sep = "\t")

week_tab <- rbind(
  out30[, .(group = "NC_cluster2_normal_case", .N), by = .(gw_week)][, setorder(.SD, gw_week)],
  rest[, .(group = "NC_cluster1_normal_case", .N), by = .(gw_week)][, setorder(.SD, gw_week)]
)
fwrite(week_tab, file.path(OUT, "gw_week_counts.tsv"), sep = "\t")

prop <- dcast(
  rbind(
    out30[!is.na(gw_week), .(group = "outlier30", gw_week)],
    rest[!is.na(gw_week), .(group = "other_normal", gw_week)]
  )[, .N, by = .(group, gw_week)],
  gw_week ~ group, value.var = "N", fill = 0
)
prop[, outlier30_pct := 100 * outlier30 / sum(outlier30)]
prop[, other_normal_pct := 100 * other_normal / sum(other_normal)]
fwrite(prop, file.path(OUT, "gw_week_percent.tsv"), sep = "\t")

fwrite(
  out30[order(gw), .(sample = Sample_ID, Group2, Gestational_Week = raw_gw,
                     gw_numeric = round(gw, 3), gw_week, cluster, Class3)],
  file.path(OUT, "outlier30_samples.tsv"), sep = "\t"
)

# also the 1 control in cluster 2
fwrite(
  m[cluster == 2, .(sample = Sample_ID, Class3, Group2,
                    Gestational_Week = raw_gw, gw_numeric = round(gw, 3), gw_week)],
  file.path(OUT, "cluster2_all31.tsv"), sep = "\t"
)

png(file.path(OUT, "gw_boxplot.png"), width = 1400, height = 1100, res = 160)
par(mar = c(8, 4.5, 3, 1))
boxplot(
  out30$gw, rest$gw, ctrl$gw,
  names = c(
    sprintf("outlier\nnormal\n(n=%d)", sum(!is.na(out30$gw))),
    sprintf("other\nnormal\n(n=%d)", sum(!is.na(rest$gw))),
    sprintf("control\n(n=%d)", sum(!is.na(ctrl$gw)))
  ),
  ylab = "Gestational week",
  main = "Gestational week: NC cluster-2 normal_case vs others",
  col = c("#E41A1C", "#377EB8", "#4DAF4A"),
  las = 1
)
stripchart(
  list(out30$gw, rest$gw, ctrl$gw),
  vertical = TRUE, method = "jitter", jitter = 0.12, add = TRUE,
  pch = 16, col = adjustcolor("black", 0.45), cex = 0.7
)
dev.off()

png(file.path(OUT, "gw_histogram_outlier30.png"), width = 1400, height = 1000, res = 160)
par(mar = c(4.5, 4.5, 3, 1))
hmax <- max(c(out30$gw_week, rest$gw_week), na.rm = TRUE)
hmin <- min(c(out30$gw_week, rest$gw_week), na.rm = TRUE)
br <- (hmin - 0.5):(hmax + 0.5)
hist(out30$gw_week, breaks = br, col = adjustcolor("#E41A1C", 0.7),
     main = "Integer week (Group3) of 30 outlier normal_case",
     xlab = "Gestational week (Group3)", ylab = "Count", las = 1)
dev.off()

png(file.path(OUT, "gw_week_percent_bars.png"), width = 1600, height = 1100, res = 160)
weeks <- sort(unique(c(out30$gw_week, rest$gw_week)))
weeks <- weeks[!is.na(weeks)]
p1 <- sapply(weeks, function(w) mean(out30$gw_week == w, na.rm = TRUE) * 100)
p2 <- sapply(weeks, function(w) mean(rest$gw_week == w, na.rm = TRUE) * 100)
barplot(rbind(p1, p2), beside = TRUE, names.arg = weeks,
        col = c("#E41A1C", "#377EB8"),
        ylab = "Percent of group", xlab = "Gestational week (Group3)",
        main = "Week distribution (%)", las = 1,
        legend.text = c("outlier 30 normal_case", "other normal_case"),
        args.legend = list(x = "topright", bty = "n"))
dev.off()

cat("=== Summary ===\n")
print(summ)
cat("\nWilcoxon / KS vs other normal_case:\n")
print(tests)
cat("\nOutlier 30 week counts (Group3):\n")
print(out30[, .N, by = gw_week][order(gw_week)])
cat("\nOutput: ", OUT, "\n", sep = "")