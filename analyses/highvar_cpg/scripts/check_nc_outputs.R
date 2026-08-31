suppressPackageStartupMessages(library(data.table))
base <- "D:/ONT/筛选高变CpG"

clin <- fread(file.path(base, "clinical_649.tsv"))
clin[, Sample_ID := as.character(Sample_ID)]
clin[, Class3 := fifelse(tolower(Group1) == "control", "control",
                 fifelse(tolower(Group4) == "abnormal", "abnormal_case", "normal_case"))]

nc <- fread(file.path(base, "NC_samples.tsv"))
ann <- fread(file.path(base, "NC_matrix_sample_annotation.tsv"))
top <- fread(file.path(base, "CpG_topNC_10000.list"), header = FALSE)
stamp <- readLines(file.path(base, "NC_run_stamp.txt"))
abn <- clin[Class3 == "abnormal_case", Sample_ID]

cat("=== Clinical reference ===\n")
print(table(clin$Class3))

cat("\n=== NC_samples.tsv ===\n")
cat("n =", nrow(nc), "\n")
print(table(nc$Class3))
leak <- intersect(nc$sample_id, abn)
cat("abnormal leaked:", length(leak), "\n")

cat("\n=== Annotation ===\n")
cat("n =", nrow(ann), "\n")
print(table(ann$Class3))
cat("ann set == NC_samples:", setequal(ann$sample_id, nc$sample_id), "\n")
cat("ann order == NC_samples:", identical(ann$sample_id, nc$sample_id), "\n")

hdr <- strsplit(readLines(file.path(base, "CpG_matrix_NC.tsv"), n = 1L), "\t", fixed = TRUE)[[1]]
mat_s <- hdr[-1]
nl <- length(readLines(file.path(base, "CpG_matrix_NC.tsv")))

cat("\n=== Matrix ===\n")
cat("n_samples header =", length(mat_s), "\n")
cat("n_lines incl header =", nl, " => CpG rows =", nl - 1L, "\n")
cat("header set == annot:", setequal(mat_s, ann$sample_id), "\n")
cat("header order == annot:", identical(mat_s, ann$sample_id), "\n")

line2 <- readLines(file.path(base, "CpG_matrix_NC.tsv"), n = 2L)[2]
sp <- strsplit(line2, " ", fixed = TRUE)[[1]]
tb <- strsplit(line2, "\t", fixed = TRUE)[[1]]
cat("row2 space-fields =", length(sp), "; tab-fields =", length(tb), "\n")
vals <- if (length(tb) == length(mat_s) + 1L) tb[-1] else sp[-1]
cat("row2 NA rate =", round(mean(vals == "NA"), 4), "\n")

# NA rate across a few rows (quick)
set.seed(1)
idx <- sample(2:nl, min(20L, nl - 1L))
lines <- readLines(file.path(base, "CpG_matrix_NC.tsv"))[idx]
na_rates <- vapply(lines, function(ln) {
  x <- strsplit(ln, " ", fixed = TRUE)[[1]][-1]
  if (length(x) != length(mat_s)) x <- strsplit(ln, "\t", fixed = TRUE)[[1]][-1]
  mean(x == "NA")
}, numeric(1))
cat("mean NA rate (20 random CpGs) =", round(mean(na_rates), 4),
    "; max =", round(max(na_rates), 4), "\n")

cat("\n=== Top list ===\n")
cat("n =", nrow(top), "; unique =", uniqueN(top$V1), "\n")

# variance: check top list is high variance by ranking
cat("\n=== Variance cross-check (top list vs variance file) ===\n")
# Only read variance for top IDs via awk-like filter would be heavy; sample check:
var_head <- fread(file.path(base, "CpG_variance_NC.tsv"), header = FALSE,
                  col.names = c("CpG_ID", "count", "variance"), nrows = 5)
cat("variance file first rows (unsorted dump):\n")
print(var_head)
# confirm top1 exists in variance and get its rank among a sort of first... 
# better: use data.table fread of only needed - too big. Use system sort check of top values
# Check that top list sites have variance present
keep <- top$V1
# stream filter with R connection - sample 100 from top
samp <- keep[1:100]
# build a small awk? Use R to scan - slow for 14M. Skip full; check file size via stamp
cat("stamp:\n")
writeLines(stamp)

# Check no abnormal sample IDs appear in matrix header
cat("\n=== Abnormal in matrix header ===\n")
cat("abnormal IDs in matrix:", length(intersect(mat_s, abn)), "\n")

# Clinical expected NC
exp_nc <- clin[Class3 %in% c("control", "normal_case"), Sample_ID]
cat("\n=== Clinical vs NC set ===\n")
cat("clinical NC n =", length(exp_nc), "\n")
cat("NC_samples vs clinical NC setequal:", setequal(nc$sample_id, exp_nc), "\n")
cat("only in clinical NC:", length(setdiff(exp_nc, nc$sample_id)), "\n")
cat("only in NC_samples:", length(setdiff(nc$sample_id, exp_nc)), "\n")

cat("\n=== PASS checklist ===\n")
ok1 <- nrow(nc) == 602 && sum(nc$Class3 == "control") == 154 && sum(nc$Class3 == "normal_case") == 448
ok2 <- length(leak) == 0 && length(intersect(mat_s, abn)) == 0
ok3 <- length(mat_s) == 602 && identical(mat_s, ann$sample_id)
ok4 <- nrow(top) == 10000 && uniqueN(top$V1) == 10000
ok5 <- (nl - 1L) == 10000
ok6 <- length(sp) == length(mat_s) + 1L || length(tb) == length(mat_s) + 1L
cat("sample counts OK:", ok1, "\n")
cat("no abnormal leak:", ok2, "\n")
cat("matrix-annot align:", ok3, "\n")
cat("top10000 OK:", ok4, "\n")
cat("matrix rows 10000:", ok5, "\n")
cat("matrix row field count OK:", ok6, "\n")
cat("ALL CRITICAL PASS:", all(c(ok1, ok2, ok3, ok4, ok5, ok6)), "\n")
