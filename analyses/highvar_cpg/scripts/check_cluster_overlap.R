# 步骤 1 与步骤 2 的样本重叠：NC k=2 小簇 = 全队列 k=3 簇 3 的非 abnormal 成员。
suppressPackageStartupMessages(library(data.table))
base <- "D:/ONT/筛选高变CpG"
nc <- fread(file.path(base, "clustering_output/NC_matrix_analysis/cluster_assignments_k2.tsv"))
glob <- fread(file.path(base, "clustering_output/cluster_vs_class3.tsv"))

cat("NC names:\n"); print(names(nc))
cat("Global names:\n"); print(names(glob))

setnames(nc, "cluster", "nc_k2", skip_absent = TRUE)
if ("cluster_k2" %in% names(nc)) setnames(nc, "cluster_k2", "nc_k2")

nc2 <- nc[nc_k2 == 2]
cat("NC cluster2 n=", nrow(nc2), "\n")
print(table(nc2$Class3))

g3 <- glob[cluster_k3 == 3]
cat("Global k3 cluster3 n=", nrow(g3), "\n")
print(table(g3$Class3))

ov <- intersect(nc2$sample, g3$sample)
cat("Overlap all samples NC2 vs G3:", length(ov), "/", nrow(nc2), "\n")

nc2_norm <- nc2[Class3 == "normal_case", sample]
g3_norm <- glob[cluster_k3 == 3 & Class3 == "normal_case", sample]
cat("Overlap normal_case:", length(intersect(nc2_norm, g3_norm)),
    " of NC2=", length(nc2_norm), " of G3=", length(g3_norm), "\n")

m <- merge(
  nc[, .(sample, Class3, nc_k2)],
  glob[, .(sample, glob_k3 = cluster_k3)],
  by = "sample", all.x = TRUE
)
cat("\nNC k2 x Global k3 (NC samples only):\n")
print(table(nc_k2 = m$nc_k2, glob_k3 = m$glob_k3, useNA = "ifany"))
cat("\nBy Class3:\n")
print(table(m$Class3, m$nc_k2, m$glob_k3, dnn = c("Class3", "nc_k2", "glob_k3")))
