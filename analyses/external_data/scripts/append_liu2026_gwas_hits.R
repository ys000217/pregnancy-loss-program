# Curate Liu 2026 Nat Genet gestational-phenotype GWAS hits from GWAS Catalog
# (not a miscarriage case–control GWAS; East Asian pregnancy context).
# Does not dump all ~4,688 independent signals.

Sys.setenv(R_LIBS = "F:/R/library")
.libPaths("F:/R/library")

raw <- "D:/ONT/analyses/external_data/data/raw"
meta <- "D:/ONT/analyses/external_data/metadata"

assoc <- read.csv(file.path(raw, "liu2026_assoc_raw.csv"), stringsAsFactors = FALSE)
assoc$rsid <- sub("-.*$", "", assoc$risk_allele)
assoc$effect_allele <- sub("^.*-", "", assoc$risk_allele)
assoc$p <- as.numeric(assoc$p)

coords <- read.csv(file.path(raw, "liu2026_snp_coords.csv"), stringsAsFactors = FALSE)
extra <- data.frame(
  rsid = c("rs4764725", "rs4838178", "rs4842131", "rs55730982",
           "rs60846954", "rs60904764", "rs61342118"),
  alleles = c("T/A/C/G", "C/A/G/T", "T/A/C", "T/C/G", "C/A/T", "G/A", "C/A"),
  chrom38 = c("12", "9", "9", "3", "14", "4", "7"),
  pos38 = c(103118536L, 124458356L, 136200833L, 157079419L,
            101474570L, 73404041L, 128121865L),
  stringsAsFactors = FALSE
)
coords <- rbind(coords[, c("rsid", "alleles", "chrom38", "pos38")], extra)

hg19 <- read.table(file.path(raw, "liu2026_hg19.bed"), stringsAsFactors = FALSE)
colnames(hg19) <- c("chr", "start", "end", "rsid")
ex19 <- read.table(file.path(raw, "liu2026_extra_hg19.bed"), stringsAsFactors = FALSE)
colnames(ex19) <- c("chr", "start", "end", "rsid")
hg19 <- rbind(hg19, ex19)
hg19$pos_hg19 <- hg19$end

study <- data.frame(
  study = c("GCST90837243", "GCST90837218", "GCST90837297", "GCST90837216",
            "GCST90837293", "GCST90837229", "GCST90837296", "GCST90837262",
            "GCST90837270", "GCST90837271", "GCST90837317"),
  phenotype = c("gestational_diabetes_mellitus", "gestational_anemia",
                "subclinical_hypothyroidism_in_pregnancy", "albumin_in_pregnancy",
                "systolic_bp_in_pregnancy", "diastolic_bp_in_pregnancy",
                "subclinical_hyperthyroxinemia_in_pregnancy",
                "isolated_hypothyroxinemia_in_pregnancy",
                "mild_anemia_in_pregnancy", "severe_or_moderate_anemia_in_pregnancy",
                "unspecified_anemia_in_pregnancy"),
  n_case = c(11670, 41134, 4961, 97444, 78535, 77868, 4539, 3634, 29438, 11696, 11240),
  n_ctrl = c(65896, 55204, 69084, NA, NA, NA, 69084, 47214, 55204, 55204, 10678),
  quantitative = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  priority = c("P0", "P0", "P0", "P0", "P1", "P1", "P0", "P0", "P0", "P0", "P0"),
  stringsAsFactors = FALSE
)

gene_map <- c(
  rs10830962 = "MTNR1B", rs11020106 = "MTNR1B", rs117464115 = "MTNR1B",
  rs271045 = "MTNR1B", rs9368222 = "CDKAL1", rs4237150 = "GLIS3",
  rs1573051 = "HHEX", rs61342118 = "SND1",
  rs9411372 = "ABO", rs1053878 = "ABO",
  rs76267595 = "HBA", rs76038336 = "HBA", rs116875075 = "HBA",
  rs3760052 = "HBA", rs181128092 = "HBA",
  rs2413450 = "HBB", rs2076086 = "HBB",
  rs737310 = "FN1", rs13015993 = "FN1",
  rs74678088 = "NTNG1", rs2125747 = "NTNG1",
  rs10983700 = "C9orf3",
  rs4764725 = "C12orf42",
  rs10857147 = "FGF5", rs13009997 = "GPD2",
  rs2077218 = "CYP17A1", rs6538195 = "ATP2B1", rs12579302 = "ATP2B1",
  rs1348004 = "NTRK3", rs10213867 = "TERB2", rs17767383 = "MAF",
  rs10443230 = "GLIS1", rs190595112 = "FGF2",
  rs4838178 = "NR5A1", rs4842131 = "ABO",
  rs145816106 = "SLCO1C1", rs60846954 = "MEG3",
  rs7883218 = "CHRDL1", rs889761 = "FANCA",
  rs60904764 = "ADAMTS3"
)

other_allele <- function(alleles, ea) {
  parts <- unlist(strsplit(alleles, "/", fixed = TRUE))
  rest <- setdiff(parts, ea)
  if (length(rest) == 0) return(NA_character_)
  rest[1]
}

keep <- merge(assoc, study, by = "study")
keep$pnum <- keep$p
# GWS only; albumin restricted to paper's gestation-specific example
keep <- keep[keep$pnum <= 5e-8, ]
keep <- keep[!(keep$study == "GCST90837216" & keep$rsid != "rs4764725"), ]

keep <- merge(keep, coords, by = "rsid", all.x = TRUE)
keep <- merge(keep, hg19[, c("rsid", "pos_hg19")], by = "rsid", all.x = TRUE)
keep$nearest_gene <- unname(gene_map[keep$rsid])
keep$nearest_gene[is.na(keep$nearest_gene)] <- "NA"
keep$other <- mapply(other_allele, keep$alleles, keep$effect_allele)

keep$hit_id <- paste0("LIU2026_", toupper(gsub("_in_pregnancy|_mellitus", "", keep$phenotype)), "_", keep$rsid)
keep$hit_id <- gsub("GESTATIONAL_DIABETES", "GDM", keep$hit_id)
keep$hit_id <- gsub("GESTATIONAL_ANEMIA", "GANEMIA", keep$hit_id)
keep$hit_id <- gsub("SUBCLINICAL_HYPOTHYROIDISM", "SCH", keep$hit_id)
keep$hit_id <- gsub("SUBCLINICAL_HYPERTHYROXINEMIA", "SCHT4", keep$hit_id)
keep$hit_id <- gsub("ISOLATED_HYPOTHYROXINEMIA", "IH", keep$hit_id)
keep$hit_id <- gsub("SYSTOLIC_BP", "SBP", keep$hit_id)
keep$hit_id <- gsub("DIASTOLIC_BP", "DBP", keep$hit_id)
keep$hit_id <- gsub("MILD_ANEMIA", "MILDANEMIA", keep$hit_id)
keep$hit_id <- gsub("SEVERE_OR_MODERATE_ANEMIA", "MODSEVANEMIA", keep$hit_id)
keep$hit_id <- gsub("UNSPECIFIED_ANEMIA", "UNSPANEMIA", keep$hit_id)
keep$hit_id <- gsub("ALBUMIN", "ALB", keep$hit_id)

keep$n_case_out <- ifelse(keep$quantitative, keep$n_case, keep$n_case)
keep$n_ctrl_out <- ifelse(keep$quantitative, NA, keep$n_ctrl)

fmt_p <- function(x) formatC(x, format = "e", digits = 1)
fmt_b <- function(x) formatC(as.numeric(x), format = "f", digits = 4)
fmt_maf <- function(x) formatC(as.numeric(x), format = "f", digits = 4)

keep$ci <- gsub("\\[|\\]", "", keep$ci)
keep$notes <- paste0(
  "Liu2026 Nat Genet gestational phenotypes (Chinese NIPT/EHR; not miscarriage GWAS). ",
  keep$study, ". Catalog-reported lead (subset of ~4688 GWS signals). ",
  ifelse(keep$rsid == "rs4764725",
         "Paper example of gestation-specific albumin locus (Ext. Data Fig. 5). ",
         ""),
  "PheWeb https://monn.pheweb.com/ ; Primary coords GRCh38 (Ensembl/GWAS Catalog); pos_hg19_legacy from hg38ToHg19 liftOver."
)

out <- data.frame(
  hit_id = keep$hit_id,
  priority = keep$priority,
  rsid = keep$rsid,
  chrom = keep$chrom38,
  pos = keep$pos38,
  genome_build = "GRCh38",
  pos_hg19_legacy = keep$pos_hg19,
  effect_allele = keep$effect_allele,
  other_allele = keep$other,
  effect_metric = "beta",
  effect_value = fmt_b(keep$beta),
  effect_ci = keep$ci,
  p_value = fmt_p(keep$pnum),
  p_note = paste0("GWAS Catalog ", keep$study, " reported"),
  maf_or_raf = fmt_maf(keep$maf),
  phenotype = keep$phenotype,
  ancestry = "Chinese_EAS",
  n_case = keep$n_case_out,
  n_ctrl = keep$n_ctrl_out,
  nearest_gene = keep$nearest_gene,
  study_short = "Liu2026",
  pmid = "42509370",
  doi = "10.1038/s41588-026-02677-w",
  citation = "Liu S et al. Nat Genet. 2026;58:1845-1854. PMID:42509370",
  notes = keep$notes,
  stringsAsFactors = FALSE
)

out <- out[order(out$priority, out$phenotype, out$p_value), ]
stopifnot(!any(is.na(out$chrom)), !any(is.na(out$pos)))

old <- read.delim(file.path(meta, "gwas_hits.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
old <- old[!grepl("^LIU2026_", old$hit_id), ]
# write.table on mixed types; coerce all to character for TSV
to_chr <- function(d) {
  for (i in seq_len(ncol(d))) d[[i]] <- as.character(d[[i]])
  d
}
old$n_ctrl <- as.character(old$n_ctrl)
old$n_case <- as.character(old$n_case)
out$n_ctrl[is.na(out$n_ctrl)] <- "NA"
out$n_case <- as.character(out$n_case)
combined <- rbind(to_chr(old), to_chr(out))
write.table(combined, file.path(meta, "gwas_hits.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat("wrote gwas_hits n=", nrow(combined), " liu rows=", nrow(out), "\n")
print(out[, c("hit_id", "rsid", "phenotype", "p_value", "chrom", "pos")])
