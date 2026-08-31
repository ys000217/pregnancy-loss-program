#!/usr/bin/env Rscript
# Run full SNV/SV burden pipeline (analyses/burden module)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (!length(file_arg)) stop("Run via: Rscript analyses/burden/scripts/run_all.R")
this_script <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
module <- normalizePath(file.path(dirname(this_script), ".."), winslash = "/", mustWork = TRUE)

data_root <- Sys.getenv("BURDEN_ROOT", unset = "D:/ONT/figure2")
Sys.setenv(
  BURDEN_ROOT = data_root,
  BURDEN_MODULE = module,
  BURDEN_OUT = file.path(module, "tables"),
  BURDEN_PLOT = file.path(module, "plots")
)

scripts <- file.path(module, "scripts")
py_script <- file.path(scripts, "compute_burden.py")
locus_script <- file.path(scripts, "compute_sv_locus_enrichment.py")
plot_script <- file.path(scripts, "analyze_and_plot.R")

py_candidates <- c(
  Sys.which("python"),
  Sys.which("python3"),
  "C:/Users/Administrator/AppData/Local/Programs/Python/Python311/python.exe"
)
python <- py_candidates[nzchar(py_candidates) & file.exists(py_candidates)][1]

to_wsl_path <- function(path) {
  path <- gsub("\\\\", "/", path)
  if (grepl("^[A-Za-z]:", path)) {
    drive <- tolower(substr(path, 1, 1))
    rest <- sub("^[A-Za-z]:", "", path)
    return(paste0("/mnt/", drive, rest))
  }
  path
}

message("Module: ", module)
message("Data root: ", data_root)

message("Step 1/3: compute_burden.py")
if (!is.na(python) && nzchar(python)) {
  status <- system2(python, args = shQuote(py_script))
} else {
  message("Local python not found; trying WSL python3...")
  status <- system2(
    "wsl",
    args = c("bash", "-lc", shQuote(sprintf(
      "export BURDEN_ROOT=%s BURDEN_OUT=%s; python3 %s",
      to_wsl_path(data_root),
      to_wsl_path(file.path(module, "tables")),
      to_wsl_path(py_script)
    )))
  )
}
if (!identical(status, 0L)) {
  stop("compute_burden.py failed with exit code ", status)
}

message("Step 2/3: compute_sv_locus_enrichment.py")
if (!is.na(python) && nzchar(python)) {
  status <- system2(python, args = shQuote(locus_script))
} else {
  status <- system2(
    "wsl",
    args = c("bash", "-lc", shQuote(sprintf(
      "export BURDEN_ROOT=%s BURDEN_OUT=%s; python3 %s",
      to_wsl_path(data_root),
      to_wsl_path(file.path(module, "tables")),
      to_wsl_path(locus_script)
    )))
  )
}
if (!identical(status, 0L)) {
  stop("compute_sv_locus_enrichment.py failed with exit code ", status)
}

message("Step 3/3: analyze_and_plot.R")
source(plot_script)

message("Done. Tables: ", file.path(module, "tables"))
message("plots: ", file.path(module, "plots"))
