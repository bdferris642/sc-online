#!/usr/bin/env Rscript

print("**************** LOADING LIBRARIES ****************")
# Detect script path when running via Rscript
args = commandArgs(trailingOnly = FALSE)
script_path = sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir = dirname(normalizePath(script_path))
  message("Script located in directory: ", script_dir)
} else {
  stop("Cannot determine script path. Are you running via Rscript?")
}

suppressMessages(suppressWarnings({
  library(dplyr)
  library(getopt)
  library(qs)
}))

spec <- matrix(c(
  'in-dir', 'i', 1, "character",
  'out',   'o', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

IN_DIR = opt[['in-dir']]
OUT_QS = opt[['out']]

if (!dir.exists(dirname(OUT_QS))) {
    dir.create(dirname(OUT_QS), recursive=TRUE)
}

if (is.null(IN_DIR) || is.null(OUT_QS)) {
    stop("Required: --in-dir <dir with per-chunk qs> --out <combined.qs>")
}

files <- list.files(IN_DIR, pattern="\\.qs$", full.names=TRUE)
if (length(files) == 0) stop(sprintf("No .qs files found in %s", IN_DIR))

read_one <- function(path) {
    x = qread(path)
    # Enforce expected structure
    if (is.null(x$summary))         x$summary <- data.frame()
    if (is.null(x$overdispersion))  x$overdispersion <- data.frame()
    if (is.null(x$convergence))     x$convergence <- numeric(0)
    if (is.null(x$algorithm))       x$algorithm <- character(0)
    x
}

lst = lapply(files, read_one)

combined = list(
  summary        = bind_rows(lapply(lst, `[[`, "summary")),
  overdispersion = bind_rows(lapply(lst, `[[`, "overdispersion")),
  convergence    = unlist(lapply(lst, `[[`, "convergence"), use.names=FALSE),
  algorithm      = unlist(lapply(lst, `[[`, "algorithm"), use.names=FALSE),
  covariance     = NULL,
  random_effect  = NULL
)

qsave(combined, OUT_QS)
message("***** Combined NEBULA results written to: ", OUT_QS)
