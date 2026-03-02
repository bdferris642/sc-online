#!/usr/bin/env Rscript
# install_r_packages.R
# Installs crumblr and its dependencies into cc-sandbox/R_libs.
# Must be run once before using the workflow.

custom_env <- Sys.getenv("R_LIBS_CUSTOM", unset = "")
if (nchar(custom_env) > 0) {
  lib_path <- custom_env
} else {
  args       <- commandArgs(trailingOnly = FALSE)
  file_arg   <- args[grep("--file=", args)]
  script_dir <- if (length(file_arg) > 0)
    dirname(normalizePath(sub("--file=", "", file_arg[1])))
  else
    getwd()
  lib_path <- normalizePath(file.path(script_dir, "..", "R_libs"),
                             mustWork = FALSE)
}

cat("Installing R packages to:", lib_path, "\n")
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_path, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkg, installer = "cran", bioc_version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    tryCatch({
      if (installer == "cran") {
        install.packages(pkg, lib = lib_path, quiet = FALSE,
                         dependencies = TRUE)
      } else if (installer == "bioc") {
        BiocManager::install(pkg, lib = lib_path, ask = FALSE,
                             update = FALSE, force = FALSE)
      } else if (installer == "github") {
        remotes::install_github(bioc_version, lib = lib_path,
                                upgrade = "never", quiet = FALSE)
      }
    }, error = function(e) {
      cat("ERROR installing", pkg, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("Already installed:", pkg, "\n")
  }
}

# ── CRAN dependencies ──────────────────────────────────────────────────────────
for (pkg in c("lme4", "compositions", "pbkrtest", "numDeriv", "minqa",
              "nloptr", "RcppEigen")) {
  install_if_missing(pkg, "cran")
}

# ── Bioconductor dependencies ──────────────────────────────────────────────────
# BiocManager / remotes already in system lib
for (pkg in c("variancePartition", "BiocParallel")) {
  install_if_missing(pkg, "bioc")
}

# ── crumblr from GitHub ────────────────────────────────────────────────────────
install_if_missing("crumblr", "github", "GabrielHoffman/crumblr")

# ── Verify ─────────────────────────────────────────────────────────────────────
cat("\n=== Verification ===\n")
for (pkg in c("lme4", "variancePartition", "crumblr", "compositions")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-25s %s\n", pkg, if (ok) "OK" else "MISSING"))
}
