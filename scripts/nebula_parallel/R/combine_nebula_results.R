#!/usr/bin/env Rscript
# Combine per-chunk NEBULA results into one result list (.qs) + summary table (.csv).
#
# Copy of ~/sc-online/scripts/combine-nebula-results.R, vendored into run_nebula_parallel.
# Changes vs upstream: self-activation below; the rank-purity branch (and its add_rank_purity.R
# source, Seurat load and --seurat-qs/--contrast-col/--covariates/--id-col/--case-value/
# --get-rank-purity options) is removed, because rank purity is computed from the h5ad by
# python/rank_purity.py and joined by add_rank_purity_to_combined.R. The merge, padj,
# -log10 and signed -log10 columns are unchanged.

# ---- self-activation: re-exec under the workflow env Rscript (env path derived from this file)
local({
    args_all = commandArgs(trailingOnly = FALSE)
    script = normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
    env_prefix = file.path(dirname(dirname(dirname(script))), "micromamba_root", "envs", "run-nebula-parallel-env")
    if (!dir.exists(env_prefix)) stop("workflow env not found: ", env_prefix, " (run run_nebula_parallel/setup.sh)")
    if (normalizePath(R.home()) != normalizePath(file.path(env_prefix, "lib", "R"))) {
        status = system2(file.path(env_prefix, "bin", "Rscript"), c(shQuote(script), shQuote(commandArgs(trailingOnly = TRUE))))
        quit(status = status, save = "no")
    }
})

suppressMessages(suppressWarnings({
  library(dplyr)
  library(getopt)
  library(qs)
  library(glue)
  g=glue::glue
}))

spec <- matrix(c(
  'in-dir', 'i', 1, "character",
  'out',   'o', 1, "character",
  'suffix', 's', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

IN_DIR = opt[['in-dir']]
OUT_QS = opt[['out']]
SUFFIX = opt[['suffix']]

if (is.null(IN_DIR) || is.null(OUT_QS)) {
    stop("Required: --in-dir <dir with per-chunk qs> --out <combined.qs>")
}
OUT_CSV = sub("\\.qs$", ".csv", OUT_QS)

if (!dir.exists(dirname(OUT_QS))) {
    dir.create(dirname(OUT_QS), recursive=TRUE)
}

files <- list.files(IN_DIR, pattern="\\.qs$", full.names=TRUE)
if (!is.null(SUFFIX) && SUFFIX != "") {
    files <- files[grepl(SUFFIX, files)]
}

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

summary_df = bind_rows(lapply(lst, `[[`, "summary"))

# put gene first
col_order = c("gene", setdiff(colnames(summary_df), "gene"))
summary_df = summary_df %>% select(all_of(col_order))

p_cols = grep("^p_", colnames(summary_df), value=TRUE)
for (col in p_cols) {
    padj_col = gsub("^p_", "padj_", col)
    if (!(padj_col %in% colnames(summary_df))) {
        summary_df[[padj_col]] = p.adjust(summary_df[[col]], method = "BH")
    }
}
p_cols = grep("^(p_|padj_)", colnames(summary_df), value=TRUE)
for (col in p_cols) {
    neg_log10_col = paste0("neg_log10_", col)
    if (!(neg_log10_col %in% colnames(summary_df))) {
        summary_df[[neg_log10_col]] = -log10(summary_df[[col]] + 1e-30)
    }
    signed_neg_log10_col = paste0("signed_neg_log10_", col)
    if (!(signed_neg_log10_col %in% colnames(summary_df))) {
        # suffix is what comes after p_ or padj_
        suffix = sub("^(p_|padj_)", "", col)
        summary_df[[signed_neg_log10_col]] = sign(summary_df[[paste0("logFC_", suffix)]]) * summary_df[[neg_log10_col]]
    }
}

combined = list(
  summary        = summary_df,
  overdispersion = bind_rows(lapply(lst, `[[`, "overdispersion")),
  convergence    = unlist(lapply(lst, `[[`, "convergence"), use.names=FALSE),
  algorithm      = unlist(lapply(lst, `[[`, "algorithm"), use.names=FALSE),
  covariance     = NULL,
  random_effect  = NULL
)
write.csv(summary_df, OUT_CSV)
message(g("***** Combined NEBULA summary written to: {OUT_CSV}"))

qsave(combined, OUT_QS)
message(g("***** Combined NEBULA results written to: {OUT_QS}"))
