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
  library(Seurat)
  library(glue)
  g=glue::glue

  source(file.path(script_dir, "../add_rank_purity.R"))
}))

spec <- matrix(c(
  'in-dir', 'i', 1, "character",
  'out',   'o', 1, "character",
  'seurat-qs', 'q', 1, "character",
  'contrast-col', 'c', 1, "character",
  'covariates', 'v', 1, "character",
  'id-col', 'd', 1, "character",
  'suffix', 's', 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

IN_DIR = opt[['in-dir']]
OUT_QS = opt[['out']]
OUT_CSV = sub("\\.qs$", ".csv", OUT_QS)
SEURAT_QS_PATH = opt[['seurat-qs']]
CONTRAST_COL = opt[['contrast-col']]
COVARIATES = strsplit(opt[['covariates']], ",")[[1]]
ID_COL = opt[['id-col']]
SUFFIX = opt[['suffix']]

if (!dir.exists(dirname(OUT_QS))) {
    dir.create(dirname(OUT_QS), recursive=TRUE)
}

if (is.null(IN_DIR) || is.null(OUT_QS)) {
    stop("Required: --in-dir <dir with per-chunk qs> --out <combined.qs>")
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

# do R equivalent of python's “%.3f”
# summary_df[] <- lapply(summary_df, function(x) if(is.numeric(x)) formatC(x, format = "f", digits = 3) else x)
logFC_cols = grep("logFC", colnames(summary_df), value=TRUE)
logFC_df = summary_df %>% select(all_of(c("gene", logFC_cols)))
colnames(logFC_df) = gsub("^logFC_", "", colnames(logFC_df))

message(g("loading seurat object from {SEURAT_QS_PATH}"))
seurat_orig = qread(SEURAT_QS_PATH)
seurat_orig@meta.data$region_SN = ifelse(grepl("SN", seurat_orig@meta.data$region), 1, 0)
seurat_orig@meta.data$region_VTA = ifelse(grepl("VTA", seurat_orig@meta.data$region), 1, 0)
message(g("computing rank purity..."))

rank_purity_df = compute_rank_purity(
    seurat_obj = seurat_orig,
    coef_df = logFC_df,
    covariates = COVARIATES,
    contrast_col = CONTRAST_COL,
    groupby_col = NULL,
    donor_col = ID_COL,
    case_value = "pd"
)

# join to summary_df on gene 
summary_df = summary_df %>%
    left_join(rank_purity_df, by = "gene")

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
