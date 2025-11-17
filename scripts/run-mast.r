# Runs MAST on a Seurat object or SCE in a .qs or .RDS file
# Notes:
# - Expects *single-cell* data (cells = columns).
# - Input counts must be raw integer-like counts (UMIs). Internally we compute log2(CPM+1) for MAST.
# - You can pseudobulk yourself before running; this script is aimed at per-cell MAST.
#
# CLI:
# --path: path to Seurat object or SCE in .qs or .RDS format
# --cov-list: a comma-separated string of covariates for the design (e.g. case_control,sex,age)
# --suffix: a suffix appended to output filenames (defaults to ISO date)
# --min-counts-gene: minimum total counts across all cells for a gene (default: 0)
# --min-frac-gene: minimum fraction of cells with nonzero counts for a gene (default: 0.01)

########################################## IMPORTS ###########################################

args = commandArgs(trailingOnly = FALSE)
script_path = sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
    script_dir = dirname(normalizePath(script_path))
    setwd(script_dir)
    message("Working directory set to: ", script_dir)
} else {
    stop("Cannot determine script path. Are you running via Rscript?")
}

print("Loading libraries...")
suppressWarnings(suppressMessages({
    library(dplyr)
    library(getopt)
    library(glue)
    library(Matrix)
    library(MAST)
    library(qs)
    library(Seurat)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(tidyr)
    # user utilities expected alongside run-deseq2.R
    source("../de.R", local = TRUE)
    source("../utils.R", local = TRUE)
}))
g = glue::glue

########################################### ARGUMENTS & CONSTANTS ###########################################

spec = matrix(c(
    "path", "p", 1, "character",
    "cov-list", "cl", 1, "character",
    "suffix", "s", 2, "character",
    "min-counts-gene", "mcg", 2, "integer",
    "min-frac-gene", "mfg", 2, "double",
    'numi-col', 'n', 2, "character",
    "out-dir", "o", 2, "character",
    "n-cores", "nc", 2, "integer"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)

print(spec)
print(opt)

PATH = opt[["path"]]
if (is.null(PATH)) stop("--path is required")

COV_LIST = unlist(strsplit(opt[["cov-list"]], split = ",", fixed = TRUE))
if (length(COV_LIST) == 0) stop("--cov-list must include at least one column name")

if (is.null(opt[["suffix"]])) {
    SUFFIX = ""
} else {
    SUFFIX = g('__{opt[["suffix"]]}')
}

if (is.null(opt[["min-counts-gene"]])) {
    MIN_COUNTS_GENE = 0L
} else {
    MIN_COUNTS_GENE = as.integer(opt[["min-counts-gene"]])
}

if (is.null(opt[["min-frac-gene"]])) {
    MIN_FRAC_GENE = 0.01
} else {
    MIN_FRAC_GENE = as.numeric(opt[["min-frac-gene"]])
}

if (is.null(opt[["numi-col"]])) {
    NUMI_COL = "nCount_RNA"
} else {
    NUMI_COL = opt[["numi-col"]]
}

if (is.null(opt[["out-dir"]])){
    # strip .qs/.rds/.RDS suffix
    SLOGAN = gsub("\\.(qs|rds)$", "", basename(PATH), ignore.case = TRUE)
    OUT_DIR = g("{dirname(PATH)}/de_results/mast/{SLOGAN}{SUFFIX}")
} else {
    OUT_DIR = opt[["out-dir"]]
}

if (!dir.exists(OUT_DIR)) {
    dir.create(OUT_DIR, recursive = TRUE)
}

PARAMS = list(
    path = PATH,
    covariates = paste(COV_LIST, collapse = ","),
    suffix = SUFFIX,
    min_counts_gene = MIN_COUNTS_GENE,
    min_frac_gene = MIN_FRAC_GENE
)
write_list(PARAMS, g("{OUT_DIR}/MAST__{SLOGAN}{SUFFIX}__params.txt"))
message("Input Params:")
print(PARAMS)

########################################### MAIN  ###########################################

message(paste0("Loading object from ", PATH))
obj = load_obj(PATH)

# Pull counts and metadata
if (inherits(obj, "Seurat")) {
    counts = GetAssayData(obj, assay = "RNA", slot = "counts")
    md = obj@meta.data
} else if (inherits(obj, "SingleCellExperiment")) {
    message("Object is a SingleCellExperiment")
    counts = SingleCellExperiment::counts(obj)
    md = as.data.frame(colData(obj))
} else {
    stop("Object is neither a Seurat object nor a SingleCellExperiment")
}

# Align metadata to counts columns
if (!all(colnames(counts) %in% rownames(md))) {
    stop("Not all count column names are present as rownames in metadata.")
}
md = md[colnames(counts), , drop = FALSE]

# Build design formula
DESIGN_FORMULA = as.formula(paste0("~ ", paste(COV_LIST, collapse = " + ")))
message(g("Running MAST with formula {DESIGN_FORMULA}..."))

# Run MAST
mast_results = run_mast(
    counts = counts,
    md = md,
    design_formula = DESIGN_FORMULA,
    min_counts_gene = MIN_COUNTS_GENE,
    min_frac_gene = MIN_FRAC_GENE,
    nUMI_col = NUMI_COL
)

message(paste0("Writing MAST results to ", OUT_DIR))
for (coef in names(mast_results)) {
    out_path = g("{OUT_DIR}/MAST__{SLOGAN}{SUFFIX}__{coef}.csv")
    write.csv(mast_results[[coef]], file = out_path, row.names = TRUE)
}
message("Done.")