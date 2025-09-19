Sys.setenv(RETICULATE_PYTHON = "/home/ferris/anaconda3/envs/r-4.3.0/bin/python")
suppressMessages(suppressWarnings(library(anndata)))
suppressMessages(suppressWarnings(library(reticulate)))

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SeuratDisk)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(source("~/sc-online/utils.R")))

#use_python("/home/ferris/anaconda3/envs/r-4.3.0/bin/python", required = TRUE) # OJO have to change this path to conda path

# Parse command line arguments
spec = matrix(c(
    'path', 'p', 1, "character",
    'min-features', 'f', 1, "numeric",
    'min-cells', 'c', 1, "numeric",
    'convert-rownames', 'r', 1, "logical",
    'use-raw', 'u', 0, "logical",
    'clobber', 'x', 0, "logical"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
PATH = opt[['path']]
MIN_FEATURES = ifelse(
    is.null(opt[['min-features']]), 
    0,
    opt[['min-features']]
)
MIN_CELLS = ifelse(
    is.null(opt[['min-cells']]), 
    0,
    opt[['min-cells']]
)

if (is.null(opt[['convert-rownames']])){
    CONVERT_ROWNAMES = FALSE
} else {
    CONVERT_ROWNAMES = opt[['convert-rownames']]
}

if (is.null(opt[['use-raw']])){
    USE_RAW = FALSE
} else {
    USE_RAW = opt[['use-raw']]
}

CLOBBER = ifelse(
    is.null(opt[['clobber']]), 
    FALSE,
    opt[['clobber']]
)


dirname = dirname(PATH)
basename = basename(PATH)
slogan = gsub(".h5ad", "", basename)

out_path = file.path(dirname, paste0(slogan, ".qs"))

if (!CLOBBER & file.exists(out_path)){
    stop(paste("Output file", out_path, "already exists. Use --clobber=T to overwrite."))
}

print(paste("Processing file:", PATH))
data <- read_h5ad(PATH)

print(paste("Loaded data from:", PATH))
print(paste("Number of cells:", nrow(data)))
print(paste("Number of features:", ncol(data)))

if (USE_RAW){
    print("Using raw data matrix(adata.raw.X)")
    data_x = data$raw$X
} else {
    print("Using data matrix (adata.x)")
    data_x = data$X
}

print(class(data_x))
print(dim(data_x))
print(data_x[100:115, 100:115])

counts_mat = t(as.matrix(data_x))

# set rownames to be gene names, and colnames to be cell names
colnames(counts_mat) = rownames(data$obs)
rownames(counts_mat) = rownames(data$var)

if (CONVERT_ROWNAMES){
    rownames(counts_mat) = convert_ensembl_to_symbol(rownames(counts_mat))
    print(paste("Converted rownames to gene symbols. Number of genes:", nrow(counts_mat)))
}
data <- CreateSeuratObject(
    counts = counts_mat, 
    meta.data = data$obs, 
    min.features = MIN_FEATURES, 
    min.cells = MIN_CELLS)
data[["RNA"]] <- as(object = data[["RNA"]], Class = "Assay")
print(paste("Created Seurat object with", ncol(data), "cells and", nrow(data), "features."))
print("Saving Seurat object...")
qsave(data, out_path)
print(paste("Saved Seurat object to:", out_path))