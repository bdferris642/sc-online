suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SeuratDisk)))
suppressMessages(suppressWarnings(library(reticulate)))
suppressMessages(suppressWarnings(library(anndata)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(library(here)))
suppressMessages(suppressWarnings(source(here("utils.R"))))

use_python("/home/ferris/anaconda3/envs/r-4.3.0/bin/python", required = TRUE) # OJO have to change this path to conda path

# Parse command line arguments
spec = matrix(c(
    'path', 'p', 1, "character",
    'min-features', 'f', 1, "numeric",
    'min-cells', 'c', 1, "numeric"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
PATH = opt[['path']]
MIN_FEATURES = ifelse(
    is.null(opt[['min-features']]), 
    500,
    opt[['min-features']]
)
MIN_CELLS = ifelse(
    is.null(opt[['min-cells']]), 
    30,
    opt[['min-cells']]
)


dirname = dirname(PATH)
basename = basename(PATH)
slogan = gsub(".h5ad", "", basename)

print(paste("Processing file:", PATH))
data <- read_h5ad(PATH)

print(paste("Loaded data from:", PATH))
print(paste("Number of cells:", nrow(data)))
print(paste("Number of features:", ncol(data)))

data_x = data$X
print(class(data_x))
print(dim(data_x))
print(data_x[100:110, 100:110])

counts_mat = t(as.matrix(data_x))
rownames(counts_mat) = convert_ensembl_to_symbol(rownames(counts_mat))
print(paste("Converted rownames to gene symbols. Number of genes:", nrow(counts_mat)))
data <- CreateSeuratObject(
    counts = counts_mat, 
    meta.data = data$obs, 
    min.features = MIN_FEATURES, 
    min.cells = MIN_CELLS)
data[["RNA"]] <- as(object = data[["RNA"]], Class = "Assay")
print(paste("Created Seurat object with", ncol(data), "cells and", nrow(data), "features."))
print("Saving Seurat object...")
qsave(data, file.path(dirname, paste0(slogan, ".qs")))
print(paste("Saved Seurat object to:", file.path(dirname, paste0(slogan, ".qs"))))