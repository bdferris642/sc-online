library(getopt)
library(qs)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(zellkonverter)

spec <- matrix(c(
    'path', 'p', 1, "character",
    'assay', 'a', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
READ_PATH = opt[["path"]]
ASSAY = ifelse(
    is.null(opt[["assay"]]),
    "RNA",
    opt[["assay"]]
)

base_dir = dirname(READ_PATH)
slogan = gsub(".qs", "", basename(READ_PATH))
write_dir = file.path(base_dir, "h5seurat")
if (!dir.exists(write_dir)) {
    dir.create(write_dir)
}

# Load the Seurat object from the qs file
sobj = qread(READ_PATH)
DefaultAssay(sobj) = ASSAY
raw_counts = GetAssayData(sobj, assay = ASSAY, slot = "counts")

# create new sce with the raw counts and the sobj metadata 
sce = SingleCellExperiment(
    assays = list(counts = as.matrix(raw_counts)),
    colData = sobj@meta.data
)

cat("Writing .h5ad\n")
zellkonverter::writeH5AD(sce, file.path(write_dir, paste0(slogan, "_orig.h5ad")))

