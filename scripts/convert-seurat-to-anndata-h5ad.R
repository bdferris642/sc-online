suppressWarnings(suppressMessages({
    library(Seurat)
    library(SeuratDisk)
    library(anndata)
    library(getopt)
    library(Matrix)
    library(qs)
}))

# Parse command line arguments
spec = matrix(c(
    'path', 'p', 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
PATH = opt[['path']]


dirname = dirname(PATH)
basename = basename(PATH)
slogan = gsub(".qs", "", basename)
h5seurat_fname = paste0(dirname, "/", slogan, ".h5Seurat")

library(Seurat)
library(SeuratDisk)

# Load Seurat object
seurat_obj <- qread(PATH)
DefaultAssay(seurat_obj) <- "RNA"

# Keep only the RNA assay
# if seurat object version is 5 or later, counts is accessed with $, 4 or earlier with @

seurat_obj[["RNA"]] = CreateAssayObject(counts = GetAssayData(seurat_obj, slot="counts", assay="RNA"))

# Remove all other assays
seurat_obj@assays <- list(RNA = seurat_obj[["RNA"]])

# Ensure metadata is properly aligned
if (!all(rownames(seurat_obj@meta.data) == colnames(seurat_obj))) {
    stop("Cell names in meta.data do not match counts matrix! Check alignment.")
}

# Save the cleaned Seurat object
SaveH5Seurat(seurat_obj, filename = h5seurat_fname, overwrite = TRUE)

# Convert to H5AD
Convert(h5seurat_fname, dest = "h5ad", overwrite = TRUE)
