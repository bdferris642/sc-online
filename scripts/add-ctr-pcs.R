print("**************** LOADING LIBRARIES ****************")
# Detect script path when running via Rscript
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir <- dirname(normalizePath(script_path))
  setwd(script_dir)
  message("Working directory set to: ", script_dir)
} else {
  stop("Cannot determine script path. Are you running via Rscript?")
}


suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(source("../utils.R")))

spec <- matrix(c(
    'ctr-path', 'p', 1, "character", 
    'sobj-path', 's', 1, "character",
    'add-ctr-pcs', 'c', 0, "logical"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
CTR_PATH = opt[['ctr-path']]
SOBJ_PATH = opt[['sobj-path']]

if (is.null(opt[['add-ctr-pcs']])) {
    ADD_CTR_PCS = FALSE
} else {
    ADD_CTR_PCS = opt[['add-ctr-pcs']]
}


ctr = qread(CTR_PATH)
gene_loadings = ctr[["pca"]]@feature.loadings
pca_genes = rownames(gene_loadings)
ctr_mat = GetAssayData(ctr, assay="RNA", slot="counts")[pca_genes, ]
ctr_means = Matrix::rowMeans(ctr_mat)
ctr_sds = apply(ctr_mat, 1, sd)


message("Processing: ", SOBJ_PATH)
sobj = qread(SOBJ_PATH)

if (! all(pca_genes %in% rownames(sobj))) {
    warning("Not all PCA genes found in combined object: ", SOBJ_PATH, "!")
}

these_genes = pca_genes[pca_genes %in% rownames(sobj)]
ctr_mat_these_genes = ctr_mat[these_genes, ]
ctr_means_these_genes = Matrix::rowMeans(ctr_mat_these_genes)
ctr_sds_these_genes = apply(ctr_mat_these_genes, 1, sd)

sobj_counts = GetAssayData(sobj, assay="RNA", slot="counts")[these_genes, ]
sobj_scaled = scale(t(as.matrix(sobj_counts)), center=ctr_means_these_genes, scale=ctr_sds_these_genes)
sobj_pca = as.data.frame(sobj_scaled %*% gene_loadings)

# remove columns starting with PC_ if they exist
sobj@meta.data = sobj@meta.data[, !grepl("^PC_", colnames(sobj@meta.data)), ]

sobj = join_df_to_sobj_metadata(
    sobj, 
    sobj_pca, 
    "row.names"
)
print("Saving Case with PCs")
qsave(sobj, SOBJ_PATH)

print("Done adding PCs to cases!")

if (ADD_CTR_PCS){
    # also add to the ctr itself
    print("Adding PCs to CTR itself...")
    ctr_counts = GetAssayData(ctr, assay="RNA", slot="counts")[pca_genes, ]
    ctr_scaled = scale(t(as.matrix(ctr_counts)), center=ctr_means, scale=ctr_sds)
    ctr_pca = as.data.frame(ctr_scaled %*% gene_loadings)

    # remove columns starting with PC_ if they exist
    ctr@meta.data = ctr@meta.data[, !grepl("^PC_", colnames(ctr@meta.data)), ]

    ctr = join_df_to_sobj_metadata(
        ctr, 
        ctr_pca, 
        "row.names"
    )
    print("Saving CTR with PCs")
    qsave(ctr, CTR_PATH)
}

print("DONE!")
