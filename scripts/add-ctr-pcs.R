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


suppressMessages(suppressWarnings({
    library(dplyr)
    library(getopt)
    library(Matrix)
    library(qs)
    library(Seurat)
    source("../utils.R")
    g = glue::glue
}))

spec <- matrix(c(
    'ctr-path', 'p', 1, "character", 
    'sobj-path', 's', 1, "character",
    'add-ctr-pcs', 'c', 0, "logical",
    'assay', 'a', 2, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
CTR_PATH = opt[['ctr-path']]
SOBJ_PATH = opt[['sobj-path']]

if (is.null(opt[['add-ctr-pcs']])) {
    ADD_CTR_PCS = FALSE
} else {
    ADD_CTR_PCS = opt[['add-ctr-pcs']]
}

if (is.null(opt[['assay']])) {
    pca_assay = "RNA"
} else {
    pca_assay = opt[['assay']]
}

ctr = qread(CTR_PATH)
gene_loadings = ctr[["pca"]]@feature.loadings
pca_genes = as.character(rownames(gene_loadings))

message("Processing: ", SOBJ_PATH)
sobj = qread(SOBJ_PATH)
if (! pca_assay %in% names(sobj@assays)) {
    stop(g("Assay '{pca_assay}' not found in target object! Exiting."))
}

# if no `data` layer in sobj, normalize  
remove_data_layer = FALSE
if (! "data" %in% names(sobj[[pca_assay]]@layers)) {
    message("No 'data' layer found in target object; normalizing data...")
    sobj = NormalizeData(sobj, assay = pca_assay)
    remove_data_layer = TRUE
}

## loadings & feature set
gene_loadings = ctr[["pca"]]@feature.loadings   # genes x PCs
pca_genes = rownames(gene_loadings)

## --- reference (ctr) matrix on the SAME scale used for PCA ---
ctr_mat = GetAssayData(ctr, assay = pca_assay, layer = "data")[pca_genes, , drop = FALSE]

## reference means/SDs (across cells)
ctr_means = Matrix::rowMeans(ctr_mat)
ctr_sds   = apply(ctr_mat, 1, sd)
ctr_sds[ctr_sds == 0 | !is.finite(ctr_sds)] = 1  # guard

## --- target object (sobj) on the SAME assay/layer ---
sobj_mat_all = GetAssayData(sobj, assay = pca_assay, layer = "data")

## restrict to shared PCA genes
these_genes = intersect(pca_genes, as.character(rownames(sobj_mat_all)))

message(g("{length(these_genes)} / {length(pca_genes)} in CTR object PCA shared with target object."))

## align everything to the same gene order
these_genes = these_genes[order(match(these_genes, pca_genes))]
sobj_mat    = sobj_mat_all[these_genes, , drop = FALSE]
ref_means   = ctr_means[these_genes]
ref_sds     = ctr_sds[these_genes]
loadings_these = gene_loadings[these_genes, , drop = FALSE]

## center/scale target cells using reference stats (genes as columns for scale())
sobj_scaled = scale(t(as.matrix(sobj_mat)), center = ref_means, scale = ref_sds)  # cells x genes

## project: (cells x genes) %*% (genes x PCs) => (cells x PCs)
sobj_pca = as.data.frame(sobj_scaled %*% as.matrix(loadings_these))

## clean old PCs and join
sobj@meta.data = sobj@meta.data[, !grepl("^PC_", colnames(sobj@meta.data)), drop = FALSE]
sobj = join_df_to_sobj_metadata(sobj, sobj_pca, "row.names")

print("Saving Case with PCs")
if (remove_data_layer){
    print("Removing 'data' layer added for normalization...")
    sobj[[pca_assay]]@layers[["data"]] = NULL
}
qsave(sobj, SOBJ_PATH)

print("Done adding PCs to cases!")

if (ADD_CTR_PCS){
    print("Adding PCs to CTR...")
    ctr_mat_target = ctr_mat[pca_genes, , drop = FALSE]     # already on correct assay/layer
    ctr_scaled = scale(t(as.matrix(ctr_mat_target)), center = ctr_means, scale = ctr_sds)
    ctr_pca = as.data.frame(ctr_scaled %*% as.matrix(gene_loadings))

    ctr@meta.data = ctr@meta.data[, !grepl("^PC_", colnames(ctr@meta.data)), drop = FALSE]
    ctr = join_df_to_sobj_metadata(ctr, ctr_pca, "row.names")

    print("Saving CTR with PCs")
    qsave(ctr, CTR_PATH)
}

print("DONE!")
