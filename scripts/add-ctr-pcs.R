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
    'sobj-path-list', 'l', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
CTR_PATH = opt[['ctr-path']]
SOBJ_PATH_LIST = opt[['sobj-path-list']]

# read the lines of the case path list file
qs_paths = readLines(SOBJ_PATH_LIST)


ctr = qread(CTR_PATH)
gene_loadings = ctr[["pca"]]@feature.loadings
pca_genes = rownames(gene_loadings)
ctr_mat = GetAssayData(ctr, assay="RNA", slot="counts")[pca_genes, ]
ctr_means = Matrix::rowMeans(ctr_mat)
ctr_sds = apply(ctr_mat, 1, sd)

for (qs_path in qs_paths) {
    message("Processing: ", qs_path)
    sobj = qread(qs_path)

    if (! all(pca_genes %in% rownames(sobj))) {
        warning("Not all PCA genes found in combined object: ", qs_path, "!")
    }

    these_genes = pca_genes[pca_genes %in% rownames(sobj)]
    ctr_mat_these_genes = ctr_mat[these_genes, ]
    ctr_means_these_genes = Matrix::rowMeans(ctr_mat_these_genes)
    ctr_sds_these_genes = apply(ctr_mat_these_genes, 1, sd)


    sobj_counts = GetAssayData(sobj, assay="RNA", slot="counts")[these_genes, ]
    sobj_scaled = scale(t(as.matrix(sobj_counts)), center=ctr_means_these_genes, scale=ctr_sds_these_genes)
    sobj_pca = as.data.frame(sobj_scaled %*% gene_loadings)

    sobj = join_df_to_sobj_metadata(
        sobj, 
        sobj_pca, 
        "row.names"
    )

    qsave(sobj, qs_path)
}
