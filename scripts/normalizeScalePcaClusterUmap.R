################# LIBRARIES #################

library(dplyr)
library(getopt)
library(Matrix)
library(Seurat)
library(qs)
library(harmony)
library(dplyr)

source("~/sc-online/clust.R")


################# PARSE ARGUMENTS #################

spec <- matrix(c(
    'path', 'p', 1, "character",
    'scaling-subset-col', 'ssc', 1, 'character',
    'var-feature-subset-col', 'vsc', 1, 'character',
    'regression-vars', 'rv', 1, 'character',
    'run-harmony', 'h', 1, 'logical',
    'harmony-group-by-vars', 'hg', 1, 'character'
), byrow = TRUE, ncol = 4)

RESOLUTIONS=c(0.2, 0.5)

opt <- getopt(spec)
PATH = opt[['path']]

if (is.null(opt[['scaling-subset-col']])) {
    SCALING_SUBSET_COL = NULL
} else {
    SCALING_SUBSET_COL = opt[['scaling-subset-col']]
}

if (is.null(opt[['var-feature-subset-col']])) {
    VAR_FEATURE_SUBSET_COL = NULL
} else {
    VAR_FEATURE_SUBSET_COL = opt[['var-feature-subset-col']]
}

if (is.null(opt[["regression-vars"]])) {
    REGRESSION_VARS = NULL
} else {
    REGRESSION_VARS = strsplit(opt[["regression-vars"]], ",")[[1]]
}

if (is.null(opt[["run-harmony"]])) {
    RUN_HARMONY = FALSE
} else {
    RUN_HARMONY = opt[["run-harmony"]]
}

if (is.null(opt[["harmony-group-by-vars"]])) {
    HARMONY_GROUP_BY_VARS = NULL
} else {
    HARMONY_GROUP_BY_VARS = strsplit(opt[["harmony-group-by-vars"]], ",")[[1]]
}

sobj = qread(PATH)

sobj = normalizeScalePcaClusterUmap(
    sobj,
    var_feature_subset_col=VAR_FEATURE_SUBSET_COL,
    scaling_subset_col=SCALING_SUBSET_COL, 
    n_hvgs_orig=2500, 
    n_dims_use=20,
    regression_vars=REGRESSION_VARS,
    resolutions=RESOLUTIONS)

if (RUN_HARMONY) {
    sobj = (sobj 
        %>% RunHarmony(group.by.vars=HARMONY_GROUP_BY_VARS, dims.use=1:20)
        %>% FindNeighbors(dims=1:20, reduction="harmony")
    )
    for (res in RESOLUTIONS) {
        sobj = sobj %>% FindClusters(resolution=res)
    }
    sobj = sobj %>% RunUMAP(dims=1:20, reduction="harmony")
}

qsave(sobj, PATH)