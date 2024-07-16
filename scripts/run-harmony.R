################# LIBRARIES #################

library(dplyr)
library(getopt)
library(harmony)
library(Matrix)
library(Seurat)
library(qs)

source("~/sc-online/clust.R")


################# PARSE ARGUMENTS #################

# TODO parametrize harmony params? 

spec <- matrix(c(
    'path', 'p', 1, "character",
    'ndims', 'nd', 1, 'integer',
    'group-by-vars', 'gbv', 1, 'character'
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)
PATH = opt[['path']]

NDIMS = ifelse(
    is.null(opt[['ndims']]), 
    20, 
    as.integer(opt[['ndims']]))

if (is.null(opt[['group-by-vars']])) {
    GROUP_BY_VARS = c("participant_id")
} else {
    GROUP_BY_VARS = strsplit(opt[['group-by-vars']], ",")[[1]]
}

RESOLUTIONS = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 1.25, 1.5)

sobj = qread(PATH)

sobj = (sobj 
    %>% RunHarmony(
        group.by.vars=GROUP_BY_VARS,
        dims.use=1:NDIMS,
        plot_convergence=F) 
    %>% FindNeighbors(dims=1:NDIMS,reduction="harmony")
)
    
for (res in RESOLUTIONS){
    sobj = sobj %>% FindClusters(resolution=res)}
    
sobj = sobj %>% RunUMAP(dims=1:20, reduction="harmony")

qsave(sobj, PATH)