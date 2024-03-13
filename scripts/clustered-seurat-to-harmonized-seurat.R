# Takes in a Seurat object that has already been clustered 
# outputs a harmonized version of each subcluster (saving with the same slogan as the read name)

### INPUTS
# `base-path`, str, absolute path to a dir containing a seurat object
# `read-basename`, str, basename of the Seurat object in `base-path`
# `cluster_col`, str, name of the column in the Seurat object that contains cluster assignments
# `split-col`, str, name of the column in the Seurat object that contains the split you want to harmonize (e.g. donor_id, participant_id, etc.)
# `ncores`, int, number of cores to use
# `nfeatures`, int, number of HVGs to calculate in each split
# `npcs`, int, number of PCs to use in Harmony

### OUTPUTS
# qsaves a seurat object for each subcluster in `cluster_col` 
# with the same slogan as the read name and the subcluster name appended to the end
# these seurat objects are harmonized by the split col and contain only the HVGs that are common to at least half of the splits


################################# LIBS #################################

library(dplyr)
library(future)
library(harmony)
library(qs)
library(getopt)
library(Seurat)

source("~/sc-online/utils.R")
source("~/sc-online/getData.R")

################################# PARSE ARGUMENTS #################################

spec <- matrix(c(
    'base-path', 'b', 1, "character",
    'read-basename', 'r', 1, "character",
    'cluster-col', 'c', 1, "character",
    'split-col', 's', 1, "character",
    'ncores', 'n', 1, "integer",
    'nfeatures', 'f', 1, "integer",
    'npcs', 'p', 1, "integer"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

read_basename = opt[['read-basename']]
ncores = ifelse(
    is.null(opt[['ncores']]), 
    1, 
    opt[['ncores']])
base_path = ifelse(
    is.null(opt[['base-path']]), 
    "/mnt/accessory/seq_data/pd_all", 
    opt[['base-path']])
cluster_col = ifelse(
    is.null(opt[['cluster-col']]), 
    "broad_cell_class", 
    opt[['cluster-col']])
split_col = ifelse(
    is.null(opt[['split-col']]), 
    "participant_id", 
    opt[['split-col']])
nfeatures = ifelse(
    is.null(opt[['nfeatures']]), 
    5000, 
    opt[['nfeatures']])
npcs = ifelse(
    is.null(opt[['npcs']]), 
    50, 
    opt[['npcs']])

################################# MAIN #################################

#plan(strategy = "multicore", workers = ncores)
options(future.globals.maxSize = 10000 * 1024^2) # Increase limit to 10,000 MiB

# remove .qs suffix from read_basename. Will use this slogan in writing the output.
slogan = gsub("\\.qs$", "", read_basename)
read_path = file.path(base_path, read_basename)
print(paste('Run Slogan =', slogan))
print(paste('Reading Seurat object at', read_path, '...'))
s_obj = qread(read_path)

cell_classes = unique(s_obj@meta.data[[cluster_col]])
print(paste0("Groups in cluster column ", cluster_col, ": "))
print(cell_classes)


for (cell_class in cell_classes) {
    write_path = file.path(base_path, paste0(slogan, "_", cell_class, ".qs"))
    print(paste(cell_class, write_path))
    s_obj_cell_class = s_obj[, s_obj[[cluster_col]] == cell_class]
    print(paste('Number of cells in subcluster', cell_class, ':', ncol(s_obj_cell_class)))
  
    ################################# FIND HVGs #################################
    print(paste0("finding HVGs in at least half of ", split_col, "'s"))
    
    subsets = Seurat::SplitObject(s_obj_cell_class, split.by = split_col)
    # only keep subsets with at least 2 cells
    subsets = Filter(function(x) ncol(x) > 1, subsets)
    n_subsets = length(subsets)
    n_subsets_hvg = floor(n_subsets/2)

    hvgs_by_subset = list()
    for (i in names(subsets)) {
        subset = subsets[[i]]
        print(paste0("finding HVGs in subset ", i, " of ", n_subsets))
        subset = FindVariableFeatures(subset, selection.method = "vst", nfeatures = nfeatures)
        hvgs_by_subset[[i]] = subset@assays$RNA@var.features
    }

    # Get list of genes occuring in some number of donors within a sort.
    hvgs = getCommonStrings(hvgs_by_subset, n_subsets_hvg)
    print(paste("Number of HVGs in common across", n_subsets_hvg, "--", length(hvgs)))
    # only keep the donors we used to grab those HVGs
    s_obj_cell_class = s_obj_cell_class[, s_obj_cell_class@meta.data[[split_col]] %in% names(subsets)]
    # set hvgs
    s_obj_cell_class@assays$RNA@var.features = hvgs

    ################################# PCA & HARMONY #################################
    print(paste('Harmonizing subcluster', cell_class, '...'))
    s_obj_cell_class = (s_obj_cell_class
        %>% RunPCA(features = hvgs, npcs = npcs)
        %>% RunHarmony(
            group.by.vars=split_col,
            dims.use = 1:npcs,
            early_stop = F,
            max_iter=25,
            plot_convergence = TRUE
        )
    )

    print(paste('Writing harmonized subcluster', cell_class, 'to', write_path, '...'))
    qsave(s_obj_cell_class, write_path)
}
