### INPUTS
# `base-path`, str, absolute path to a dir containing a seurat object
# `read-basename`, str, basename of the Seurat object in `base-path`. Used to name the output file with the _markers suffix
# `ncores`, int, number of cores to use
# `cluster-col`, str, name of the column in the Seurat object that contains cluster assignments
# `logfc-threshold`, float, log fold change threshold for Seurat's FindAllMarkers
# `min-pct`, float, minimum percent of cells in a cluster that must express a gene for it to be considered a marker

### OUTPUTS
# creates a dataframe using Seurat's FindAllMarkers, saves that dataframe as a csv to `marker_path`

library(future)
library(qs)
library(getopt)
library(Seurat)

spec <- matrix(c(
  'base-path', 'b', 1, "character",
  'cluster-col', 'c', 1, "character",
  'ncores', 'n', 1, "integer",
  'read-basename', 'r', 1, "character",
  'log-fc-thresh', 'f', 1, "numeric",
  'min-pct', 'p', 1, "numeric",
  'only-pos', 'o', 1, "logical",
  'default-assay', 'a', 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

read_basename = opt[['read-basename']]
ncores = opt[['ncores']]
base_path = ifelse(
    is.null(opt[['base-path']]), 
    "/mnt/accessory/seq_data/pd_all", 
    opt[['base-path']])
cluster_col = ifelse(
    is.null(opt[['cluster-col']]), 
    "seurat_clusters", 
    opt[['cluster-col']])
logfc_threshold = ifelse(
    is.null(opt[['log-fc-thresh']]), 
    0.1, 
    opt[['log-fc-thresh']])
min_pct = ifelse(
    is.null(opt[['min-pct']]), 
    0.01, 
    opt[['min-pct']]
  )
ONLY_POS = ifelse(
    is.null(opt[['only-pos']]), 
    TRUE, 
    opt[['only-pos']]
  )
DEFAULT_ASSAY = ifelse(
    is.null(opt[['default-assay']]), 
    "RNA", 
    opt[['default-assay']]
  )

MARKER_SUBDIR = "markers"

plan(strategy = "multicore", workers = ncores)
options(future.globals.maxSize = 5000 * 1024^2) # Increase limit to 5000 MiB


# remove .qs suffix from read_basename. Will use this slogan in writing the output.
slogan = gsub("\\.qs$", "", read_basename)
read_path = file.path(base_path, read_basename)
dir.create(file.path(base_path, MARKER_SUBDIR), showWarnings = FALSE)
marker_path = file.path(base_path, MARKER_SUBDIR, paste0(slogan, "_markers_", cluster_col, ".csv"))
print(paste('Run Slogan =', slogan))
print(paste('Reading Seurat object at', read_path, '...'))
s_obj = qread(read_path)

DefaultAssay(s_obj) = DEFAULT_ASSAY

Idents(s_obj) = s_obj[[cluster_col]]

print(paste('Collecting markers and saving to', marker_path, '...'))
# Find all markers; adjust parameters as needed
markers = FindAllMarkers(
  s_obj, 
  logfc.threshold=logfc_threshold, min.pct=min_pct,
  only.pos=ONLY_POS, verbose=TRUE, return.thresh=0.05)

write.csv(markers, marker_path)