###################### IMPORTS ######################

library(getopt)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(tidyr)

source("~/code/sconline_code.R") # <-- TODO inherit from pseudocells.R
source("~/sc-online/utils.R")

###################### PARSE ARGUMENTS ######################

spec <- matrix(c(
    'base-path', 'bp', 1, "character",
    'basename', 'bn', 1, 'character',
    'pseudocell-size', 'ps', 1, 'numeric',
    'cluster-col', 'cc', 1, 'character',
    'grouping-col', 'gc', 1, 'character',
    'embedding', 'e', 1, 'character',
    'suffix', 's', 1, 'character'
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

PSEUDOCELL_SIZE = ifelse(
    is.null(opt[['pseudocell-size']]), 
    30, 
    as.numeric(opt[['pseudocell-size']]))
CLUSTER_COL = ifelse(
    is.null(opt[['cluster-col']]), 
    "cell_class", 
    opt[['cluster-col']])
GROUPING_COL = ifelse(
    is.null(opt[['grouping-col']]), 
    "participant_id", 
    opt[['grouping-col']]
)
EMBD = ifelse(
    is.null(opt[['embedding']]), 
    "pca", 
    opt[['embedding']]
)
ADDITIONAL_SUFFIX = ifelse(
    is.null(opt[['suffix']]), 
    "", 
    paste0('__', opt[['suffix']])
)

BASE_PATH = opt[['base-path']]
BASENAME = opt[['basename']]
if (is.null(BASE_PATH) || is.null(BASENAME)) {
    stop("Please provide both a base-path and basename")
}

# Hard Code For Now
NCORES = 10
N_PCS = 30
MIN_SIZE_LIMIT = 10

###################### MAIN ######################
suffix = paste0('_pseudocells__split_by_',  CLUSTER_COL, '__grouped_by_', GROUPING_COL, '__mean_size_', PSEUDOCELL_SIZE, ADDITIONAL_SUFFIX, '.qs')

read_path = file.path(BASE_PATH, BASENAME)
write_path = gsub(".qs", suffix, read_path)

print(read_path) 
print(write_path)

# load object
s_obj = qread(read_path)
print(paste("Seurat Object Dimensions:", unlist(dim(s_obj))))

# Extracting the cell embedding data
embedding_data = Embeddings(s_obj, EMBD)[,1:N_PCS]
print(paste("Embedding Data Dimensions:", unlist(dim(embedding_data))))

data_sce = as.SingleCellExperiment(s_obj)

pseudocells=suppressWarnings(
    .sconline.PseudobulkGeneration(argList=NULL, 
    parsing.col.names = c(CLUSTER_COL, GROUPING_COL),
    pseudocell.size=PSEUDOCELL_SIZE,
    inputExpData=data_sce,
    min_size_limit=MIN_SIZE_LIMIT,
    inputPhenoData=as.data.frame(colData(data_sce)),
    inputEmbedding=embedding_data,
    nPCs=N_PCS,
    ncores=NCORES,
    rand_pseudobulk_mod=F,
    organism="Human")
)

print("Number of Pseudocells:", length(pseudocells))

#Transforming the nGene and nUMI
pseudocells$QC_Gene_total_log=log2(pseudocells$QC_Gene_total_count)
pseudocells$QC_Gene_unique_log=log2(pseudocells$QC_Gene_unique_count)

pseudocells_list=.mySplitObject(pseudocells, CLUSTER_COL)

#Writing the pseudocells
qsave(pseudocells_list, write_path)
