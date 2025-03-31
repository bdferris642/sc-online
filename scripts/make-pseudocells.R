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
    'path', 'p', 1, "character",
    'pseudocell-size', 'ps', 1, 'numeric',
    'min-size', 'ms', 1, 'numeric',
    'cluster-col', 'cc', 1, 'character',
    'grouping-cols', 'gc', 1, 'character',
    'embedding', 'e', 1, 'character',
    'rand', 'r', 1, 'logical',
    'suffix', 's', 1, 'character',
    'default_assay', 'da', 1, 'character',
    "filter-str", "fs", 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

PSEUDOCELL_SIZE = ifelse(
    is.null(opt[['pseudocell-size']]), 
    30, 
    as.numeric(opt[['pseudocell-size']]))
MIN_SIZE_LIMIT = ifelse(
    is.null(opt[['min-size']]), 
    10, 
    as.numeric(opt[['min-size']]))

CLUSTER_COL = ifelse(
    is.null(opt[['cluster-col']]), 
    "cell_class", 
    opt[['cluster-col']])
if (is.null(opt[['grouping-cols']])){
    GROUPING_COLS = c("participant_id")
    GROUPING_COL_STR = "participant_id"
} else {
    GROUPING_COL_STR = sub(',', "_and_", opt[['grouping-cols']])
    GROUPING_COLS = strsplit(opt[['grouping-cols']], ",")[[1]]
}

RAND = ifelse(
    is.null(opt[['rand']]), 
    TRUE, 
    opt[['rand']]
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

DEFAULT_ASSAY = ifelse(
    is.null(opt[['default_assay']]), 
    "RNA", 
    opt[['default_assay']]
)

if (is.null(opt[['filter-str']])) {
    FILTER_STR = NULL
} else {
    FILTER_STR = opt[['filter-str']]
    clean_filter_str <- FILTER_STR
    print(clean_filter_str)
    clean_filter_str <- gsub("|", "or", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    clean_filter_str <- gsub("!", "not", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    clean_filter_str <- gsub("&", "and", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    clean_filter_str <- gsub(">", "gt", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    clean_filter_str <- gsub("<", "lt", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    clean_filter_str <- gsub("=", "", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    clean_filter_str <- gsub(" ", "_", clean_filter_str, fixed = TRUE)
    print(clean_filter_str)
    ADDITIONAL_SUFFIX = paste0(ADDITIONAL_SUFFIX, "__", clean_filter_str)
}


PATH = opt[['path']]
if (is.null(PATH)){
    stop("Path is required")
}
BASE_PATH = dirname(PATH)
BASENAME = basename(PATH)

# Hard Code For Now
NCORES = 28
N_PCS = 20

print(paste("Path:", PATH))
print(paste("Pseudocell Size:", PSEUDOCELL_SIZE))
print(paste("Min Size Limit:", MIN_SIZE_LIMIT))
print(paste("Cluster Column:", CLUSTER_COL))
print(paste("Grouping Columns:", GROUPING_COLS))
print(paste("Embedding:", EMBD))
print(paste("Randomize:", RAND))
print(paste("Additional Suffix:", ADDITIONAL_SUFFIX))
print(paste("Default Assay:", DEFAULT_ASSAY))
print(paste("Filter String:", FILTER_STR))


###################### MAIN ######################
if (toupper(CLUSTER_COL) %in% c("NULL", "NONE", "NA")){
    parsing_cols = GROUPING_COLS
    suffix = paste0('_pseudocells__unsplit__grouped_by_', GROUPING_COL_STR, '__mean_size_', PSEUDOCELL_SIZE, "__rand_", RAND, ADDITIONAL_SUFFIX, '.qs')
} else {
    parsing_cols = c(CLUSTER_COL, GROUPING_COLS)
    suffix = paste0('_pseudocells__split_by_',  CLUSTER_COL, '__grouped_by_', GROUPING_COL_STR, '__mean_size_', PSEUDOCELL_SIZE, "__rand_", RAND, ADDITIONAL_SUFFIX, '.qs')
}


read_path = file.path(BASE_PATH, BASENAME)
write_path = gsub(".qs", suffix, read_path)

print(read_path) 
print(write_path)

# load object
s_obj = qread(read_path)

if (!is.null(FILTER_STR)) {
    logical_vector = with(s_obj@meta.data, eval(parse(text = FILTER_STR)))
    s_obj = s_obj[, logical_vector]
}

DefaultAssay(s_obj) = DEFAULT_ASSAY
print(paste("Seurat Object Dimensions:", unlist(dim(s_obj))))

# Extracting the cell embedding data
embedding_data = Embeddings(s_obj, EMBD)[,1:N_PCS]
print(paste("Embedding Data Dimensions:", unlist(dim(embedding_data))))

data_sce = as.SingleCellExperiment(s_obj)

pseudocells=suppressWarnings(
    .sconline.PseudobulkGeneration(argList=NULL, 
    parsing.col.names = parsing_cols,
    pseudocell.size=PSEUDOCELL_SIZE,
    inputExpData=data_sce,
    min_size_limit=MIN_SIZE_LIMIT,
    inputPhenoData=as.data.frame(colData(data_sce)),
    inputEmbedding=embedding_data,
    nPCs=N_PCS,
    ncores=NCORES,
    rand_pseudobulk_mod=RAND,
    organism="Human")
)

#Transforming the nGene and nUMI
pseudocells$QC_Gene_total_log=log2(pseudocells$QC_Gene_total_count)
pseudocells$QC_Gene_unique_log=log2(pseudocells$QC_Gene_unique_count)

if (toupper(CLUSTER_COL) %in% c("NULL", "NONE", "NA")){
    pseudocells_list=list(pseudocells)
} else {
    pseudocells_list=.mySplitObject(pseudocells, CLUSTER_COL)
}

#Writing the pseudocells
qsave(pseudocells_list, write_path)
