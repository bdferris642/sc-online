# INPUTS: A path to a seurat object and some pseudobulking parameters
# OUTPUTS: A list of pseudobulks (Seurat or SingleCellExperiment objects) saved to disk as .qs files

# The number of pseudobulks in the list will correspond to the number of distinct values of CLUSTER_COL in the Seurat object
# Each Pseudobulk in the list will have at least as many columns as there are distinct groupings of GROUPING_COLS in the input data. 

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


suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

suppressMessages(suppressWarnings(source("../de.R")))
suppressMessages(suppressWarnings(source("../utils.R")))
suppressMessages(suppressWarnings(source("../getData.R")))

print("**************** PARSING ARGUMENTS ****************")
spec <- matrix(c(
    'path', 'p', 1, "character", # seurat qs 
    'output-type', 'o', 1, 'character', # sce or seurat (default sce)
    'contrast-col', 'cc', 1, 'character', # optional
    'cluster-col', 'cl', 1, 'character', # optional, if provided, will pseudobulk within each cluster and output a list of pseudobulks
    'sample-col', 'sc', 1, 'character', # required, column in metadata that indicates sample identity
    'grouping-cols', 'g', 1, 'character', # required, comma-separated list of columns in metadata to group by when pseudobulking
    'cols-to-mean', 'c', 1, 'character', # optional, comma-separated list of columns in metadata to average when pseudobulking
    'cols-to-weighted-avg', 'w', 1, 'character', # optional, comma-separated list of columns in metadata to average weighted by number of UMIs when pseudobulking
    'cols-to-median', 'cm', 1, 'character', # optional, comma-separated list of columns in metadata to take median when pseudobulking
    'only-clusters', 'oc', 1, 'character', # optional, comma-separated list of clusters to pseudobulk, if not provided will pseudobulk all clusters
    'filter-samples', 'fs', 1, 'character', # optional, comma-separated list of samples to filter out before pseudobulking
    'assay', 'a', 1, 'character', # optional, default RNA
    'min-n-cells', 'mnc', 1, 'numeric', # optional, default 10 (minimum number of cells per pseudobulk)
    'min-counts-gene', 'mcg', 1, 'numeric', # optional, (minimum number of counts in a gene for it to be included in a pseudobulk)
    'min-frac-gene', 'mfg', 1, 'numeric', # optional, (minimum fraction of pseudobulks that must have a count for a gene to be included)
    'suffix', 's', 1, 'character',
    'sva-cols', 'sva', 1, 'character',
    'sva-ctr-cols', 'svac', 1, 'character',
    'n-svs', 'nsv', 1, 'integer',
    'drop-na', 'dn', 1, 'logical'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]

GROUPING_COLS = strsplit(opt[['grouping-cols']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
SAMPLE_COL = opt[['sample-col']]

cat(paste("Reading Seurat object from: ", PATH, "\n"))
sobj = load_obj(PATH)
DefaultAssay(sobj) = "RNA"

# OJO: fixing a naming issue in a simulation run
sobj@meta.data$Row.names = NULL
sobj@meta.data$Row.names = colnames(sobj)

if (is.null(opt[['filter-samples']])){
    FILTER_SAMPLES = NULL
} else {
    FILTER_SAMPLES = strsplit(opt[['filter-samples']], ",")[[1]]
    message(paste("\nFiltering samples: ", paste(FILTER_SAMPLES, collapse=", ")))
    sobj = sobj[,!sobj@meta.data[[SAMPLE_COL]] %in% FILTER_SAMPLES] 
}

cat(paste("Seurat Object Loaded with", nrow(sobj), "genes and", ncol(sobj), "cells\n"))
md = sobj@meta.data

md_cols = colnames(md)

if (! SAMPLE_COL %in% colnames(md)){
    stop(paste("Error: Sample column not found in metadata: ", SAMPLE_COL))
}
md[[SAMPLE_COL]] = paste0("sample_", md[[SAMPLE_COL]])
md[[SAMPLE_COL]] = gsub("-", "_", md[[SAMPLE_COL]])

if (is.null(opt[['cluster-col']])){
    CLUSTER_COL = NULL
} else {
    CLUSTER_COL = opt[['cluster-col']]
    if (! CLUSTER_COL %in% colnames(md)){
        stop(paste("Error: Cluster column not found in metadata: ", CLUSTER_COL))
    }
}

if (! is.null(CLUSTER_COL)){
    md[[CLUSTER_COL]] = gsub(" ", "_", md[[CLUSTER_COL]])

    if( !(CLUSTER_COL %in% GROUPING_COLS)){
        print(paste("Adding cluster column to grouping columns: ", CLUSTER_COL))
        GROUPING_COLS = c(GROUPING_COLS, CLUSTER_COL)
    }
}

if (! is.null(CONTRAST_COL)){
    md[[CONTRAST_COL]] = factor(md[[CONTRAST_COL]], levels = sort(unique(md[[CONTRAST_COL]])))
}

if ("age" %in% colnames(md)){
    md$age = as.numeric(md$age)
}
sobj@meta.data = md

OUTPUT_TYPE = if(is.null(opt[['output-type']])){
    "sce"
} else {
    opt[['output-type']]
}

ONLY_CLUSTERS = if(is.null(opt[['only-clusters']])){
    ONLY_CLUSTERS = NULL
} else {
    ONLY_CLUSTERS = strsplit(opt[['only-clusters']], ",")[[1]]
}

if (is.null(opt[['assay']])){
    ASSAY = "RNA"
} else {
    ASSAY = opt[['assay']]
}
if(is.null(opt[['min-n-cells']])){
    MIN_N_CELLS = 10
} else {
    MIN_N_CELLS = opt[['min-n-cells']]
}
if(is.null(opt[['min-counts-gene']])){
    MIN_COUNTS_GENE = 0
} else {
    MIN_COUNTS_GENE = opt[['min-counts-gene']]
}
if(is.null(opt[['min-frac-gene']])){
    MIN_FRAC_GENE = 0
} else {
    MIN_FRAC_GENE = opt[['min-frac-gene']]
}

if(is.null(opt[['suffix']])){
    suffix = ""
} else {
    suffix = opt[['suffix']]
}

if (is.null(opt[['cols-to-mean']])){
    COLS_TO_MEAN = NULL
} else {
    COLS_TO_MEAN = strsplit(opt[['cols-to-mean']], ",")[[1]]
}

if (is.null(opt[['cols-to-weighted-avg']])){
    COLS_TO_WEIGHTED_AVG = NULL
} else {
    COLS_TO_WEIGHTED_AVG = strsplit(opt[['cols-to-weighted-avg']], ",")[[1]]
}

if (is.null(opt[['cols-to-median']])){
    COLS_TO_MEDIAN = NULL
} else {
    COLS_TO_MEDIAN = strsplit(opt[['cols-to-median']], ",")[[1]]
}

slogan = gsub(".qs", "", basename(PATH))
slogan = gsub(".rds", "", slogan)
output_dir = file.path(dirname(PATH), "pseudobulk", slogan)
if (! dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}

if (! is.null(CLUSTER_COL)){
    clusters = unique(sobj@meta.data[[CLUSTER_COL]])
    if (! is.null(ONLY_CLUSTERS)){    
        clusters = clusters[clusters %in% ONLY_CLUSTERS]
        sobj = sobj[, sobj@meta.data[[CLUSTER_COL]] %in% clusters]
    }
}

if (! is.null(opt[['n-svs']])){
    N_SVS = opt[['n-svs']]
} else {
    N_SVS = 0
    RUN_SVA = FALSE
    SVA_COLS=NULL
    SVA_CTR_COLS = NULL
}

if (N_SVS > 0){
    RUN_SVA = TRUE
    if (is.null(opt[['sva-cols']])){
        stop("Error: Please provide sva-cols if N_SVS > 0")
    } else {
        SVA_COLS = strsplit(opt[['sva-cols']], ",")[[1]]
    }

    if (! is.null(opt[['sva-ctr-cols']])){
        SVA_CTR_COLS = strsplit(opt[['sva-ctr-cols']], ",")[[1]]
    } else {
        SVA_CTR_COLS = NULL
    }
} else {
    RUN_SVA = FALSE
    SVA_COLS=NULL
    SVA_CTR_COLS = NULL
}
all_cols = sort(unique(c(GROUPING_COLS,CONTRAST_COL,CLUSTER_COL,SAMPLE_COL,SVA_COLS,SVA_CTR_COLS)))
if (! all(all_cols %in% colnames(md))){
    stop(paste("Error: One or more grouping columns not found in metadata:\n", 
        paste(all_cols[! all_cols %in% colnames(md)], collapse=", ")))
}

# NA handling
if (any(colSums(is.na(md[all_cols])) > 0)){
    warning(paste(
    "\nWARNING! The following columns contain NA values:\n", 
    paste(all_cols[colSums(is.na(md[all_cols])) > 0], collapse = "\n")))
} else {
    print("No missing values found in grouping columns")
}

# if true, drop any samples that are missing values
if (is.null(opt[['drop-na']])){
    DROP_NA = TRUE
} else {
    DROP_NA = opt[['drop-na']]
}
if (DROP_NA){
    md_complete = md[, all_cols]
    md_complete = md_complete[complete.cases(md_complete), ]

    print(dim(md))
    print(dim(md_complete))

    rownames_to_drop = rownames(md)[!rownames(md) %in% rownames(md_complete)]
    if (length(rownames_to_drop) == 0){
        print("No cells with missing values to drop")
        rownames_to_drop = NULL
    } else {
        print(head(rownames_to_drop))
        warning(paste(
            "\nWARNING! Dropping", 
            length(rownames_to_drop), 
            "cells with missing values in any of",  
            paste(all_cols, collapse=", ")))
        sobj = sobj[, rownames(md_complete)]
    }

}

print("**************** RUNNING PSEUDOBULK CODE ****************")

counts_list = list()
metadata_list = list()

if (is.null(CLUSTER_COL)){
    cluster = "all"
    cat("Psuedobulking all cells\n")
    pb_list = pseudobulk_seurat(
        sobj=sobj, 
        grouping_cols=GROUPING_COLS, 
        assay=ASSAY, 
        min_n_cells=MIN_N_CELLS, 
        min_counts_gene=MIN_COUNTS_GENE,
        min_frac_gene=MIN_FRAC_GENE,
        contrast_col=CONTRAST_COL,
        cols_to_mean=COLS_TO_MEAN,
        cols_to_weighted_avg=COLS_TO_WEIGHTED_AVG,
        cols_to_median=COLS_TO_MEDIAN)

    counts_list[[cluster]] = pb_list[["counts"]] 
    metadata_list[[cluster]] = pb_list[["metadata"]]
} else {
    for (cluster in clusters){
        cat(paste("Processing cluster: ", cluster, "\n"))
        sobj_cluster = sobj[, sobj@meta.data[[CLUSTER_COL]] == cluster]
        pb_list = pseudobulk_seurat(
            sobj=sobj_cluster,
            grouping_cols=GROUPING_COLS, 
            assay=ASSAY, 
            min_n_cells=MIN_N_CELLS, 
            min_counts_gene=MIN_COUNTS_GENE,
            min_frac_gene=MIN_FRAC_GENE,
            contrast_col=CONTRAST_COL,
            cols_to_mean=COLS_TO_MEAN,
            cols_to_weighted_avg=COLS_TO_WEIGHTED_AVG,
            cols_to_median=COLS_TO_MEDIAN)
        counts_list[[cluster]] = pb_list[["counts"]]
        metadata_list[[cluster]] = pb_list[["metadata"]]
        print(head(pb_list[["metadata"]]))
    }
}

all_cols = unique(c(
    GROUPING_COLS,
    CLUSTER_COL,
    CONTRAST_COL, 
    SAMPLE_COL, 
    SVA_COLS, 
    SVA_CTR_COLS
))

for (cluster in names(metadata_list)){
    print(paste("Sanitizing metadata for cluster", cluster))
    md = metadata_list[[cluster]]
    for (col in all_cols){
        md = sanitize_col_values(md, col)
    }
    metadata_list[[cluster]] = md
}

obj_list = list()
if (OUTPUT_TYPE == "seurat"){
    for (cluster in names(counts_list)){
        cat("Creating Seurat object for Cluster: ", cluster, "\n")
        counts = counts_list[[cluster]]
        metadata = metadata_list[[cluster]]
        metadata$cluster = cluster
        sobj_pb = CreateSeuratObject(
            counts = counts, 
            meta.data = metadata)
        obj_list[[cluster]] = sobj_pb
    }
    merged_obj = mergeSeuratListWithMetaData(obj_list)

    if (RUN_SVA){
        # remove pre-existing cols starting with SV
        cd = merged_obj@meta.data
        cd = cd[, !grepl("^SV", toupper(colnames(cd)))]

        counts = GetAssayData(merged_obj, slot="counts", assay="RNA") # always want to run SVA on raw counts
        counts_low_expr_removed = counts[rowMeans(counts) > 1e-6, ] # remove genes with low expression
        counts_low_expr_removed = as.matrix(counts_low_expr_removed)
        # run SVA
        cd = get_df_with_svs(
            edata=counts_low_expr_removed,
            df=cd,
            cols=SVA_COLS,
            ctr_cols=SVA_CTR_COLS,
            n=N_SVS)
        
        # store the columns used to calculated SVs and the Null model
        cd$SVA_cols = paste0(SVA_COLS, collapse=",")
        cd$SVA_ctr_cols = paste0(SVA_CTR_COLS, collapse=",")
        
        merged_obj@meta.data = cd
    }

    cat("Saving Seurat object for Cluster: ", cluster, "\n")
    qsave(merged_obj, file.path(output_dir, paste0(slogan, "_pb_seurat__", cluster, "__", suffix, ".qs")))

} else if (OUTPUT_TYPE == "sce"){
    for (cluster in names(counts_list)){
        cat ("Creating SingleCellExperiment object for Cluster: ", cluster, "\n")
        counts = counts_list[[cluster]]
        metadata = metadata_list[[cluster]]
        sce = SingleCellExperiment(
            assays = list(counts = counts),
            colData = metadata)
        obj_list[[cluster]] = sce
    }

    if (RUN_SVA){

        for (cluster in names (obj_list)){
            cat("Running SVA for Cluster: ", cluster, "\n")
            sce = obj_list[[cluster]]
            cd = as.data.frame(colData(sce))
            cd = cd[, !grepl("^SV", toupper(colnames(cd)))]
            print(dim(sce))
            print(dim(counts(sce)))
            print(dim(cd))
            
            counts = as.matrix(counts(sce)) # always want to run SVA on raw counts
            counts_low_expr_removed = counts[rowMeans(counts) > 1e-6, ] # remove genes with low expression, otherwise SVA fails
            counts_low_expr_removed
            
            # run SVA
            cd = get_df_with_svs(
                edata=counts_low_expr_removed,
                df=cd,
                cols=SVA_COLS,
                ctr_cols=SVA_CTR_COLS,
                n=N_SVS)
            
            # store the columns used to calculated SVs and the Null model
            cd$sva_cols = paste0(SVA_COLS, collapse=",")
            cd$sva_ctr_cols = paste0(SVA_CTR_COLS, collapse=",")
            colData(sce) = DataFrame(cd)
            obj_list[[cluster]] = sce
        }
    }

    for (cluster in names(obj_list)){
        sce = obj_list[[cluster]]
        print(dim(sce))

        if (suffix == ""){
            output_basename = paste0(slogan, "__", cluster, "__pb_sce.qs")
        } else {
            output_basename = paste0(slogan, "__", cluster, "__pb_sce__", suffix, ".qs")
        }
        cat("Saving SingleCellExperiment object for Cluster: ", cluster, "\nto path: ", 
            file.path(output_dir, output_basename), "\n")
        qsave(sce, file.path(output_dir, output_basename))
    }
}