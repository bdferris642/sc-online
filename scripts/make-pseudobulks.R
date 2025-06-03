# INPUTS: A path to a seurat object and some pseudobulking parameters
# OUTPUTS: A list of pseudobulks (Seurat or SingleCellExperiment objects) saved to disk as .qs files

# The number of pseudobulks in the list will correspond to the number of distinct values of CLUSTER_COL in the Seurat object
# Each Pseudobulk in the list will have at least as many columns as there are distinct groupings of GROUPING_COLS in the input data. 

print("**************** LOADING LIBRARIES ****************")
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

suppressMessages(suppressWarnings(source("/home/ferris/sc-online/de.R")))
suppressMessages(suppressWarnings(source("/home/ferris/sc-online/utils.R")))
suppressMessages(suppressWarnings(source("/home/ferris/sc-online/getData.R")))

print("**************** PARSING ARGUMENTS ****************")
spec <- matrix(c(
    'path', 'p', 1, "character",
    'output-type', 'o', 1, 'character',
    'contrast-col', 'cc', 1, 'character',
    'cluster-col', 'cl', 1, 'character',
    'sample-col', 'sc', 1, 'character',
    'grouping-cols', 'g', 1, 'character',
    'only-clusters', 'oc', 1, 'character',
    'filter-samples', 'fs', 1, 'character',
    'assay', 'a', 1, 'character',
    'min-n-cells', 'mnc', 1, 'numeric',
    'min-counts-gene', 'mcg', 1, 'numeric',
    'min-frac-gene', 'mfg', 1, 'numeric',
    'suffix', 's', 1, 'character',
    'sva-cols', 'sva', 1, 'character',
    'sva-ctr-cols', 'svac', 1, 'character',
    'n-svs', 'nsv', 1, 'integer',
    'numi-col', 'numi', 1, 'character',
    'mito-col', 'mito', 1, 'character',
    'intronic-col', 'intr', 1, 'character'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]

GROUPING_COLS = strsplit(opt[['grouping-cols']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
SAMPLE_COL = opt[['sample-col']]

cat(paste("Reading Seurat object from: ", PATH, "\n"))
sobj = load_obj(PATH)
DefaultAssay(sobj) = "RNA"

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

if(is.null(opt[['numi-col']])){
    NUMI_COL = "nCount_RNA"
} else {
    NUMI_COL = opt[['numi-col']]
}
if (! NUMI_COL %in% colnames(md)){
    stop(paste("Error: nUMI column not found in metadata: ", NUMI_COL))
}
md$nUMI = md[[NUMI_COL]]

if(is.null(opt[['mito-col']])){
    MITO_COL = "pct_mito"
} else {
    MITO_COL = opt[['mito-col']]
}
if (! MITO_COL %in% colnames(md)){
    warning(
        paste("WARNING: Mitochondrial percentage column not found in metadata: ", MITO_COL,
        "\nCalculating pct_mito from gene names starting with MT-"))
    
    mito_genes = rownames(sobj)[grepl("^MT-", toupper(rownames(sobj)))]

    if (length(mito_genes) == 0){
        warning("WARNING: No mitochondrial genes found in the Seurat object!\nSetting pct_mito to 0 for all cells!")
        md$pct_mito = 0
    } else{
        counts = GetAssayData(sobj, slot="counts", assay="RNA")
        md$pct_mito = 100 * colSums(counts[mito_genes, ]) / colSums(counts)
    } 
} else {
    md$pct_mito = md[[MITO_COL]]
}

if(is.null(opt[['intronic-col']])){
    INTRONIC_COL = "pct_intronic"
} else {
    INTRONIC_COL = opt[['intronic-col']]
}
if (! INTRONIC_COL %in% colnames(md)){
    warning(
        paste("WARNING: Intronic percentage column not found in metadata: ", INTRONIC_COL, 
              "\nSetting pct_intronic to 0 for all cells!"))
    md$pct_intronic = 0
} else {
    md$pct_intronic = md[[INTRONIC_COL]]
}

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
    MIN_COUNTS_GENE = 10
} else {
    MIN_COUNTS_GENE = opt[['min-counts-gene']]
}
if(is.null(opt[['min-frac-gene']])){
    MIN_FRAC_GENE = 0.01
} else {
    MIN_FRAC_GENE = opt[['min-frac-gene']]
}

if(is.null(opt[['suffix']])){
    suffix = ""
} else {
    suffix = opt[['suffix']]
}

slogan = gsub(".qs", "", basename(PATH))
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
}
all_cols = sort(unique(c(GROUPING_COLS,CONTRAST_COL,CLUSTER_COL,SAMPLE_COL,SVA_COLS,SVA_CTR_COLS)))
if (! all(all_cols %in% colnames(md))){
    stop(paste("Error: One or more grouping columns not found in metadata:\n", paste(all_cols[! all_cols %in% colnames(md)], collapse=", ")))
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
        contrast_col=CONTRAST_COL)

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
            contrast_col=CONTRAST_COL)
        counts_list[[cluster]] = pb_list[["counts"]]
        metadata_list[[cluster]] = pb_list[["metadata"]]
        print(head(pb_list[["metadata"]]))
    }
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

        # run SVA
        cd = get_df_with_svs(
            edata=as.matrix(GetAssayData(merged_obj, slot="counts", assay="RNA")), # always want to run SVA on raw counts
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
            
            # run SVA
            cd = get_df_with_svs(
                edata=as.matrix(counts(sce)),
                df=cd,
                cols=SVA_COLS,
                ctr_cols=SVA_CTR_COLS,
                n=N_SVS)
            
            # store the columns used to calculated SVs and the Null model
            cd$sva_cols = paste0(SVA_COLS, collapse=",")
            cd$sva_ctr_cols = paste0(SVA_CTR_COLS, collapse=",")
            colData(sce) = DataFrame(cd)
            obj_list[[cluster]] = sce

            cat("Saving SingleCellExperiment object for Cluster: ", cluster, "\n")
            if (suffix == ""){
                output_basename = paste0(slogan, "__", cluster, "__pb_sce.qs")
            } else {
                output_basename = paste0(slogan, "__", cluster, "__pb_sce__", suffix, ".qs")
            }
            qsave(sce, file.path(output_dir, output_basename))
        }
    }

    # cat("Saving SingleCellExperiment object for Cluster: ", cluster, "\n")
    # if (suffix == ""){
    #     output_basename = paste0(slogan, "__pb_sce.qs")
    # } else {
    #     output_basename = paste0(slogan, "__pb_sce_", suffix, ".qs")
    # }
    # qsave(merged_obj, file.path(output_dir, output_basename))
}