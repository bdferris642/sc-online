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
    'assay', 'a', 1, 'character',
    'min-n-cells', 'mnc', 1, 'numeric',
    'min-counts-gene', 'mcg', 1, 'numeric',
    'min-frac-gene', 'mfg', 1, 'numeric',
    'suffix', 's', 1, 'character'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]

GROUPING_COLS = strsplit(opt[['grouping-cols']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
SAMPLE_COL = opt[['sample-col']]

cat(paste("Reading Seurat object from: ", PATH, "\n"))
sobj = qread(PATH)
md = sobj@meta.data
md[[SAMPLE_COL]] = paste0("sample_", md[[SAMPLE_COL]])
md[[SAMPLE_COL]] = gsub("-", "_", md[[SAMPLE_COL]])

if (is.null(opt[['cluster-col']])){
    CLUSTER_COL = NULL
} else {
    CLUSTER_COL = opt[['cluster-col']]
}

if (! is.null(CLUSTER_COL)){
    md[[CLUSTER_COL]] = gsub(" ", "_", md[[CLUSTER_COL]])
}
md[[CONTRAST_COL]] = factor(md[[CONTRAST_COL]], levels = sort(unique(md[[CONTRAST_COL]])))
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
output_dir = file.path(dirname(PATH), "pseudobulk")
if (! dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}

if (! is.null(ONLY_CLUSTERS)){
    clusters = unique(sobj@meta.data[[CLUSTER_COL]])
    clusters = clusters[clusters %in% ONLY_CLUSTERS]
    sobj = sobj[, sobj@meta.data[[CLUSTER_COL]] %in% clusters]
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
    merged_obj = .mycBindFn(obj_list)
    cat("Saving SingleCellExperiment object for Cluster: ", cluster, "\n")
    qsave(merged_obj, file.path(output_dir, paste0(slogan, "_pb_sce__", cluster, "__", suffix, ".qs")))
}