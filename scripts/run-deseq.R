# Runs DE and, optionally, GSEA on a Seurat object in a .qs file
# DE Requirements:
# path points to a .qs of a seurat object containing counts with integer values, 
# along with metadata columns nUMI, pct_mito, and pct_intronic 
# (these can be dummy columns filled with 0s but right now the code will fail if they are not included, TODO: fix this in the future)

# --path: a path to a seurat object
# --contrast-col: column name to compare ( e.g. case_control, status)
# --cluster-col: column containing groups to run DE on (e.g. if cluster-col=cell_class, DE will be run once per unique value in this column)
# --sample-col: column containing labels for individuals (e.g. participant_id, sample_id)
# --covariates: a comma-separated string of covariates to use in the analysis (e.g. age,sex,pct_mito,pct_intronic,cell_type)
# --grouping-cols: a comma-separated string of column names that are the same within each participant 
#   (e.g. age,sex . pct_mito and pct_intronic can be different from cell to cell, so you don't want to group by them when you bulk).
#   TODO: in the future, compute these automatically
# --suffix: a suffix to append to end of file names just before .csv, defaulting to date
# --only-clusters: a comma-separated string of clusters within cluster_col. If not null, only runs DE on these clusters (e.g. oligo,opc)
# --assay: the Seurat assay to use. Defaults to RNA
# --min-n-cells: the minimum number of cells in a bulk in order for it to be included. Defaults to 10
# --min-counts-gene: the minimum counts in a gene in order for it to be analyzed. Defaults to 10 
#   (OJO: this is low, maybe compute this automatically in the future to be sensitive to dataset size / sequencing depth)
# --min-frac-gene: the minimum fraction of cells in which a must be is expressed in order to be analyzed. Defaults to 0.01
# --run-gsea: if TRUE, runs GSEA on the DE results. Defaults to TRUE


suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(qs)))
suppressMessages(suppressWarnings(library(Seurat)))

spec <- matrix(c(
    'path', 'p', 1, "character",
    'contrast-col', 'cc', 1, 'character',
    'cluster-col', 'cl', 1, 'character',
    'sample-col', 'sc', 1, 'character',
    'covariates', 'c', 1, 'character',
    'grouping-cols', 'g', 1, 'character',
    'suffix', 's', 1, 'character',
    'only-clusters', 'oc', 1, 'character',
    'assay', 'a', 1, 'character',
    'min-n-cells', 'mnc', 1, 'numeric',
    'min-counts-gene', 'mcg', 1, 'numeric',
    'min-frac-gene', 'mfg', 1, 'numeric',
    'run-gsea', 'rg', 0, 'logical'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
COVARIATES = strsplit(opt[['covariates']], ",")[[1]]
GROUPING_COLS = strsplit(opt[['grouping-cols']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
CLUSTER_COL = opt[['cluster-col']]
SAMPLE_COL = opt[['sample-col']]
MIN_SAMPLES_PER_CONDITION = 3 # hard code for now, a cluster needs at least this many samples in each condition with at least MIN_N_CELLS

cat(paste("Reading Seurat object from: ", PATH, "\n"))
sobj = qread(PATH)
md = sobj@meta.data
md[[SAMPLE_COL]] = paste0("sample_", md[[SAMPLE_COL]])
md[[SAMPLE_COL]] = gsub("-", "_", md[[SAMPLE_COL]])
if (! is.null(CLUSTER_COL)){
    md[[CLUSTER_COL]] = gsub(" ", "_", md[[CLUSTER_COL]])
}
md[[CONTRAST_COL]] = factor(md[[CONTRAST_COL]], levels = sort(unique(md[[CONTRAST_COL]])))
if ("age" %in% colnames(md)){
    md$age = as.numeric(md$age)
}
sobj@meta.data = md

ONLY_CLUSTERS = if(is.null(opt[['only-clusters']])){
    ONLY_CLUSTERS = NULL
} else {
    ONLY_CLUSTERS = strsplit(opt[['only-clusters']], ",")[[1]]
}
if (is.null(opt[['suffix']])){
    SUFFIX = format(Sys.Date(), "%y_%m_%d") 
} else {
    SUFFIX = opt[['suffix']]
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

if(is.null(opt[['run-gsea']])){
    RUN_GSEA = TRUE
} else {
    RUN_GSEA = opt[['run-gsea']]
}

DESIGN_FORMULA = as.formula(paste0("~", paste(COVARIATES, collapse="+"), "+", CONTRAST_COL))

cat(paste("path = ", PATH, "\n"))
cat(paste("contrast-col = ", CONTRAST_COL, "\n"))
cat(paste("cluster-col = ", CLUSTER_COL, "\n"))
cat("covariates:\n", COVARIATES, "\n")
cat("grouping-cols:\n", GROUPING_COLS, "\n")
cat(paste("sample-col = ", SAMPLE_COL, "\n"))
cat(paste("suffix = ", SUFFIX, "\n"))
cat("only-clusters:\n", ONLY_CLUSTERS, "\n")
cat(paste("assay = ", ASSAY, "\n"))
cat(paste("min-n-cells = ", MIN_N_CELLS, "\n"))
cat(paste("min-counts-gene = ", MIN_COUNTS_GENE, "\n"))
cat(paste("min-frac-gene = ", MIN_FRAC_GENE, "\n"))
cat(paste("design-formula = ", DESIGN_FORMULA, "\n"))

create_df_with_contrast_col = function() {
  setNames(data.frame(factor(levels = sort(unique(md[[CONTRAST_COL]])))), CONTRAST_COL)
}

pseudobulk_seurat = function(sobj, grouping_cols=GROUPING_COLS,
    assay="RNA", min_n_cells = 10, min_counts_gene = 10, min_frac_gene = 0.01){
    
    df = sobj@meta.data
    if (! all(GROUPING_COLS %in% colnames(df))){
        stop(paste0("Not all grouping columns are present in the Seurat object!\n",
            "Missing Columns: ", paste(GROUPING_COLS[!GROUPING_COLS %in% colnames(df)], collapse=", ")))
    }
    df$grouping = apply(df[grouping_cols], 1, function(row) paste(row, collapse = "_"))    
    sobj$grouping = df$grouping
    grouping_cols = c(grouping_cols, "grouping")

    # group by participant_id
    # sum nUMI 
    # weighted average of pct_intronic and pct_mito (weighted by nUMI)
    # take log10 of summed nUMI
    # the following columns should be identical for all cells in the same participant:
    # age, sex, brain_bank, case_control
    cat(paste("Grouping By Columns: ", paste(grouping_cols, collapse=", "), "\n"))

    df_bulk = (df 
        %>% group_by(across(all_of(grouping_cols)))
        %>% summarise(
            n_cells = n(),
            sum_nUMI = sum(nUMI),
            pct_mito = weighted.mean(pct_mito, nUMI),
            pct_intronic = weighted.mean(pct_intronic, nUMI))
        %>% mutate(log10_nUMI = log10(sum_nUMI))
    )
    df_bulk = df_bulk[df_bulk$n_cells >= min_n_cells,]
    if (nrow(df_bulk) == 0){
        return(list(counts=NULL,  metadata=create_df_with_contrast_col()))
    }

    sobj = sobj[,df$grouping %in% df_bulk$grouping]
    
    counts_bulk = AggregateExpression(sobj, group.by = "grouping")[[assay]]
    
    counts_orig = sobj@assays[[assay]]@counts
    counts_bulk = counts_bulk[rowSums(counts_orig)>= min_counts_gene & rowMeans(counts_orig > 0) >= min_frac_gene,]
    
    colnames(counts_bulk) = gsub("-", '_', colnames(counts_bulk))
    df_bulk = df_bulk[match(colnames(counts_bulk), df_bulk$grouping),]
    
    if (nrow(df_bulk) == 0){
        return(list(counts=NULL, metadata=create_df_with_contrast_col()))
    }

    return(
        list(counts=counts_bulk, metadata=df_bulk)
    )
}
run_deseq = function(sobj, design_formula=DESIGN_FORMULA){
    cat("Pseudobulking Seurat Object\n")
    pseudobulk_list = pseudobulk_seurat(sobj)

    # check whether there are at least MIN_SAMPLES_PER_CONDITION samples in each condition of contrast col
    t = table(pseudobulk_list$metadata[[CONTRAST_COL]])
    if (any(t < MIN_SAMPLES_PER_CONDITION)){
        cat(paste("Skipping DESeq2 for this contrast because there are not enough samples in each condition!"))
        print(t)
        return(NULL)
    }

    cat(paste("Running DESeq2 with design formula: ", as.character(design_formula), "\n"))
    dds = DESeqDataSetFromMatrix(
        countData = pseudobulk_list$counts,
        colData = pseudobulk_list$metadata,
        design = design_formula
    )
    dds = DESeq(dds)
    cat("DESeq Complete. Coefficient Names:\n")
    cat(resultsNames(dds))

    result_list = list()
    for (coef in resultsNames(dds)){
        if (tolower(coef) == "intercept"){
            next
        }
        cat("\n")
        cat(paste("Running lfcShrink for coef: ", coef, "\n"))
        res_shrunk = suppressMessages(suppressWarnings(lfcShrink(dds, coef=coef, type="apeglm")))
        res_shrunk = res_shrunk[!is.na(res_shrunk$padj),]
        res_shrunk = res_shrunk[order(-res_shrunk$log2FoldChange),]

        res_shrunk$gene = rownames(res_shrunk)
        res_shrunk$logFC = res_shrunk$log2FoldChange
        
        res_shrunk$log2FoldChange = NULL
        rownames(res_shrunk) = NULL

        res_shrunk = res_shrunk[,c("gene", "logFC", "padj", "baseMean", "lfcSE", "pvalue")]

        result_list[[coef]] = res_shrunk
    }
    return(result_list)
}

write_list = function(x, fil){ 
    z = deparse(substitute(x))
    cat(z, "\n", file=fil)
    nams=names(x) 
    for (i in seq_along(x) ){ 
        cat(nams[i], "\t",  x[[i]], "\n", file=fil, append=TRUE)
    }
}

GSEA_SCRIPT_PATH = "/home/ferris/sc-online/scripts/run-gsea.R"

################################# MAIN #################################

dirname = dirname(PATH)
basename = basename(PATH)
de_dir = file.path(dirname, "de/deseq", gsub(".qs", "", basename))
if (!dir.exists(de_dir)){
    dir.create(de_dir, recursive = TRUE)
}

if (is.null(CLUSTER_COL)){
    cat("Running DESeq2 for all cells\n")
    result_list = run_deseq(sobj)
    for (coef in names(result_list)){
        cat(paste("Writing results for coef: ", coef, '\n'))
        cat(paste("Number of Genes in Result: ", nrow(result_list[[coef]]), '\n'))
        
        suffix_path = file.path(de_dir, 'all_cells', coef, SUFFIX, paste0("deseq__", coef, "__", SUFFIX, ".csv"))
        latest_path = file.path(de_dir, 'all_cells', coef, "latest", paste0("deseq__", coef, "__latest.csv"))

        if(! dir.exists(dirname(suffix_path))){
            dir.create(dirname(suffix_path), recursive = TRUE)
        }
        if(! dir.exists(dirname(latest_path))){
            dir.create(dirname(latest_path), recursive = TRUE)
        }
        write.csv(result_list[[coef]], suffix_path, row.names=FALSE)
        write.csv(result_list[[coef]], latest_path, row.names=FALSE)
        
        if (RUN_GSEA){
            # Run the GSEA script if the coef is age or the contrast column
            if (tolower(coef) == "age" | grepl(CONTRAST_COL, coef)){
                system(paste0("Rscript ", GSEA_SCRIPT_PATH, " --path=", suffix_path))
            }
        }
    }
} else {
    clusters = sort(unique(md[[CLUSTER_COL]]))
    if (!is.null(ONLY_CLUSTERS)){
        clusters = clusters[clusters %in% ONLY_CLUSTERS]
    }
    for (cluster in clusters){
        if (sum(md[[CLUSTER_COL]] == cluster) == 0){
            cat(paste("Skipping cluster ", cluster, " because it has no cells!\n"))
            next
        }

        cat(paste("Running DESeq2 for cluster: ", cluster, '\n'))
        sobj_sub = sobj[,md[[CLUSTER_COL]] == cluster]


        result_list = run_deseq(sobj_sub)
        if (is.null(result_list)){
            next
        }
        for (coef in names(result_list)){
            cat(paste("Writing results for coef: ", coef, '\n'))
            cat(paste("Number of Genes in Result: ", nrow(result_list[[coef]]), '\n'))
            
            suffix_path = file.path(de_dir, CLUSTER_COL, coef, SUFFIX, paste0("deseq__", CLUSTER_COL, "__", cluster, "__", coef, "__", SUFFIX, ".csv"))
            latest_path = file.path(de_dir, CLUSTER_COL, coef, "latest", paste0("deseq__", CLUSTER_COL, "__", cluster, "__", coef, "__latest.csv"))
            if(! dir.exists(dirname(suffix_path))){
                dir.create(dirname(suffix_path), recursive = TRUE)
            }
            if (! dir.exists(dirname(latest_path))){
                dir.create(dirname(latest_path), recursive = TRUE)
            }
            cat("Writing to: ", suffix_path, "\n")
            write.csv(result_list[[coef]], suffix_path, row.names=FALSE)

            cat("Writing to: ", latest_path, "\n")
            write.csv(result_list[[coef]], latest_path, row.names=FALSE)

            if (RUN_GSEA){
                # Run the GSEA script if the coef is age or the contrast column
                if (tolower(coef) == "age" | grepl(CONTRAST_COL, coef)){
                    system(paste0("Rscript ", GSEA_SCRIPT_PATH, " --path=", suffix_path), ignore.stdout = TRUE)
                }
            }
        }
    }
}


if (is.null(CLUSTER_COL)){
   param_cluster_col = "all_cells"
} else {
    param_cluster_col = CLUSTER_COL
}
params = list(
    sobj_path = PATH,
    contrast_col = CONTRAST_COL,
    cluster_col = param_cluster_col,
    covariates = paste(COVARIATES, collapse=", "),
    sample_col = SAMPLE_COL,
    grouping_cols = paste(GROUPING_COLS, collapse=", "),
    suffix = SUFFIX,
    only_clusters = paste(ONLY_CLUSTERS, collapse=", "),
    assay = ASSAY,
    min_n_cells = MIN_N_CELLS,
    min_counts_gene = MIN_COUNTS_GENE,
    min_frac_gene = MIN_FRAC_GENE,
    n_case_cells = sum(! tolower(md[[CONTRAST_COL]]) %in% c("control", 'ctr', 'ctrl')),
    n_control_cells = sum(tolower(md[[CONTRAST_COL]]) %in% c("control", 'ctr', 'ctrl'))
)
write_list(params, file.path(de_dir, paste0("params__", param_cluster_col, "__", SUFFIX, ".txt")))
# remove the latest file
if (file.exists(file.path(de_dir, paste0("params__", param_cluster_col, "__latest.txt")))){
    file.remove(file.path(de_dir, paste0("params__", param_cluster_col, "__latest.txt")))
}
write_list(params, file.path(de_dir, paste0("params__", param_cluster_col, "__latest.txt")))