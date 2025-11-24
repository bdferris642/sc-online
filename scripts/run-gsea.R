print("**************** LOADING LIBRARIES ****************")
# Detect script path when running via Rscript
args = commandArgs(trailingOnly = FALSE)
script_path = sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir = dirname(normalizePath(script_path))
  message("Script located in directory: ", script_dir)
} else {
  stop("Cannot determine script path. Are you running via Rscript?")
}

# a script that, given a .csv of DE outputs at --path with foldchanges in a column --rank-col  
# runs GSEA on the gene sets in the --gene-sets variable, saving these csvs into a `gsea_outdir`
# and plotting results to a `figure_outdir`.
# OJO: some arguments to .sconline.GSEA.readGMT (min.gs.size=15,max.gs.size=250) and runGSEA (abs=FALSE, desc=TRUE) are hard coded. TODO: parametrize!
# OJO: requires `gene_sets` to be a list of paths to gene sets in GMT format. TODO: parametrize!

# PATH is the path to a single CSV file containing differential expression results
# RANK_COL is the column in the CSV file that contains the rank values for GSEA
# SEGMENT_BY is an optional argument, a the column in the CSV file that contains the segment values for GSEA 
#   (in case you want to run GSEA on different segments of the data)

# add the paths to the gene sets you want to use here
gene_sets = list(
    kegg_2021_human = file.path(script_dir, "../genesets/KEGG_2021_Human.txt"),
    go_process = file.path(script_dir, "../genesets/GO_Biological_Process_2021.txt"),
    go_function = file.path(script_dir, "../genesets/GO_Molecular_Function_2021.txt"),
    msigdb_hallmark = file.path(script_dir, "../genesets/MSigDB_Hallmark_2020.txt"),
    syngo_ontologies = file.path(script_dir, "../genesets/syngo_ontologies.txt"),
    cort_and_snap = file.path(script_dir, "../genesets/cort_and_snap.txt")   
)

# latest_me_pd = "~/genesets/latest_me_pd_de.txt",
# cell_class_markers_1 = "~/genesets/dapi_nurr_merged_seurat_clean_ctr_markers__cell_type__avg_log2FC_gt_0.5_expr_ratio_gt_1.33.txt",
# cell_type_markers_1 = "~/genesets/dapi_nurr_merged_seurat_clean_ctr_markers__cell_class__avg_log2FC_gt_0.5_expr_ratio_gt_1.33.txt"
# handsaker_genes = "~/genesets/handsaker_genes.txt"

suppressWarnings(suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(glue)
    library(getopt)
    library(Matrix)
    library(qs)
    library(Seurat)
    source(file.path(script_dir, "../plot.R"))
    source(file.path(script_dir, "../gsea.R"))
}))

spec <- matrix(c(
    'path', 'p', 1, "character",
    'rank-col', 'r', 1, "character",
    'segment-by', 's', 1, "character",
    'abs', 'a', 0, "logical",
    'min-size', 'min', 1, "integer",
    'max-size', 'max', 1, "integer"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]

RANK_COL = ifelse(
    is.null(opt[['rank-col']]), 
    "t", 
    opt[['rank-col']]
)

if (is.null(opt[['segment-by']])){
    SEGMENT_BY = NULL
} else {
    SEGMENT_BY = opt[['segment-by']]
}

if (is.null(opt[['abs']])){
    ABS = FALSE
} else {
    ABS = opt[['abs']]
}

if (ABS){
    abs_str = "Absolute_Value"
} else {
    abs_str = "Signed"
}

if (is.null(opt[['min-size']])){
    MIN_SIZE = 15
} else {
    MIN_SIZE = opt[['min-size']]
}

if (is.null(opt[['max-size']])){
    MAX_SIZE = 250
} else {
    MAX_SIZE = opt[['max-size']]
}

basename = basename(PATH)
de_slogan = gsub(".csv", paste0("_", RANK_COL), basename)
base_path = dirname(PATH)

gsea_outdir = file.path(base_path, "gsea", tolower(abs_str))
figure_outdir = file.path(base_path, "gsea", tolower(abs_str), "figures/png")

if(!dir.exists(figure_outdir)){
    dir.create(figure_outdir, recursive = TRUE)
}

de_df = read.csv(PATH)
if (is.null(SEGMENT_BY)){
    dataDE_list = list(all=prep_df_for_gsea(de_df))
} else {
    dataDE_list = list()
    for (segment in unique(de_df[[SEGMENT_BY]])){
        print(segment)
        if(nrow(de_df[de_df[[SEGMENT_BY]]==segment,]) == 0){
            next
        }
        dataDE_list[[segment]] = prep_df_for_gsea(de_df[de_df[[SEGMENT_BY]]==segment,])
    }
}

for (segment in names(dataDE_list)){
    dataDE = dataDE_list[[segment]]
    if (!is.null(SEGMENT_BY)){
        de_slogan = segment
    }
    cat(paste0("\nRunning GSEA for ", de_slogan, "\n"))
    gsea_list = list()
    for (gene_set in names(gene_sets)){
        geneset_path = gene_sets[[gene_set]]

        # ojo: specific carve-outs for gene sets which have large pathway sizes 
        # that we nevertheless feel to be relevant
        if (gene_set %in% c("handsaker_genes", "latest_me_pd")){
            this_max_size = 50000
        } else {
            this_max_size = MAX_SIZE
        }

        gs=.sconline.GSEA.readGMT(
            file=geneset_path,
            bkg_genes=dataDE$gene_short_name,
            min.gs.size=MIN_SIZE, max.gs.size=this_max_size)
        message(paste0("Gene set list: ", gene_set, " has ", length(gs), " gene sets"))
        gsea = runGSEA(dataDE, gs, rank_col=RANK_COL, abs=ABS, desc=TRUE, min_size=MIN_SIZE, max_size=this_max_size)
        gsea$gene_set = gene_set
        #gsea = gsea[which(gsea$padj<0.05),]

        if (nrow(gsea) > 0){
            gsea_list[[gene_set]] = gsea
        } else {
            message(paste0("No genes in gene set: ", gene_set))
        }
    }
    gsea_df = do.call(rbind, gsea_list)

    cat(paste0("\nWriting GSEA Results for ", de_slogan, "\nto ", de_slogan, "_gsea.csv\n"))
    write.csv(gsea_df, file.path(gsea_outdir, paste0(de_slogan, "_gsea.csv")), row.names=FALSE)

    for (gene_set in unique(gsea_df$gene_set)){
        cat(paste0("\nPlotting GSEA Results for ", de_slogan, " ", gene_set, "\n"))
        gsea_df_subset = gsea_df[gsea_df$gene_set==gene_set & gsea_df$padj < 0.05,]

        if (nrow(gsea_df_subset) == 0){
            next
        }

        max_nes = ceiling(max(c(max(gsea_df_subset$NES), 0)))
        min_nes = floor(min(c(min(gsea_df_subset$NES), 0)))

        plot_gsea_result_hdot(
            gsea_df_subset, 
            title = paste("GSEA NES by Pathway:", gene_set, gsub("_", " ", abs_str)),
            xlim=c(min_nes, max_nes),
            leading_edge_n=10,
            leading_edge_linebreak_n=5,
            top_n=10,
            fig_filename=file.path(figure_outdir, paste0(de_slogan, "_", gene_set, "_gsea.png"))
        )
    }
}