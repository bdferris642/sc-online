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
# gene_sets = list(
#     kegg_2021_human = file.path(script_dir, "../genesets/KEGG_2021_Human.txt"),
#     go_process = file.path(script_dir, "../genesets/GO_Biological_Process_2021.txt"),
#     go_function = file.path(script_dir, "../genesets/GO_Molecular_Function_2021.txt"),
#     msigdb_hallmark = file.path(script_dir, "../genesets/MSigDB_Hallmark_2020.txt"),
#     syngo_ontologies = file.path(script_dir, "../genesets/syngo_ontologies.txt"))
#  #   cort_and_snap = file.path(script_dir, "../genesets/cort_and_snap.txt")   

gene_sets = list(
    kegg_2026_human = file.path(script_dir, "../genesets/latest/KEGG_2026.gmt"),
    reactome_2024 = file.path(script_dir, "../genesets/latest/Reactome_Pathways_2024.gmt"),
    gwas_2025 = file.path(script_dir, "../genesets/latest/GWAS_Catalog_2025.gmt"),
    jensen_disease_2025 = file.path(script_dir, "../genesets/latest/Jensen_DISEASES_Curated_2025.gmt"),
    trrust_2019 = file.path(script_dir, "../genesets/latest/TRRUST_Transcription_Factors_2019.gmt"),
    jaspar_tf_2025 = file.path(script_dir, "../genesets/latest/JASPAR_PWM_Human_2025.gmt"),
    go_bp_2025 = file.path(script_dir, "../genesets/latest/GO_Biological_Process_2025.gmt"),
    go_mf_2025 = file.path(script_dir, "../genesets/latest/GO_Molecular_Function_2025.gmt")
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

sanitize_column_name = function(colname){
    colname = gsub("\\.", "_", colname)
    colname = gsub("-", "_", colname)
    colname = gsub(" ", "_", colname)
    colname = gsub("\\.", "_", colname)
    colname = gsub("\\+", "_and_", colname)
    colname = gsub("\\&", "_and_", colname)
    colname = gsub("\\:", "_x_", colname)
    colname = gsub("\\?", "", colname)
    colname = gsub("\\/", "", colname)
    return(colname)
}

RANK_COL = ifelse(
    is.null(opt[['rank-col']]), 
    "t", 
    sanitize_column_name(opt[['rank-col']])
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
if ("signed_neg_log10_p_case_controlpd" %in% colnames(de_df)) {
    de_df = de_df %>% filter(is.finite(signed_neg_log10_p_case_controlpd))
}

colnames(de_df) = sapply(colnames(de_df), sanitize_column_name)

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

    # Redundancy filtering: applied independently per gene set collection
    filtered_list = list()
    for (gene_set in names(gsea_list)) {
        cat(paste0("\nApplying redundancy filter to gene set: ", gene_set, "\n"))
        gs_df = gsea_list[[gene_set]]
        gs_df_filtered = filter_redundant_gene_sets(gs_df)
        n_raw = sum(abs(gs_df$NES) >= 1.0 & gs_df$padj < 0.05)
        n_kept = sum(abs(gs_df_filtered$NES) >= 1.0 & gs_df_filtered$padj < 0.05)
        cat(paste0("  Significant sets: ", n_raw, " -> ", n_kept, " after pruning\n"))
        filtered_list[[gene_set]] = gs_df_filtered
    }
    gsea_df_filtered = do.call(rbind, filtered_list)

    cat(paste0("\nWriting Redundancy-Filtered GSEA Results for ", de_slogan,
               "\nto ", de_slogan, "_gsea_redundancy_filtered.csv\n"))
    write.csv(gsea_df_filtered,
              file.path(gsea_outdir, paste0(de_slogan, "_gsea_redundancy_filtered.csv")),
              row.names=FALSE)

    figure_outdir_filtered = file.path(base_path, "gsea", tolower(abs_str), "figures/png/redundancy_filtered")
    if (!dir.exists(figure_outdir_filtered)){
        dir.create(figure_outdir_filtered, recursive = TRUE)
    }

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

        # Plot redundancy-filtered results
        cat(paste0("\nPlotting Redundancy-Filtered GSEA Results for ", de_slogan, " ", gene_set, "\n"))
        gsea_df_subset_filtered = gsea_df_filtered[gsea_df_filtered$gene_set==gene_set & gsea_df_filtered$padj < 0.05,]

        if (nrow(gsea_df_subset_filtered) == 0){
            next
        }

        max_nes_f = ceiling(max(c(max(gsea_df_subset_filtered$NES), 0)))
        min_nes_f = floor(min(c(min(gsea_df_subset_filtered$NES), 0)))

        plot_gsea_result_hdot(
            gsea_df_subset_filtered,
            title = paste("GSEA NES by Pathway (Redundancy Filtered):", gene_set, gsub("_", " ", abs_str)),
            xlim=c(min_nes_f, max_nes_f),
            leading_edge_n=10,
            leading_edge_linebreak_n=5,
            top_n=10,
            fig_filename=file.path(figure_outdir_filtered, paste0(de_slogan, "_", gene_set, "_gsea_redundancy_filtered.png"))
        )
    }
}