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
    kegg_2021_human = "~/genesets/KEGG_2021_Human.txt",
    go_process = "~/genesets/GO_Biological_Process_2021.txt",
    go_function = "~/genesets/GO_Molecular_Function_2021.txt",
    msigdb_hallmark = "~/genesets/MSigDB_Hallmark_2020.txt",
    syngo_ontologies = "~/genesets/syngo_ontologies.txt"
)

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(qs)))
suppressWarnings(suppressMessages(library(Seurat)))

suppressWarnings(suppressMessages(source("~/sc-online/plot.R")))
suppressWarnings(suppressMessages(source("~/sc-online/gsea.R")))

spec <- matrix(c(
    'path', 'p', 1, "character",
    'rank-col', 'r', 1, "character",
    'segment-by', 's', 1, "character",
    'abs', 'a', 0, "logical"
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

        gs=.sconline.GSEA.readGMT(
            file=geneset_path,
            bkg_genes=dataDE$gene_short_name,
            min.gs.size=15,max.gs.size=250)
        gsea = runGSEA(dataDE, gs, rank_col=RANK_COL, abs=ABS, desc=TRUE)
        gsea$gene_set = gene_set
        #gsea = gsea[which(gsea$padj<0.05),]

        if (nrow(gsea) > 0){
            gsea_list[[gene_set]] = gsea
        }
    }
    gsea_df = do.call(rbind, gsea_list)

    cat(paste0("\nWriting GSEA Results for ", de_slogan, "\n"))
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