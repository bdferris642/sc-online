
library(data.table)
library(dplyr)
library(ggthemes)
library(ggplot2)
library(glue)
library(harmony)
library(patchwork)
library(purrr)
library(qs)
library(Matrix)
library(Seurat)
library(rcna)

source("~/sc-online/plot.R")
source("~/sc-online/gsea.R")

BASE_PATH = "/mnt/accessory/seq_data/pd_all/240514"
PATH_LIST = c(
    "nurr_da_double_clean.qs"
    # "dapi_mg_clean_ctr.qs"
    # "dapi_mg_clean.qs"
    # "dapi_astro_clean.qs",
    # "dapi_endo_clean.qs",
    # "dapi_da_clean_hi_numi.qs",
    # "dapi_oligo_seurat.qs",
    # "dapi_opc_seurat.qs",
    # "dapi_seurat_clean.qs",
    # "nurr_seurat_nonda_clean.qs"
)

REDUCTION_TO_USE = "harmony"


figure_outdir = file.path(BASE_PATH, "cna/figures/png")
gsea_outdir = file.path(BASE_PATH, "cna/gsea")

dir.create(figure_outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(gsea_outdir, recursive=TRUE, showWarnings=FALSE)

gene_sets = list(
    # add the paths to the gene sets you want to use here
    kegg_2021_human = "~/genesets/KEGG_2021_Human.txt",
    go_process = "~/genesets/GO_Biological_Process_2021.txt",
    go_function = "~/genesets/GO_Molecular_Function_2021.txt",
    msigdb_hallmark = "~/genesets/MSigDB_Hallmark_2020.txt",
    syngo_ontologies = "~/genesets/syngo_ontologies.txt"
)

for (path in PATH_LIST){

    cat("\nProcessing Object", path, "\n")
    
    obj = qread(file.path(BASE_PATH, path))
    slogan = sub(".qs$", "", path)

    if (length(obj@graphs) == 0){
        if (REDUCTION_TO_USE %in% names(obj@reductions)){
            obj = obj %>% FindNeighbors(dims = 1:20, reduction=REDUCTION_TO_USE)
        } else if ("pca" %in% names(obj@reductions)){
            print(paste(REDUCTION_TO_USE, " not found -- using PCA"))
            obj = obj %>% FindNeighbors(dims = 1:20, reduction="pca")
        } else {
            stop("No reduction found!")
        }
    }

    # obj$status_val = as.numeric(factor(obj$case_control, c('ctr', 'pd')))
    obj$status_val = as.numeric(factor(obj$source, c('GTEx', 'Calico')))
    obj$sex_val = as.numeric(factor(obj$sex, c('Female', 'Male')))
    
    cat("Running CNA against sex\n")
    obj = association.Seurat(
        seurat_object = obj, 
        test_var = 'sex_val', 
        samplem_key = 'participant_id', 
        graph_use = 'RNA_nn', 
        verbose = TRUE,
        batches = NULL,  ## no batch variables to include
        covs = c("status_val", "age"))
    
    obj[["cna_sex"]] = CreateDimReducObject(
        embeddings = obj[["cna"]]@cell.embeddings, 
        loadings = obj[["cna"]]@feature.loadings, 
        stdev = obj[["cna"]]@stdev,  
        key = "cna_sex",
        misc = obj[["cna"]]@misc
    )

    obj$cna_ncorrs_sex = obj$cna_ncorrs
    obj$cna_ncorrs_fdr10_sex = obj$cna_ncorrs_fdr10
    obj$cna_ncorrs_fdr05_sex = obj$cna_ncorrs_fdr05

    cat("Running CNA against age\n")
    obj = association.Seurat(
        seurat_object = obj, 
        test_var = 'age', 
        samplem_key = 'participant_id', 
        graph_use = 'RNA_nn', 
        verbose = TRUE,
        batches = NULL,  ## no batch variables to include
        covs = c("sex_val", "status_val"))

    obj[["cna_age"]] = CreateDimReducObject(
        embeddings = obj[["cna"]]@cell.embeddings, 
        loadings = obj[["cna"]]@feature.loadings, 
        stdev = obj[["cna"]]@stdev,  
        key = "cna_age",
        misc = obj[["cna"]]@misc
    )

    obj$cna_ncorrs_age = obj$cna_ncorrs
    obj$cna_ncorrs_fdr10_age = obj$cna_ncorrs_fdr10
    obj$cna_ncorrs_fdr05_age = obj$cna_ncorrs_fdr05

    cat("Running CNA against status_val\n")
    obj = association.Seurat(
        seurat_object = obj, 
        test_var = 'status_val', 
        samplem_key = 'participant_id', 
        graph_use = 'RNA_nn', 
        verbose = TRUE,
        batches = NULL,  ## no batch variables to include
        covs = c("sex_val", "age"))

    obj[["cna_status_val"]] = CreateDimReducObject(
        embeddings = obj[["cna"]]@cell.embeddings, 
        loadings = obj[["cna"]]@feature.loadings, 
        stdev = obj[["cna"]]@stdev,  
        key = "cna_status_val",
        misc = obj[["cna"]]@misc
    )

    obj$cna_ncorrs_status_val = obj$cna_ncorrs
    obj$cna_ncorrs_fdr10_status_val = obj$cna_ncorrs_fdr10
    obj$cna_ncorrs_fdr05_status_val = obj$cna_ncorrs_fdr05

    cat("Saving Data\n")
    qsave(obj, file.path(BASE_PATH, path))

    options(repr.plot.width=14, repr.plot.height=4)

    cat("\nPlotting Case Control CNA Results")
    p1 = FeaturePlot(obj, features = c('cna_ncorrs_status_val'))[[1]] + 
        scale_color_gradient2_tableau() + 
        labs(
            title = 'CNA disease association', color = 'Correlation',
            subtitle = sprintf('global p=%0.3f', obj@reductions$cna_status_val@misc$p)
        ) + 
    FeaturePlot(obj, features = c('cna_ncorrs_fdr10_status_val'))[[1]] + 
        scale_color_gradient2_tableau() + 
        labs(title = 'CNA disease association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
    DimPlot(obj, group.by = 'status_val')[[1]] +
        scale_color_tableau() + labs(title = 'Disease status') + 
    plot_layout(nrow = 1) + 
    plot_annotation(tag_levels = 'a')
    ggsave(
        filename=file.path(figure_outdir, paste0(slogan, "_cna_status_val.png")), 
        plot=p1, width = 14, height = 4, dpi=500, bg="white")

    cat("\nPlotting Age CNA Results")
    p2 = FeaturePlot(obj, features = c('cna_ncorrs_age'))[[1]] + 
        scale_color_gradient2_tableau() + 
        labs(
            title = 'CNA age association', color = 'Correlation',
            subtitle = sprintf('global p=%0.3f', obj@reductions$cna_age@misc$p)
        ) + 
    FeaturePlot(obj, features = c('cna_ncorrs_fdr10_age'))[[1]] + 
        scale_color_gradient2_tableau() + 
        labs(title = 'CNA age association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
    FeaturePlot(obj, features = c('age'))[[1]] +
        scale_color_gradient2_tableau() + 
        labs(
            title = 'Age', color = 'Age'
        ) + 
    plot_layout(nrow = 1) + 
    plot_annotation(tag_levels = 'a')
    ggsave(
        filename=file.path(figure_outdir, paste0(slogan, "_cna_age.png")), 
        plot=p2, width = 14, height = 4, dpi=500, bg="white")

    cat("\nPlotting Sex CNA Results")
    p3 = FeaturePlot(obj, features = c('cna_ncorrs_sex'))[[1]] + 
        scale_color_gradient2_tableau() + 
        labs(
            title = 'CNA sex association', color = 'Correlation',
            subtitle = sprintf('global p=%0.3f', obj@reductions$cna_sex@misc$p)
        ) + 
    FeaturePlot(obj, features = c('cna_ncorrs_fdr10_sex'))[[1]] + 
        scale_color_gradient2_tableau() + 
        labs(title = 'CNA sex association', subtitle = 'Filtered for FDR<0.10', color = 'Correlation') + 
    DimPlot(obj, group.by = 'sex')[[1]] +
        scale_color_tableau() + labs(title = 'Sex') + 
    plot_layout(nrow = 1) + 
    plot_annotation(tag_levels = 'a')    
    ggsave(
        filename=file.path(figure_outdir, paste0(slogan, "_cna_sex.png")), 
        plot=p3, width = 14, height = 4, dpi=500, bg="white")


    cat("Densifying Data Matrix\n")
    obj_data = t(as.matrix(obj@assays$RNA@data))

    cat("Correlating Gene Expression to Neighborhood Coefficients and Saving\n")

    ncorrs_status_val = cor(obj_data, obj$cna_ncorrs_status_val, method='spearman')
    status_val_DE = prep_df_for_gsea(data.frame(gene = rownames(ncorrs_status_val),r = ncorrs_status_val))
    status_val_DE = status_val_DE[!is.na(status_val_DE$r),]

    ncorrs_age = cor(obj_data, obj$cna_ncorrs_age, method='spearman')
    age_DE = prep_df_for_gsea(data.frame(gene = rownames(ncorrs_age), r = ncorrs_age))
    age_DE = age_DE[!is.na(age_DE$r),]

    ncorrs_sex = cor(obj_data, obj$cna_ncorrs_sex, method='spearman')
    sex_DE = prep_df_for_gsea(data.frame(gene = rownames(ncorrs_sex),r = ncorrs_sex))
    sex_DE = sex_DE[!is.na(sex_DE$r),]

    qsave(ncorrs_status_val, file.path(BASE_PATH, 'cna', paste0(slogan, "_ncorrs_status_val.qs")))
    qsave(ncorrs_age, file.path(BASE_PATH, 'cna', paste0(slogan, "_ncorrs_age.qs")))
    qsave(ncorrs_sex, file.path(BASE_PATH, 'cna', paste0(slogan, "_ncorrs_sex.qs")))

    # Run GSEA
    cat("\nPrepping Data for GSEA")
    de_files = list(status_val_DE, age_DE, sex_DE)

    names(de_files) = c(
        paste0("gsea__", slogan, "__ncorrs_status_val"),
        paste0("gsea__", slogan, "__ncorrs_age"),
        paste0("gsea__", slogan, "__ncorrs_sex")
    )

    for (de_slogan in names(de_files)){

        cat(paste0("\nRunning GSEA on ", de_slogan))

        dataDE = de_files[[de_slogan]]
        anno=.extraHumanGeneAnnoAdderFn(inputGeneNames=dataDE$gene)
        anno=anno[match(dataDE$gene,anno$gene_id),]
        dataDE=cbind(dataDE,anno)
        dataDE=dataDE[which(dataDE$gene_biotype=="protein_coding"),]

        gsea_list = list()

        for (gene_set in names(gene_sets)){
            geneset_path = gene_sets[[gene_set]]

            gs=.sconline.GSEA.readGMT(
                file=geneset_path,
                bkg_genes=dataDE$gene_short_name,
                min.gs.size=8,max.gs.size=500)
            gsea = runGSEA(dataDE, gs, rank_col="r", abs=FALSE, desc=TRUE)
            gsea$gene_set = gene_set
            gsea = gsea[which(gsea$padj<0.05),]

            if (nrow(gsea) > 0){
                gsea_list[[gene_set]] = gsea
            }
        }
        gsea_df = do.call(rbind, gsea_list)

        cat(paste0("\nWriting GSEA Results for ", de_slogan))
        write.table(gsea_df, file.path(gsea_outdir, paste0(de_slogan, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)

        for (gene_set in unique(gsea_df$gene_set)){
            cat(paste0("\nPlotting GSEA Results for ", de_slogan, " ", gene_set))
            gsea_df_subset = gsea_df[gsea_df$gene_set==gene_set,]

            if (nrow(gsea_df_subset) == 0){
                next
            }

            max_nes = ceiling(max(c(max(gsea_df_subset$NES), 0)))
            min_nes = floor(min(c(min(gsea_df_subset$NES), 0)))

            plot_gsea_result_hdot(
                gsea_df_subset, 
                title = paste0("GSEA NES by Pathway: ", gene_set),
                xlim=c(min_nes, max_nes),
                leading_edge_n=10,
                leading_edge_linebreak_n=5,
                top_n=10,
                fig_filename=file.path(figure_outdir, paste0(de_slogan, "_", gene_set, "_gsea.png"))
            )
        }
    }
}