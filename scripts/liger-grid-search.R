library(RSpectra)
library(rliger)
library(Matrix)
library(irlba)
library(qs)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
source("~/sc-online/utils.R")

lambdas = c(1, 2.5, 5, 7.5, 10)
ks = c(10, 15, 20, 25, 30)

BASE_PATH = "/mnt/accessory/seq_data/pd_all/20240410"
out_dir = file.path(BASE_PATH, "liger_objs")
dir.create(out_dir)

print("Loading data")
nurr_da_clean_sce = qread(file.path(BASE_PATH, "sce_nurr_da_clean_harmonized_240416.qs"))
nurr_da_clean_ctr_sce = nurr_da_clean_sce[, nurr_da_clean_sce$case_control == "ctr"]
nurr_da_clean_ctr_sce_list = .mySplitObject(nurr_da_clean_ctr_sce, "participant_id")
gte_100_names = names(nurr_da_clean_ctr_sce_list)[sapply(nurr_da_clean_ctr_sce_list, ncol) >= 100]

print("Creating Liger Object")
nurr_da_clean_ctr_gte100_liger = createLiger(nurr_da_clean_ctr_sce_list[gte_100_names])

print("Normalizing, Selecting Variable Genes, Scaling")
nurr_da_clean_ctr_gte100_liger = (nurr_da_clean_ctr_gte100_liger
    %>% rliger::normalize()
    %>% rliger::selectGenes()
    %>% rliger::scaleNotCenter()
)


for (this_lambda in lambdas){
    for (this_k in ks){
        cat("\n\n")
        print(paste0("Running iNMF with k=", this_k, " and lambda=", this_lambda))
        
        suppressWarnings({
            nurr_da_clean_ctr_gte100_liger = runIntegration(nurr_da_clean_ctr_gte100_liger, k = this_k, lambda=this_lambda, method="iNMF")
            nurr_da_clean_ctr_gte100_liger = quantileNorm(nurr_da_clean_ctr_gte100_liger)
        })

        # get gene loadings
        # remove if they contain at least 3 bad genes
        # as defined by (1) starting with "RPL", "RPS", or "MT-" or (2) being in the preset list of bad genes
        dims_to_exclude = c()        
        w = getMatrix(nurr_da_clean_ctr_gte100_liger, "W")

        bad_genes = c(            
            "GAPDH", "HSP90A", "HSP90B", "UBB", "TUBB2A",
            rownames(w)[grepl("^MT-", rownames(w))], 
            rownames(w)[grepl("^RPL", rownames(w))],
            rownames(w)[grepl("^RPS", rownames(w))]
        )
        
        for (i in 1:this_k) {
            loadings <- w[,i]

            # Sort the loadings to find the highest values
            sorted_loadings <- sort(loadings, decreasing = TRUE)

            # Select the top 10 genes
            top_genes <- names(sorted_loadings)[1:15]

            # Check if any of the top genes are in the bad genes list
            if (sum(top_genes %in% bad_genes) >= 3) {
                dims_to_exclude <- c(dims_to_exclude, i)
            }
        }

        # run clustering and umap on the remaining dimensions
        dims_to_use = setdiff(1:this_k, dims_to_exclude)
        
        print("Running clustering and UMAP")
        suppressWarnings({
            nurr_da_clean_ctr_gte100_liger = runCluster(nurr_da_clean_ctr_gte100_liger, useDims=dims_to_use, resolution = 0.6, nNeighbors = 30)
            runUMAP(nurr_da_clean_ctr_gte100_liger, useDims=dims_to_use, n_neighbors=30, min_dist=0.3)
        })

        print("Saving object")
        qsave(
            nurr_da_clean_ctr_gte100_liger, 
            file.path(
                out_dir, 
                paste0("liger_nurr_da_clean_ctr_gte100_k_", this_k, "_lambda_", this_lambda, ".qs")
        ))
    }
}

