BASE_DIR = "/mnt/accessory/seq_data/calico"
LT_OUTPUT_DIR_NAME = "lt_subtype_nurr_20240130"

SOURCE_PATH = "/mnt/accessory/seq_data/tkamath/nurr_ctr_seurat.qs"
SOURCE_DS_NAME = "tk"

TARGET_PATH = "/mnt/accessory/seq_data/calico/seurat_nurr_merged_initial_harmonized_20240130.qs"
TARGET_DS_NAME = "me"

#LT_COL = "cell_class"
LT_COL = "liger_cluster"

CELL_CLASSES = list(
    #"astro", "endo", 
    #"da", "immune", 
    #"oligo", "opc", 
    "nonda"
)

MIN_NUM_CELLS = 30

source("~/sc-online/extraGenes.R")
source("~/sc-online/myExpSetCreatorFn.R")
source("~/sc-online/utils.R")
source("~/sc-online/getData.R")
source("~/sc-online/plot.R")
source("~/sc-online/clust.R")
source("~/sc-online/labelTransfer.R")

library("DropletUtils")
library(future)
library(future.apply)
library(dplyr)
library(ggplot2)
library(harmony)
library(Matrix)
library(qs)
library(rlang)
library(SingleCellExperiment)
library(Seurat)
library(tidyr)

print("LOADING DATA")
source_seurat_obj_all = qread(SOURCE_PATH)
target_seurat_obj_all = qread(TARGET_PATH)

source_seurat_obj_all$dataset = SOURCE_DS_NAME
target_seurat_obj_all$dataset = TARGET_DS_NAME

print("ITERATING THROUGH CELL CLASSES")
for (cc in CELL_CLASSES){

    print(paste('CELL CLASS:',cc))

    source_seurat_obj = source_seurat_obj_all[, source_seurat_obj_all$cell_class == cc]
    target_seurat_obj = target_seurat_obj_all[, target_seurat_obj_all$cell_class == cc]

    # combine all donors with fewer than MIN_NUM_CELLS into a combined superdonor
    target_seurat_obj$donor_id_adj = target_seurat_obj$donor_id
    donor_counts = table(target_seurat_obj$donor_id)
    donor_ids_to_combine = names(donor_counts[donor_counts < MIN_NUM_CELLS])
    target_seurat_obj$donor_id_adj[target_seurat_obj$donor_id %in% donor_ids_to_combine] = "combined"

    print('GETTING VARIABLE FEATURES')
    # TODO: parametrize var features args?
    source_vf = getSeuratVarFeatureIntersectByCol(
        source_seurat_obj,
        subset_col="donor_id",
        original_nfeatures=5000,
        n_subsets_to_cover=NULL
    )

    target_vf = getSeuratVarFeatureIntersectByCol(
        target_seurat_obj,
        subset_col="donor_id",
        original_nfeatures=5000,
        n_subsets_to_cover=NULL
    )

    merged_vf = union(source_vf, target_vf)
    print(paste("Number of merged HVGs:", length(merged_vf)))


    options(future.globals.maxSize= 20000*1024^2) # 20gb max size
    
    #plan(strategy = "sequential")
    plan(strategy = "multicore", workers = 6)

    # Define the function to run for each donor_id_adj
    # TODO: parametrize what to iterate through?
    # TODO: parametrize args to RunPCA & RunHarmony?
    process_donor = function(donor_id) {
        print(paste("DONOR ID:", donor_id))
        target_seurat_subset = target_seurat_obj[, target_seurat_obj$donor_id_adj == donor_id]

        if (length(unique(target_seurat_subset$donor_id_adj)) > 1) {
            print("ERROR: target_seurat_subset BEFORE MERGE has more than one donor_id:")
            print(unique(target_seurat_subset$donor_id_adj))
            return(NULL)
        }
        if(ncol(target_seurat_subset) < MIN_NUM_CELLS){return(NULL)}

        merged_seurat_obj = mergeSeuratListWithMetadata(
            setNames(list(target_seurat_subset[merged_vf, ], source_seurat_obj[merged_vf, ]), 
            c(TARGET_DS_NAME, SOURCE_DS_NAME)),
            project="merged")
        print(paste(donor_id, ": Number of cells in source_seurat_obj:", ncol(source_seurat_obj)))
        print(paste(donor_id, "Number of cells in target_seurat_subset:", ncol(target_seurat_subset)))
        print(paste(donor_id, "Number of cells in merged_seurat_obj:", ncol(merged_seurat_obj)))
        merged_seurat_obj = (merged_seurat_obj 
            %>% NormalizeData()
            %>% ScaleData(verbose = FALSE, features = merged_vf, split.by="donor_id")
            %>% RunPCA(verbose = FALSE, features = merged_vf, split.by="donor_id", npcs = 50)
            %>% RunHarmony(
                    group.by.vars="dataset",
                    dims.use = 1:50,
                    early_stop = FALSE,
                    max_iter=10,
                    plot_convergence = FALSE
                )
        )

        source_seurat_subset = merged_seurat_obj[ , merged_seurat_obj$dataset == SOURCE_DS_NAME]
        target_seurat_subset = merged_seurat_obj[ , merged_seurat_obj$dataset == TARGET_DS_NAME]
        if (length(unique(target_seurat_subset$donor_id_adj)) != 1) {
            print("ERROR: target_seurat_subset AFTER MERGE has more than one donor_id:")
            print(unique(target_seurat_obj$donor_id_adj))
            return(NULL)
        }

        print(paste("RUNNING LABEL TRANSFER FOR", ncol(target_seurat_subset), "CELLS FOR DONOR", donor_id, "IN CELL CLASS", cc))

        lt = .myLabelTransfer_aligned(
            pca_source=source_seurat_subset@reductions$harmony@cell.embeddings,
            meta_source=source_seurat_subset@meta.data,
            pca_target=target_seurat_subset@reductions$harmony@cell.embeddings,
            meta_target=target_seurat_subset@meta.data,
            source_label_col=LT_COL,
            source_data=source_seurat_subset,
            target_data=target_seurat_subset,
            target_label_col=NULL,
            nPCs=50,
            n.adaptiveKernel=5,
            nPropIter=3,
            return_seurat_obj=TRUE)
        
        qsave(lt, file=file.path(BASE_DIR, LT_OUTPUT_DIR_NAME, paste0("lt_", cc, "_", donor_id, ".qs")))
    }

    print('RUNNING LABEL TRANSFER IN PARALLEL')
    # Apply the function in parallel for each donor_id_adj
    donor_ids = unique(target_seurat_obj$donor_id_adj)
    future_lapply(donor_ids, process_donor)
}