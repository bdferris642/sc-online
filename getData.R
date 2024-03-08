library(SingleCellExperiment)
library(Seurat)

source("~/sc-online/utils.R")

getSeuratVarFeatures = function(sobj){
    # the magical incantation that returns the slot of the attribute of the slot that actually holds the list of variable feature -_-
    return(sobj@assays$RNA@var.features)
}

getSeuratVarFeatureIntersectByCol = function(
    seurat_obj,
    subset_col,
    original_nfeatures,
    n_subsets_to_cover=NULL
){
    
    hvgs = list()

    unique_ids = unique(seurat_obj@meta.data[[subset_col]])

    if (is.null(n_subsets_to_cover)){
        n_subsets_to_cover = floor(length(unique_ids)/2)
    }

    i=1
    for (id in unique_ids) {
        print(paste("Subset", id, "--", i, "of", length(unique_ids)))
        i = i + 1
        
        seurat_subset = seurat_obj[, seurat_obj[[subset_col]] == id]
        
        if (ncol(seurat_subset) < 2){next}
        
        suppressWarnings({
            seurat_subset = FindVariableFeatures(
                seurat_subset, nfeatures=original_nfeatures, verbose = FALSE)
        })
        
        hvgs[[id]] = getSeuratVarFeatures(seurat_subset)
        
    }

    common_hvgs = getCommonStrings(hvgs, n_subsets_to_cover)
    print(paste("Number of HVGs in common across", n_subsets_to_cover, "--", length(common_hvgs)))

    return(common_hvgs)
}



harmonizeSeuratObjList = function(
    seurat_obj_list,
    hvgs=NULL
){
    # Given a list of seurat objects whose names correspond to different datasets...
    
    # (1) FindVariableFeatures for all, union across a subset of donors (at least half)
    # (2) ScaleData, grouping by "donor_id"
    # (3) RunPCA
    # (4) RunHarmony
    
    # Return the merged, harmonized Seurat objects


    # TODO fix this, it doesn't change the original dataset entries in the metadata
    dataset_names = names(seurat_obj_list)
    for (name in dataset_names) {
        seurat_obj = seurat_obj_list[[name]]
        seurat_obj@meta.data$dataset = name
    }

    # can skip the time-consuming hvg step if they have been precomputed and provided as an argument. 
    # Otherwise...
    if (is.null(hvgs)){
        dataset_hvgs = list()
        for (name in dataset_names) {
            seurat_obj = seurat_obj_list[[name]]
            seurat_obj$dataset = name

            dataset_hvgs[[name]] = getSeuratVarFeatureIntersectByCol(
                seurat_obj=seurat_obj,
                subset_col="donor_id",
                original_nfeatures=2500)
        }
        
        # Union dataset_hvgs
        hvgs = Reduce(union, dataset_hvgs)
    }

    # Combine
    seurat_merged = (
        mergeSeuratListWithMetadata(seurat_obj_list) %>%
        NormalizeData() %>%
        ScaleData(split.by = "donor_id", features=hvgs, verbose=F) %>% # assuming different donor_ids between datasets....
        RunPCA(verbose=F, npcs=50, features=hvgs)
    )

    # Harmonize
    options(repr.plot.width=7, repr.plot.heigh=7)
    seurat_merged = RunHarmony(
        object=seurat_merged, 
        group.by.vars="donor_id",
        dims.use = 1:50,
        early_stop = F,
        max_iter=20,
        plot_convergence = TRUE
    )
    
    # Return the merged, harmonized Seurat objects
    return(seurat_merged)
    
}


loadCbSceList = function(
    dir_list=calico_libs,
    filter_out_unclean=TRUE,
    filter_out_unassignable=TRUE,
    pct_mt_max=10,
    pct_intronic_min=20,
    log10_nUMI_threshold_list=NULL,
    short_lib_names=NULL,
    dapi_nurr=NULL,
    base_path="/mnt/accessory/seq_data/calico",
    cb_sce_basename="ingested_data/cb_data_sce_FPR_0.01.rds",
    vireo_donor_ids_basename="vireo_outs/donor_list/donor_ids.tsv",
    unclean_cells_basename = "ingested_data/unclean_cells.csv",
    n_donors_list = NULL,
    num_cells_cutoff = 15
    ){

    if (is.null(log10_nUMI_threshold_list)){
        log10_nUMI_threshold_list = list()
        for (name in dir_list) {
            log10_nUMI_threshold_list[[name]] = 3
        }
    }

    # load cb sce's, filter out unassignable cells
    cb_sce_list = list()
    for (name in dir_list) {

        cb_sce = readRDS(file.path(base_path, name, cb_sce_basename))
        cd = as.data.frame(colData(cb_sce))

        # add / update some columns 
        cd$log10_nUMI = log10(cd$nUMI + 1)    
        
        cd$QC_MT.pct[is.na(cd$QC_MT.pct)] = 0
        
        cd$pct_intronic = cd$pct_intronic * 100
        
        cd$prob_max = cd$prob_donor 
        
        cd$library = name
        
        if (is.null(short_lib_names)){
            cd$short_library = name
        } else {
            cd$short_library = as.character(short_lib_names[cd$library])
        }

        if (is.null(dapi_nurr)){
            cd$sort = NULL
        } else {
            cd$sort = dapi_nurr[[name]]
        }

        cd$short_donor_id = gsub(".*[-_]", "", cd$donor_id)

        log10_nUMI_threshold = log10_nUMI_threshold_list[[name]]
        cd$log10_nUMI_threshold = log10_nUMI_threshold
        
        
        # add is_clean and is_assignable columns
        cd$is_clean = 1
        cd$is_assignable = 1

        cd$is_clean[cd$log10_nUMI < log10_nUMI_threshold] = 0
        cd$is_clean[cd$pct_intronic < pct_intronic_min] = 0
        cd$is_clean[cd$QC_MT.pct >= pct_mt_max] = 0

        v = read.table(file.path(base_path, name, vireo_donor_ids_basename), header = TRUE, sep = "\t")
        v = v %>% filter(cell %in% colnames(cb_sce))
        
        v_assignable = v %>% filter(! donor_id %in% c("unassigned", "doublet"))

        if (!is.null(n_donors_list)){
            # get top donors, and subset to those with num assignable cells > num_cells_cutoff
            n_donors = n_donors_list[[name]]
            v_donor_counts = (v_assignable 
                %>% group_by(donor_id) 
                %>% summarize(n_cells = n()) 
                %>% filter(n_cells >= num_cells_cutoff)
                %>% arrange(desc(n_cells))
            )
            v_donor_counts = v_donor_counts[1:n_donors,]
            v_assignable = v_assignable %>% filter(donor_id %in% v_donor_counts$donor_id)
    
        }
        cd$is_assignable[! colnames(cb_sce) %in% v_assignable$cell] = 0
        cd$is_assignable_and_clean = cd$is_assignable & cd$is_clean

        # create new colnames of donor_id_barcode in the coldata
        # also add the pure barcode as a column
        cd = (cd %>% 
            mutate(
                barcode = colnames(cb_sce),
                donor_id_barcode = paste(donor_id, colnames(cb_sce), sep="_")
            )
        )
        rownames(cd) = cd$donor_id_barcode
        colnames(cb_sce) = cd$donor_id_barcode
        
        if(filter_out_unclean){
            #keep track of who you're filtering out
            cd_to_toss = cd[cd$is_clean==0,]
            write.csv(cd_to_toss, file.path(base_path, name, unclean_cells_basename))

            cb_sce=cb_sce[, cd$is_clean==1]
            cd = cd[cd$is_clean==1,]
        }
        if(filter_out_unassignable){
            cb_sce=cb_sce[, cd$is_assignable==1]
            cd = cd[cd$is_assignable==1,]
        }

        colData(cb_sce) = DataFrame(cd)
        cb_sce_list[[name]] = cb_sce
    }
    return(cb_sce_list)
}

mergeSeuratListWithMetadata = function(seurat_obj_list, cell_ids=NULL, project=NULL){

    # Harmonize metadata columns
    all_colnames = unique(unlist(lapply(seurat_obj_list, function(x) colnames(x@meta.data))))
    seurat_obj_list = lapply(seurat_obj_list, function(x) {
        missing_cols = setdiff(all_colnames, colnames(x@meta.data))
        if(length(missing_cols) > 0){
            x@meta.data[missing_cols] = NA
        }
        return(x)
    })

    if (is.null(project)){
        project = "merged"
    }
    
    if (is.null(cell_ids)){
        seurat_merged = Reduce(function(x, y) merge(x, y, project=project), 
            seurat_obj_list)
    } else {
        seurat_merged = Reduce(function(x, y) merge(x, y, project=project), 
            add.cell.ids=cell_ids, 
            seurat_obj_list)
    }
    
    md = lapply(seurat_obj_list, function(x){
        x@meta.data$orig.row.names = rownames(x@meta.data)
        x@meta.data
    })
    
    md = do.call(rbind, md)
    rownames(md) = md$orig.row.names
    md$orig.row.names = NULL
    seurat_merged@meta.data = md
    return(seurat_merged)
}




rawSceToHarmonizedSeurat = function(
    sce_list,
    n_donors_hvg = 12,
    split_col = "donor_id",
    n_var_features = 5000,
    n_dims_use = 50,
    res = 1.0,
    harmony_group_by_vars = "donor_id",
    project = NULL
){

    if (is.null(n_donors_hvg)){
        n_donors_hvg = floor(length(sce_list)/2)
    }

    # find variable genes within each donor
    donor_hvgs = list()
    donor_seurat_objs = list()
    
    for (donor_id in names(sce_list)) {
                
        donor_sce = sce_list[[donor_id]]

        if (ncol(donor_sce) < 2){next}

        donor_seurat = sceToSeurat(donor_sce, donor_id)
        donor_seurat = (
            Seurat::NormalizeData(object=donor_seurat[rowSums(donor_seurat) > 0, ]) %>%
            FindVariableFeatures(nfeatures = n_var_features)
        )
        donor_seurat_objs[[donor_id]] = donor_seurat
        donor_hvgs[[donor_id]] = donor_seurat@assays$RNA@var.features
    }

    # Get list of genes occuring in some number of donors within a sort.
    hvgs = getCommonStrings(donor_hvgs, n_donors_hvg)
    print(paste("Number of HVGs in common across", n_donors_hvg, "--", length(hvgs)))

    # Combine  
    seurat_merged = mergeSeuratListWithMetadata(donor_seurat_objs, project)
    
    # At the combined level
    # Subset to HVGs and Scale
    # remove second normalization. See:
    #https://docs.google.com/document/d/1p22HGFJK86-scFkPY817iHdnxurdWIaZ8HWBuXdJNxg/edit#bookmark=id.oefjnhz4ot00
    seurat_merged = (seurat_merged %>%
        #Seurat::NormalizeData() %>% 
        ScaleData(split.by = split_col, features = hvgs)
    )

    # (6) PCA
    seurat_merged = RunPCA(object=seurat_merged, features=hvgs, verbose=F, npcs=n_dims_use)


    # (7) Harmony (integrating donor)
    seurat_merged = RunHarmony(
        object=seurat_merged, 
        group.by.vars=harmony_group_by_vars,
        dims.use = 1:n_dims_use,
        early_stop = F,
        max_iter=25,
        plot_convergence = TRUE
    )

    # (8) Cluster
    # (9) UMAP
    seurat_merged = (seurat_merged %>%
        FindNeighbors(dims = 1:n_dims_use, reduction="harmony") %>%
        FindClusters(resolution = res) %>%
        RunUMAP(dims = 1:n_dims_use, verbose=F, reduction="harmony")
    )

    return(seurat_merged)
}

sceToSeurat = function(sce, project=NULL){
  sobj = CreateSeuratObject(counts = assays(sce)$counts, project=project)
  cd = as.data.frame(colData(sce))
  new_metadata = merge(sobj@meta.data, cd, by = "row.names")
  rownames(new_metadata) = new_metadata$Row.names
  new_metadata$Row.names = NULL
  sobj@meta.data = new_metadata 
  return(sobj)
}
