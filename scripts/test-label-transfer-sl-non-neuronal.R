
################# LIBRARIES #################

source("~/sc-online/utils.R")
source("~/sc-online/plot.R")
source("~/sc-online/clust.R")
source("~/sc-online/getData.R")
source("~/code/sconline_code.R")

library(dplyr)
library(ggplot2)
library(Matrix)
library(qs)
library(Seurat)
library(tidyr)
library(harmony)
library(stringr)

################# CONSTANTS #################


BASE_PATH="/mnt/accessory/seq_data/pd_all/240514"

astro = "astro"
endo = "endo"
mg = "mg"
da = "da"
nonda = "nonda"
oligo = "oligo"
fibro = "fibro"
opc = "opc"
mix = "mix"
mito = "mito"
none = "none"
immune = "immune"
ependymal = "ependymal"



################# FUNCTIONS #################

sum_duplicate_rownames_of_dgc_matrix=function(dgc){   
    # some umi_data has repeated rownames which causes problems down the road
    # we will sum the umi counts in these rows
    # collect rows with duplicate rownames as a dataframe
    n_occur = data.frame(table(rownames(dgc)))
    repeated_rownames = n_occur[n_occur$Freq > 1,]$Var1
    duplicated_dgc = as.data.frame(dgc[rownames(dgc) %in% repeated_rownames,])
    duplicated_dgc['symbol'] = sub("\\.\\d$", "", rownames(duplicated_dgc))

    # sum across numeric columns
    numeric_cols = setdiff(colnames(duplicated_dgc), "symbol")
    summed_duped_umi_rows = duplicated_dgc %>%
        group_by(symbol) %>%
        summarise_at(vars(numeric_cols), sum, na.rm = TRUE) %>%
        ungroup() %>%
        as.data.frame()
    # set rownames to be symbol column, then remove it
    rownames(summed_duped_umi_rows) = summed_duped_umi_rows$symbol
    summed_duped_umi_rows = summed_duped_umi_rows[, !(names(summed_duped_umi_rows) %in% c('symbol'))]

    # Hack
    # each round of the loop will remove ONE of each duplicate row 
    for(i in 1:max(n_occur$Freq)){
        dgc = dgc[!rownames(dgc) %in% repeated_rownames,]
    }

    # finally, replace the deduped rows with the summed rows via an rbind
    dgc = rbind(dgc, as(summed_duped_umi_rows, "sparseMatrix"))

    return(dgc)
}

.extraHumanGeneAnnoAdderFn=function(inputGeneNames=NULL){
  #require(EnsDb.Hsapiens.v75)
  require(EnsDb.Hsapiens.v86)
  
  if(!dir.exists("~/serverFiles")){
    dir.create("~/serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(EnsDb.Hsapiens.v86))
  gns$gene_short_name=gns$gene_name
  gns$symbol=toupper(gns$symbol)
  gns$ensembl_gene_id=row.names(gns)
  tst=unique(as.character(gns$gene_biotype))
  tst=setdiff(tst,c("protein_coding","lincRNA"))
  tst=factor(as.character(gns$gene_biotype),levels=c("protein_coding","lincRNA",tst))
  gns=gns[order(tst,decreasing = F),]
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    psCols=c("gene_short_name","ensembl_gene_id")
    slCounts=0
    slCol=""
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    {
      
      if(!file.exists("~/serverFiles/human_map_to_ensembl.rda")){
        stop("Error! file ~/serverFiles/human_map_to_ensembl.rda is missing")
      }
      
      load("~/serverFiles/human_map_to_ensembl.rda")
    }
    
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    if(!file.exists("~/serverFiles/human_mapping_hg19.rda")){
      stop("Error! file ~/serverFiles/human_mapping_hg19.rda is missing")
    }
    
    load("~/serverFiles/human_mapping_hg19.rda")
    human_hg19$source=toupper(human_hg19$source)
    
    if(sum(toupper(rwNames) %in% human_hg19$source) > sum(toupper(rwNames) %in% map_to_ensmbl$source)){
      map_to_ensmbl=rbind(human_hg19,map_to_ensmbl)
      map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    } else {
      map_to_ensmbl=rbind(map_to_ensmbl,human_hg19)
      map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    }
    
    c_map=data.frame(source=gns$gene_short_name,target=gns$ensembl_gene_id,stringsAsFactors = F)
    c_map=rbind(c_map,data.frame(source=gns$ensembl_gene_id,target=gns$ensembl_gene_id,stringsAsFactors = F))
    c_map=c_map[!is.na(c_map$ensembl_gene_id),]
    map_to_ensmbl=rbind(map_to_ensmbl,c_map)
    map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns$gene_id=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  return(gns)
}

parse.sconline.outs = function(outs, slot="combined_labels") {
  outs = outs[[slot]]
  col.names = outs %>% colnames()
  is.label.col = stringr::str_detect(colnames(outs), "inferred_")
  df_label_cols = outs[,is.label.col]
  colnames(df_label_cols) = stringr::str_remove(colnames(df_label_cols), "inferred_")
  
  pred.prob.mat = df_label_cols %>% as.matrix() %>% t()
  cell.names = pred.prob.mat %>% colnames()
  cell.names.removed.prefix = stringr::str_remove(cell.names, "^query_")
  colnames(pred.prob.mat) = cell.names.removed.prefix

  pred.prob.df = as.data.frame(t(pred.prob.mat))
  
  max_values = apply(pred.prob.df, 1, max)
  second_max_values = apply(pred.prob.df, 1, function(x) sort(x, decreasing = TRUE)[2]) 
  max_cols = colnames(pred.prob.df)[apply(pred.prob.df, 1, which.max)]

  pred.prob.df$max_lt_prob = max_values
  pred.prob.df$second_highest_lt_prob = second_max_values
  pred.prob.df$log2_first_to_second_lt_ratio = log2(pred.prob.df$max_lt_prob / (pred.prob.df$second_highest_lt_prob + 1E-6) + 1E-6)
  pred.prob.df$best_lt_label = max_cols

  return(pred.prob.df)
}




################# MAIN #################

# read in data
# gns = .extraHumanGeneAnnoAdderFn()
# dapi = qread(file.path(BASE_PATH, "dapi_seurat_clean.qs"))
# sl_non_neuronal_orig = readRDS("/mnt/accessory/seq_data/sl_ref/sl_non_neuronal.rds")


# # subset both objects / count matrices to the same genes
# dapi_gns = gns[gns$gene_name %in% rownames(dapi),]
# sl_counts = sl_non_neuronal_orig@assays$RNA@data
# sl_counts_in_dapi = sl_counts[rownames(sl_counts) %in% dapi_gns$gene_id,]

# # additionally, we want to remove rows of sl_counts_in_dapi with duplicated gene names
# dupe_gns = gns[duplicated(gns$gene_name),]
# sl_counts_in_dapi = sl_counts_in_dapi[!rownames(sl_counts_in_dapi) %in% dupe_gns$gene_id,]

# # now subset dapi to the same genes as sl_counts_in_dapi
# sl_gns = gns[gns$gene_id %in% rownames(sl_counts_in_dapi),]
# dapi_in_sl = dapi[rownames(dapi) %in% sl_gns$gene_name,]

# # rename rows of sl_counts to be gene names
# ens_to_name_map = setNames(
#     gns$gene_name,
#     gns$gene_id
# )
# rownames(sl_counts_in_dapi) = as.character(ens_to_name_map[rownames(sl_counts_in_dapi)])
# sl_non_neuronal = CreateSeuratObject(
#     counts = sl_counts_in_dapi,
#     project = "sl",
#     meta.data = sl_non_neuronal_orig@meta.data)

# sl_non_neuronal = normalizeScalePcaClusterUmap(
#     sl_non_neuronal,
#     var_feature_subset_col="donor_id", 
#     scaling_subset_col="donor_id",
#     n_hvgs_orig=2500, 
#     n_dims_use=20,
#     resolutions = c(0.2, 0.4, 0.6),
#     regression_vars = NULL)


# qsave(sl_non_neuronal, "/mnt/accessory/seq_data/sl_ref/sl_non_neuronal_in_dapi.qs")

print("Reading in data")
sl_non_neuronal = qread("/mnt/accessory/seq_data/sl_ref/sl_non_neuronal_in_dapi.qs")
sl_non_neuronal$participant_id = sl_non_neuronal$donor_id
sl_non_neuronal$cell_class = sl_non_neuronal$cell_type
hvgs = sl_non_neuronal@assays$RNA@var.features

print(paste("NUM HVGS:", length(hvgs)))
cat(head(hvgs))


dapi = qread(file.path(BASE_PATH, "dapi_seurat_clean.qs"))
dapi = dapi[rownames(dapi) %in% rownames(sl_non_neuronal), !dapi$cell_class %in% c("da", "nonda")]
dapi = dapi[match(rownames(sl_non_neuronal), rownames(dapi)),]

print("Scaling data")
dapi = ScaleData(dapi, features=hvgs, split.by='participant_id')

print("Merging datasets")

#merged = mergeSeuratListWithMetadata(list(sl_non_neuronal, dapi))

merged = merge(
    sl_non_neuronal[hvgs,],
    y=dapi[hvgs,],
    project = "merged"
)

ref_metadata = sl_non_neuronal@meta.data[, c("participant_id", "cell_class", "cluster_id", "subcluster_id")]
q_metadata = dapi@meta.data[, c("participant_id", "cell_class")]
q_metadata$cluster_id = "NA"
q_metadata$subcluster_id = "NA"

print("Concatenating meta.data")
merged@meta.data = rbind(
  ref_metadata, q_metadata
)
merged$dataset = ifelse(
    colnames(merged) %in% colnames(sl_non_neuronal),
    "reference",
    "query"
)

sl_scaled_data = sl_non_neuronal@assays$RNA@scale.data 
dapi_scaled_data = dapi@assays$RNA@scale.data

print("Concatenating scale.data")
merged@assays$RNA@scale.data = cbind(
    sl_scaled_data,
    dapi_scaled_data
)

print("Creating PCA")
feature_loadings = sl_non_neuronal[["pca"]]@feature.loadings
feature_loadings = feature_loadings[rownames(merged@assays$RNA@scale.data),]

merged_cell_embeddings = t(merged@assays$RNA@scale.data) %*% feature_loadings

merged[["pca"]] = CreateDimReducObject(
    embeddings = merged_cell_embeddings,
    loadings = feature_loadings,
    key="PC_", assay="RNA"
)

print("Running Harmony")
merged = (merged
    %>% RunHarmony(
        group.by.vars=c("dataset"),
        dims.use=1:20,
        plot_convergence=F,
        reference_values="reference"
    )
)

#print("Clustering")
#merged = (merged
    #%>% FindNeighbors(reduction="harmony")
    # %>% FindClusters(resolution=0.1)
    #%>% FindClusters(resolution=0.2)
    # %>% FindClusters(resolution=0.3)
    # %>% FindClusters(resolution=0.4)
    # %>% FindClusters(resolution=0.5)
    # %>% FindClusters(resolution=0.6)
    # %>% FindClusters(resolution=0.8)
    # %>% FindClusters(resolution=1)
    # %>% FindClusters(resolution=1.25)
    # %>% FindClusters(resolution=1.5)
    # %>% FindClusters(resolution=1.75)
    # %>% FindClusters(resolution=2)
    #%>% RunUMAP(reduction="harmony", dims=1:20)
#)

qsave(merged, "/mnt/accessory/seq_data/sl_ref/sl_non_neuronal_dapi_merged.qs")

print("Label Transfer")
query.embeddings = merged@reductions$harmony@cell.embeddings[merged$dataset == "query",]
ref.embeddings=merged@reductions$harmony@cell.embeddings[merged$dataset == "reference",]
query.metadata=merged@meta.data[merged$dataset == "query",]
ref.metadata=merged@meta.data[merged$dataset == "reference",]
target.col = "cluster_id" 
# Perform label transfer and convert output into probability matrix
sconline.res = .myLabelTransfer_aligned(pca_source = ref.embeddings, source_label_col = target.col, meta_source = ref.metadata,
                                        pca_target = query.embeddings, meta_target = query.metadata, return_seurat_obj = F)

print("Parsing Output")

prob_mat_combined = parse.sconline.outs(sconline.res)
merged@meta.data = cbind(merged@meta.data, prob_mat_combined)

qsave(merged, "/mnt/accessory/seq_data/sl_ref/sl_non_neuronal_dapi_merged.qs")