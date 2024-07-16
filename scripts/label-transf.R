
################# LIBRARIES #################

library(dplyr)
library(getopt)
library(ggplot2)
library(Matrix)
library(qs)
library(Seurat)
library(tidyr)
library(harmony)
library(stringr)

source("~/sc-online/utils.R")
source("~/sc-online/plot.R")
source("~/sc-online/clust.R")
source("~/sc-online/getData.R")
source("~/code/sconline_code.R")


################# CONSTANTS #################

print("READING IN ARGUMENTS")

spec <- matrix(c(
    'reference-base-path', 'rb', 1, "character",
    'reference-basename', 'rn', 1, "character",
    'query-base-path', 'qb', 1, "character",
    'query-basename', 'qn', 1, "character",
    'lt-column', 'lt', 1, "character",
    'merged-outname', 'mo', 1, "character",
    'var-feature-subset-col', 'vf', 1, "character",
    'scaling-subset-col', 'sc', 1, "character",
    'regression-vars', 'rv', 1, "character",
    'convert-gene-names', 'cg', 1, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)



REFERENCE_BASE_PATH = ifelse(
    is.null(opt[['reference-base-path']]), 
    "/mnt/accessory/seq_data/sl_ref", 
    opt[['reference-base-path']])
print(paste("REFERENCE_BASE_PATH:", REFERENCE_BASE_PATH))

REFERENCE_BASENAME = ifelse(
    is.null(opt[['reference-basename']]), 
    "sl_non_neuronal.qs", 
    opt[['reference-basename']])
print(paste("REFERENCE_BASENAME:", REFERENCE_BASENAME))

QUERY_BASE_PATH = ifelse(
    is.null(opt[['query-base-path']]), 
    "/mnt/accessory/seq_data/pd_all/240514", 
    opt[['query-base-path']])
print(paste("QUERY_BASE_PATH:", QUERY_BASE_PATH))

QUERY_BASENAME = ifelse(
    is.null(opt[['query-basename']]), 
    "dapi_seurat_clean.qs", 
    opt[['query-basename']]
)
print(paste("QUERY_BASENAME:", QUERY_BASENAME))

LT_COL = ifelse(
    is.null(opt[['lt-column']]), 
    "cluster_id", 
    opt[['lt-column']]
)
print(paste("LT_COL:", LT_COL))

MERGED_OUTNAME = ifelse(
    is.null(opt[['merged-outname']]), 
    paste0(
      "lt/label_transfer_prob_mat__query_",
      sub(".qs", "", QUERY_BASENAME), "__ref_",
       sub("\\..*", "", REFERENCE_BASENAME), "__", LT_COL, ".qs"),
    opt[['merged-outname']]
)
print(paste("MERGED_OUTNAME:", MERGED_OUTNAME))

if (is.null(opt[['convert-gene-names']])){
    CONVERT_GENE_NAMES = FALSE
} else {
    CONVERT_GENE_NAMES = opt[['convert-gene-names']]
}
print(paste("CONVERT_GENE_NAMES:", CONVERT_GENE_NAMES))

if (is.null(opt[['var-feature-subset-col']])){
    VAR_FEATURE_SUBSET_COL = NULL
} else {
    VAR_FEATURE_SUBSET_COL = opt[['var-feature-subset-col']]
}
print(paste("VAR_FEATURE_SUBSET_COL:", VAR_FEATURE_SUBSET_COL))

if (is.null(opt[['scaling-subset-col']])){
    SCALING_SUBSET_COL = NULL
} else {
    SCALING_SUBSET_COL = opt[['scaling-subset-col']]
}
print(paste("SCALING_SUBSET_COL:", SCALING_SUBSET_COL))

if (is.null(opt[['regression-vars']])) {
    REGRESSION_VARS = NULL
} else {
    REGRESSION_VARS = strsplit(opt[['regression-vars']], ",")[[1]]
}

print("REGRESSION_VARS:")
cat(REGRESSION_VARS)

################# FUNCTIONS #################

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

query = qread(file.path(QUERY_BASE_PATH, QUERY_BASENAME))
if (grepl(".qs$", REFERENCE_BASENAME)) {
    ref = qread(file.path(REFERENCE_BASE_PATH, REFERENCE_BASENAME))
} else if (grepl(".rds$", REFERENCE_BASENAME)){
    ref = readRDS(file.path(REFERENCE_BASE_PATH, REFERENCE_BASENAME))
} else {
    stop("Reference file must be either .qs or .rds")
}

if (CONVERT_GENE_NAMES){
    print("Converting gene names")
    
    #read in data
    gns = .extraHumanGeneAnnoAdderFn()

    # subset both objects / count matrices to the same genes
    # need to get references from ensemble ids to gene names
    # todo: parametrize this (some references will have gene names as row names)
    query_gns = gns[gns$gene_name %in% rownames(query),]
    ref_counts = ref@assays$RNA@data
    ref_counts_in_query = ref_counts[rownames(ref_counts) %in% query_gns$gene_id,]

    # additionally, we want to remove rows of ref_counts_in_query with duplicated gene names
    dupe_gns = gns[duplicated(gns$gene_name),]
    ref_counts_in_query = ref_counts_in_query[!rownames(ref_counts_in_query) %in% dupe_gns$gene_id,]

    # now subset the query to the same genes as ref_counts_in_query
    ref_gns = gns[gns$gene_id %in% rownames(ref_counts_in_query),]
    query_in_ref = query[rownames(query) %in% ref_gns$gene_name,]

    # rename rows of ref counts matrix to be gene names
    ens_to_name_map = setNames(
        gns$gene_name,
        gns$gene_id
    )

    rownames(ref_counts_in_query) = as.character(ens_to_name_map[rownames(ref_counts_in_query)])

    ref_renamed = CreateSeuratObject(
        counts = ref_counts_in_query,
        project = "sl",
        meta.data = ref@meta.data)
    
} else {
    
    joint_rownames = sort(intersect(rownames(query), rownames(ref)))
    query = query[joint_rownames,]
    ref_renamed = ref[joint_rownames,]

}

if ((! "cell_class" %in% colnames(ref_renamed@meta.data)) & 
    ("cell_type" %in% colnames(ref_renamed@meta.data))) {
        ref_renamed$cell_class = ref_renamed$cell_type
    }

if ((! "participant_id" %in% colnames(ref_renamed@meta.data)) & 
    ("donor_id" %in% colnames(ref_renamed@meta.data))) {
        ref_renamed$participant_id = ref_renamed$donor_id
    }

if ((! "pct_mito" %in% colnames(ref_renamed@meta.data)) & 
    ("fraction_mitochondrial" %in% colnames(ref_renamed@meta.data))) {
        ref_renamed$pct_mito = 100*ref_renamed$fraction_mitochondrial
    }

if ((! "pct_intronic" %in% colnames(ref_renamed@meta.data)) & 
    ("fraction_unspliced" %in% colnames(ref_renamed@meta.data))) {
        ref_renamed$pct_intronic = 100*ref_renamed$fraction_unspliced
    }


if ((! "cell_class" %in% colnames(query@meta.data)) & 
    ("cell_type" %in% colnames(query@meta.data))) {
        query$cell_class = query$cell_type
    }

if ((! "participant_id" %in% colnames(query@meta.data)) & 
    ("donor_id" %in% colnames(query@meta.data))) {
        query$participant_id = query$donor_id
    }


ref_renamed = normalizeScalePcaClusterUmap(
    ref_renamed,
    var_feature_subset_col=VAR_FEATURE_SUBSET_COL,
    scaling_subset_col=SCALING_SUBSET_COL,
    n_hvgs_orig=2500, 
    n_dims_use=20,
    resolutions = c(0.2),
    regression_vars = REGRESSION_VARS)
hvgs = ref_renamed@assays$RNA@var.features

print(paste("NUM HVGS:", length(hvgs)))
cat(head(hvgs), "...")

query_counts = query@assays$RNA@data
query_counts = query_counts[rownames(query_counts) %in% rownames(ref_renamed),]
query_counts = query_counts[match(rownames(ref_renamed), rownames(query_counts)),]

query_renamed = CreateSeuratObject(
    counts = query_counts,
    project = "query",
    meta.data = query@meta.data)

print("Scaling data")
query_renamed = query_renamed %>% NormalizeData() %>% ScaleData(features=hvgs, split.by=SCALING_SUBSET_COL, vars.to.regress=REGRESSION_VARS)

print("Merging datasets")

merged = merge(
    ref_renamed[hvgs,],
    y=query_renamed[hvgs,],
    project = "merged"
)

ref_metadata = ref_renamed@meta.data[, c("participant_id", "cell_class", LT_COL)]
q_metadata = query_renamed@meta.data[, c("participant_id", "cell_class")]
q_metadata[[LT_COL]] = "NA"

print("Concatenating meta.data")
merged@meta.data = rbind(
  ref_metadata, q_metadata
)
merged$dataset = ifelse(
    colnames(merged) %in% colnames(ref_renamed),
    "reference",
    "query"
)

ref_scaled_data = ref_renamed@assays$RNA@scale.data 
query_scaled_data = query_renamed@assays$RNA@scale.data

print("Concatenating scale.data")
merged@assays$RNA@scale.data = cbind(
    ref_scaled_data,
    query_scaled_data
)

print("Creating PCA")
feature_loadings = ref_renamed[["pca"]]@feature.loadings
feature_loadings = feature_loadings[rownames(merged@assays$RNA@scale.data),]

merged_cell_embeddings = t(merged@assays$RNA@scale.data) %*% feature_loadings

merged[["pca"]] = CreateDimReducObject(
    embeddings = merged_cell_embeddings,
    loadings = feature_loadings,
    key="PC_", assay="RNA"
)

print("Running Harmony")
# todo: parametrize the group.by.vars?
merged = (merged
    %>% RunHarmony(
        group.by.vars=c("dataset"),
        dims.use=1:20,
        plot_convergence=F,
        reference_values="reference"
    )
)

#qsave(merged, file.path(REFERENCE_BASE_PATH, MERGED_OUTNAME))

print("Label Transfer")
query.embeddings = merged@reductions$harmony@cell.embeddings[merged$dataset == "query",]
ref.embeddings=merged@reductions$harmony@cell.embeddings[merged$dataset == "reference",]
query.metadata=merged@meta.data[merged$dataset == "query",]
ref.metadata=merged@meta.data[merged$dataset == "reference",]

# Perform label transfer and convert output into probability matrix
sconline.res = .myLabelTransfer_aligned(pca_source = ref.embeddings, source_label_col = LT_COL, meta_source = ref.metadata,
                                        pca_target = query.embeddings, meta_target = query.metadata, return_seurat_obj = F)

print("Parsing Output")

prob_mat_combined = parse.sconline.outs(sconline.res)
merged@meta.data = cbind(merged@meta.data, prob_mat_combined)

prob_mat_out = merged@meta.data[,
    c("dataset", "participant_id", "cell_class", LT_COL, colnames(prob_mat_combined))
]

prob_mat_out_query = prob_mat_out[prob_mat_out$dataset == "query", colnames(prob_mat_combined)]
prob_mat_out_query = prob_mat_out_query[match(colnames(query_renamed), rownames(prob_mat_out_query)),]

# paste inferred_<LT_COL>_ to the column names
new_colnames = paste0("inferred_", LT_COL, "__", colnames(prob_mat_out_query))
colnames(prob_mat_out_query) = new_colnames

# append this to the actual query object
query_renamed@meta.data = cbind(query_renamed@meta.data, prob_mat_out_query)
qsave(query_renamed, file.path(QUERY_BASE_PATH, QUERY_BASENAME))

qsave(prob_mat_out, file.path(QUERY_BASE_PATH, MERGED_OUTNAME))