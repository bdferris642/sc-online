############## LIBRARIES ##############

source("~/sc-online/extraGenes.R")
source("~/sc-online/myExpSetCreatorFn.R")
source("~/sc-online/extraGenes.R")
source("~/sc-online/myExpSetCreatorFn.R")
source("~/sc-online/utils.R")
source("~/sc-online/plot.R")
source("~/sc-online/clust.R")
source("~/sc-online/labelTransfer.R")
source("~/sc-online/getData.R")

library(caret)
library("DropletUtils")
library(dplyr)
library(getopt)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(harmony)
library(lisi)
library(Matrix)
library(patchwork)
library(pheatmap)
library(qs)
library(RColorBrewer)
library(SingleCellExperiment)
library(rhdf5)
library(rlang)

library(Seurat)
library(tidyr)
library(xml2)

library(viridis)
library(viridisLite)

############## ARGUMENTS ##############

spec = matrix(c(
    'write-slogan', 's', 1, "character",
    'res', 'r', 1, "numeric",
    'base-path', 'b', 1, "character",
    'cb-sce-basename', 'cb', 1, "character",
    'pct-intronic-min', 'pi', 1, "numeric",
    'pct-mt-max', 'pm', 1, "numeric",
    'summary-basename', 'sb', 1, "character",
    'calico-libs', 'cl', 1, "character",
    'gtex-libs', 'gl', 1, "character",
    'var-adj-pca', 'vp', 1, "logical",
    'harmonize', 'h', '1', 'logical',
    'lib-info-path', 'li', 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

WRITE_SLOGAN = opt[['write-slogan']]
if (is.null(WRITE_SLOGAN)) {
    stop("MUST PROVIDE AN OUTPUT WRITE SLOGAN")
}

BASE_PATH = ifelse(
    is.null(opt[['base-path']]), 
    "/mnt/accessory/seq_data", 
    opt[['base-path']])

CB_SCE_BASENAME = ifelse(
    is.null(opt[['cb-sce-basename']]), 
    "ingested_data_no_subset/cb_data_sce_FPR_0.01.rds", 
    opt[['cb-sce-basename']])

PCT_INTRONIC_MIN = ifelse(
    is.null(opt[['pct-intronic-min']]), 
    20, 
    opt[['pct-intronic-min']])

PCT_MT_MAX = ifelse(
    is.null(opt[['pct-mt-max']]), 
    10, 
    opt[['pct-mt-max']])

SUMMARY_BASENAME = ifelse(
    is.null(opt[['summary-basename']]),
    'vireo_outs/no_subset/summary.tsv',
    opt[['summary-basename']]
)

CALICO_LIBS_LONG_PATH = ifelse(
    is.null(opt[['calico-libs']]),
    "~/calico-libs-long.txt",
    opt[['calico-libs']]
)

GTEX_LIBS_LONG_PATH = ifelse(
    is.null(opt[['gtex-libs']]),
    "~/gtex-libs-long.txt",
    opt[['gtex-libs']]
)

HARMONIZE = ifelse(
    is.null(opt[['harmonize']]),
    TRUE,
    opt[['harmonize']]
)

VAR_ADJ_PCA = ifelse(
    is.null(opt[['var-adj-pca']]),
    FALSE,
    opt[['var-adj-pca']]
)

LIB_INFO_PATH = ifelse(
    is.null(opt[['lib-info-path']]),
    "~/sc-online/notebook_data/pd_lib_info_20240510.csv",
    opt[['lib-info-path']]
)


############## CONSTANTS ##############

OUTPUT_DIR = "pd_all"
dir.create(file.path(BASE_PATH, OUTPUT_DIR, WRITE_SLOGAN), showWarnings = FALSE)

write_basename_dapi = file.path(OUTPUT_DIR, WRITE_SLOGAN, "cb_sce_donor_list_dapi.qs")
write_basename_dapi_harmonized = file.path(OUTPUT_DIR, WRITE_SLOGAN, "seurat_dapi_merged.qs")

write_basename_nurr = file.path(OUTPUT_DIR, WRITE_SLOGAN, "cb_sce_donor_list_nurr.qs") 
write_basename_nurr_harmonized = file.path(OUTPUT_DIR, WRITE_SLOGAN, "seurat_nurr_merged.qs") 

EARLY_LIBS = readLines("~/early-libs.txt")

############## GRAB LIBRARY LISTS ##############

print("Loading library lists")

calico_libs_long = readLines(CALICO_LIBS_LONG_PATH)
calico_libs_long = calico_libs_long[calico_libs_long != ""]
calico_libs = lapply(calico_libs_long, function(x) {
    split = strsplit(x, split = "_")[[1]]
    return(paste(split[2:length(split)], collapse = "_"))
})
names(calico_libs_long) = calico_libs

gtex_libs_long = readLines(GTEX_LIBS_LONG_PATH)
gtex_libs_long = gtex_libs_long[gtex_libs_long != ""]
gtex_libs = lapply(gtex_libs_long, function(x) {
    split = strsplit(x, split = "_")[[1]]
    return(paste(split[2:length(split)], collapse = "_"))
})
names(gtex_libs_long) = gtex_libs

libs_long = c(calico_libs_long, gtex_libs_long)
libs = c(calico_libs, gtex_libs)

lib_info = read.table(LIB_INFO_PATH, sep = ",", header = TRUE )
lib_info_calico = lib_info[lib_info$source == "Calico",]
lib_info_gtex = lib_info[lib_info$source == "GTEx",]

dapi_nurr_calico = setNames(lib_info_calico$sort, lib_info_calico$library)
dapi_nurr_gtex = setNames(lib_info_gtex$sort, lib_info_gtex$library)

# TODO can change this later
log10_nUMI_thresholds_calico = rep(2.5, length(calico_libs))
names(log10_nUMI_thresholds_calico) = calico_libs

log10_nUMI_thresholds_gtex = rep(2.5, length(gtex_libs))
names(log10_nUMI_thresholds_gtex) = gtex_libs

############## LOAD CALICO SCEs ##############

print("Loading calico SCEs")

libs2skip = c()
for (i in c(calico_libs, gtex_libs)){
    # if the file doesn't exist, append the lib to the skip list
    if(i %in% calico_libs){
        base_path = "/mnt/accessory/seq_data/calico"
    } else {
        base_path = "/mnt/accessory/seq_data/gtex"
    }
    if( (!file.exists(file.path(base_path, i, "vireo_outs/no_subset/donor_ids.tsv"))) 
        | (!file.exists(file.path(base_path, i, "ingested_data_no_subset/cb_data_sce_FPR_0.01.rds")))){
        libs2skip = c(libs2skip, i)
    }
}

n_donors_calico = list()
for (lib in calico_libs[!calico_libs %in% libs2skip]) {
    donor_list = readLines(file.path(BASE_PATH, 'calico', lib, "vireo_outs/donor_list/donor_list.txt"))
    donor_list = donor_list[donor_list != ""]
    n_donors_calico[[lib]] = length(donor_list)
}

# # load calico sce lists
# # TODO: parametrize cb sce basename
# calico_sce_list = loadCbSceList(calico_libs[!calico_libs %in% libs2skip],
#     base_path="/mnt/accessory/seq_data/calico",
#     vireo_donor_ids_basename="vireo_outs/no_subset/donor_ids.tsv",
#     cb_sce_basename="ingested_data_no_subset/cb_data_sce_FPR_0.01.rds",
#     pct_mt_max=PCT_MT_MAX,
#     pct_intronic_min= PCT_INTRONIC_MIN,
#     log10_nUMI_threshold_list=log10_nUMI_thresholds_calico,
#     dapi_nurr=dapi_nurr_calico,
#     n_donors_list = n_donors_calico
# )

# print("Num Calico Donors:")
# print(length(calico_sce_list))

# ############## LOAD GTEX SCES ##############

# print("Loading gtex SCEs")

# # TODO: parametrize cb sce basename
# gtex_sce_list = loadCbSceList(gtex_libs[!gtex_libs %in% libs2skip],
#     base_path="/mnt/accessory/seq_data/gtex",
#     vireo_donor_ids_basename="vireo_outs/no_subset/donor_ids.tsv",
#     cb_sce_basename="ingested_data_no_subset/cb_data_sce_FPR_0.01.rds",
#     pct_mt_max=PCT_MT_MAX,
#     pct_intronic_min= PCT_INTRONIC_MIN,
#     log10_nUMI_threshold_list=log10_nUMI_thresholds_gtex,
#     dapi_nurr = dapi_nurr_gtex,
#     n_donors_list = NULL
# )

# print("Num GTEx Donors:")
# print(length(gtex_sce_list))

# ############## COMBINE LIBRARIES, SPLIT BY SORT, SCALE BY PARTICIPANT, AND HARMONIZE ##############
# print("Combining libraries")

# # join calico and gtex. 
# calico_all = .mycBindFn(calico_sce_list)
# gtex_all = .mycBindFn(gtex_sce_list)
# # first, we want to append our best understanding of participant IDs to each of the SCEs

# sces_all = .mycBindFn(list(calico_all, gtex_all))
# print(sort(unique(colData(sces_all)$participant_id)))

# #we would like to append metadata about library_prep to each library 
# sces_all$early = sces_all$library %in% EARLY_LIBS

# print("Splitting by sort and participant")

# # Split by sort and then by donor. Save.
# sce_nurr = sces_all[, colData(sces_all)$sort == 'nurr']
# sce_dapi = sces_all[, colData(sces_all)$sort == 'dapi']
# print(sort(unique(colData(sce_dapi)$participant_id)))
# print(sort(unique(colData(sce_nurr)$participant_id)))

# cb_sce_nurr_donor_list = .mySplitObject(sce_nurr, "participant_id")
# cb_sce_nurr_donor_list = cb_sce_nurr_donor_list[!is.na(names(cb_sce_nurr_donor_list))]

# cb_sce_dapi_donor_list = .mySplitObject(sce_dapi, "participant_id")
# cb_sce_dapi_donor_list = cb_sce_dapi_donor_list[!is.na(names(cb_sce_dapi_donor_list))]

# print("Saving SCE lists")

# print("Num DAPI Participants")
# print(length(cb_sce_dapi_donor_list))
# print(sort(names(cb_sce_dapi_donor_list)))

# print("Num NURR Participants")
# print(length(cb_sce_nurr_donor_list))
# print(sort(names(cb_sce_nurr_donor_list)))

# qsave(cb_sce_nurr_donor_list, file.path(BASE_PATH, write_basename_nurr))
# qsave(cb_sce_dapi_donor_list, file.path(BASE_PATH, write_basename_dapi))

cb_sce_nurr_donor_list = qread(file.path(BASE_PATH, write_basename_nurr))
cb_sce_dapi_donor_list = qread(file.path(BASE_PATH, write_basename_dapi))


# Harmonize. Takes 90 mins for each!
print("Generating harmonized NURR Seurat")
seurat_nurr_merged = rawSceToHarmonizedSeurat(
    cb_sce_nurr_donor_list,
    split_col = "participant_id",
    harmony_group_by_vars = "participant_id",
    n_donors_hvg = NULL, # have to be in half of the participants
    n_var_features = 2500,
    n_dims_use = 20,
    harmonize=HARMONIZE,
    var_adj_pca = VAR_ADJ_PCA)
print("Saving harmonized NURR Seurat")
qsave(seurat_nurr_merged, file.path(BASE_PATH, write_basename_nurr_harmonized))

print("Generating harmonized DAPI Seurat")
seurat_dapi_merged = rawSceToHarmonizedSeurat(
    cb_sce_dapi_donor_list,
    split_col = "participant_id",
    harmony_group_by_vars = "participant_id",
    n_donors_hvg = NULL, # have to be in half of the participants
    n_var_features = 2500,
    n_dims_use = 20,
    harmonize=HARMONIZE,
    var_adj_pca = VAR_ADJ_PCA)

print("Saving harmonized DAPI Seurat")
qsave(seurat_dapi_merged, file.path(BASE_PATH, write_basename_dapi_harmonized))
