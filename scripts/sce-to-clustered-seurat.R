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
    'lib-info-path', 'li', 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

WRITE_SLOGAN = opt[['write-slogan']]
if (is.null(WRITE_SLOGAN)) {
    stop("MUST PROVIDE AN OUTPUT WRITE SLOGAN")
}

RES = ifelse(
    is.null(opt[['res']]), 
    0.5, 
    opt[['res']])
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
    'vireo_outs/donor_list/summary.tsv',
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

VAR_ADJ_PCA = ifelse(
    is.null(opt[['var-adj-pca']]),
    FALSE,
    opt[['var-adj-pca']]
)

LIB_INFO_PATH = ifelse(
    is.null(opt[['lib-info-path']]),
    "~/sc-online/notebook_data/pd_lib_info_20240402.csv",
    opt[['lib-info-path']]
)


############## CONSTANTS ##############

OUTPUT_DIR = "pd_all"
dir.create(file.path(BASE_PATH, OUTPUT_DIR, WRITE_SLOGAN), showWarnings = FALSE)

write_basename_dapi = file.path(OUTPUT_DIR, WRITE_SLOGAN, "cb_sce_donor_list_dapi.qs")
write_basename_dapi_harmonized = file.path(OUTPUT_DIR, WRITE_SLOGAN, "seurat_dapi_harmonized.qs")

write_basename_nurr = file.path(OUTPUT_DIR, WRITE_SLOGAN, "cb_sce_donor_list_nurr.qs") 
write_basename_nurr_harmonized = file.path(OUTPUT_DIR, WRITE_SLOGAN, "seurat_nurr_harmonized.qs") 

EARLY_LIBS = list(
    "pPDsHSrSNiPoold221126A1"
    ,"pPDsHSrSNiPoold221126B1"
    ,"pPDsHSrSNiPoold221126C1"
    ,"pPDsHSrSNiPoold221126D1"
    ,"pPDsHSrSNiPoold221126E1"
    ,"pPDsHSrVTAiPoold221126F1"
    ,"pPDsHSrVTAiPoold221126G1"
    ,"pPDsHSrVTAiPoold221126H1"
    ,"pPDsHSrVTAiPoold221126A2"
    ,"pPDsHSrVTAiPoold221126B2"
    ,"pPDsHSrSNiPoold221219A3"
    ,"pPDsHSrSNiPoold221219B3"
    ,"pPDsHSrVTAiPoold221219C3"
    ,"pPDsHSrVTAiPoold221219D3"
    ,"pPDsHSrSNiPoold221219E3"
    ,"pPDsHSrSNiPoold221219F3"
    ,"pPDsHSrSNiPoold221219G3"
    ,"pPDsHSrVTAiPoold221219H3"
    ,"pPDsHSrVTAiPoold221219A4"
    ,"pPDsHSrVTAiPoold221219B4"
    ,"pCalicoPDsHSrSNDAPI1iPoold230214"
    ,"pCalicoPDsHSrSNDAPI2iPoold230214"
    ,"pCalicoPDsHSrVTADAPI1iPoold230214"
    ,"pCalicoPDsHSrVTADAPI2iPoold230214"
    ,"pCalicoPDsHSrRRFDAPI1iPoold230214"
    ,"pCalicoPDsHSrRRFDAPI2iPoold230214"
    ,"pCalicoPDsHSrRRFNURR1iPoold230214"
    ,"pCalicoPDsHSrRRFNURR2iPoold230214"
    ,"pCalicoPDsHSrRRFNURR3iPoold230214"
    ,"pCalicoPDsHSrVTANURR1iPoold230214"
    ,"pCalicoPDsHSrVTANURR2iPoold230214"
    ,'pCalicoPDsHSrVTANURR3iPoold230214'
    ,"pCalicoPDsHSrSNNURR1iPoold230214"
    ,"pCalicoPDsHSrSNNURR2iPoold230214"
    ,"pCalicoPDsHSrSNNURR3iPoold230214"
)


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

n_donors_calico = list()
for (lib in calico_libs) {
    donor_list = readLines(file.path(BASE_PATH, 'calico', lib, "vireo_outs/donor_list/donor_list.txt"))
    donor_list = donor_list[donor_list != ""]
    n_donors_calico[[lib]] = length(donor_list)
}

# load calico sce lists
# TODO: parametrize cb sce basename
calico_sce_list = loadCbSceList(calico_libs,
    base_path="/mnt/accessory/seq_data/calico",
    vireo_donor_ids_basename="vireo_outs/no_subset/donor_ids.tsv",
    cb_sce_basename="ingested_data_no_subset/cb_data_sce_FPR_0.01.rds",
    pct_mt_max=PCT_MT_MAX,
    pct_intronic_min= PCT_INTRONIC_MIN,
    log10_nUMI_threshold_list=log10_nUMI_thresholds_calico,
    dapi_nurr=dapi_nurr_calico,
    n_donors_list = n_donors_calico
)

############## LOAD GTEX SCES ##############

print("Loading gtex SCEs")

# TODO: remove this "skipped" block once all gtex libraries have been processed
skipped = c()
for (d in gtex_libs){
    if (d == ""){next}
    donor_ids_path = file.path(BASE_PATH, "gtex", d, "vireo_outs", "donor_list", "donor_ids.tsv")
    if (!file.exists(donor_ids_path)){
        skipped = c(skipped, d)
    }
}
# these libraries are bad and should probably jsut be removed from the gtex list
skipped = c(skipped
    , "pCalico_GTExsHSrSNA11iNURRd231120"
    ,"pCalico_GTExsHSrSNB11iNURRd231120"
    ,"pCalico_GTExsHSrSNC11iDAPId231120"
    ,"pCalico_GTExsHSrSND11iNURRd231120"
    ,"pCalico_GTExsHSrSNE11iNURRd231120"
    ,"pCalico_GTExsHSrSNF11iDAPId231120"
    ,""
)

# TODO: parametrize cb sce basename
gtex_sce_list = loadCbSceList(gtex_libs[!gtex_libs %in% skipped],
    base_path="/mnt/accessory/seq_data/gtex",
    vireo_donor_ids_basename="vireo_outs/donor_list/donor_ids.tsv",
    cb_sce_basename="ingested_data/cb_data_sce_FPR_0.01.rds",
    pct_mt_max=PCT_MT_MAX,
    pct_intronic_min= PCT_INTRONIC_MIN,
    log10_nUMI_threshold_list=log10_nUMI_thresholds_gtex,
    dapi_nurr = dapi_nurr_gtex,
    n_donors_list = NULL
)

# TODO! remove this block for non SN/VTA
# in select libraries (pCalico_GTExsHSrSNA11iNURRd231120, pCalico_GTExsHSrSNB11iNURRd231120, and pCalico_GTExsHSrSNC11iDAPId231120)
# we want to remove unwanted donors (GTEX-QDT8, GTEX-RVPU, GTEX-ZAJG, GTEX-R55F, GTEX-RNOR)
for (name in names(gtex_sce_list)){
    if (name %in% c("pCalico_GTExsHSrSNA11iNURRd231120", "pCalico_GTExsHSrSNB11iNURRd231120", "pCalico_GTExsHSrSNC11iDAPId231120")){
        sce = gtex_sce_list[[name]]
        sce = sce[, !sce$donor_id %in% c(
            "GTEX-QDT8", "GTEX-RVPU", "GTEX-ZAJG", "GTEX-R55F", "GTEX-RNOR"
        )]
        gtex_sce_list[[name]] = sce
    }
}

############## COMBINE LIBRARIES, SPLIT BY SORT, SCALE BY PARTICIPANT, AND HARMONIZE ##############
print("Combining libraries")

# join calico and gtex. 
calico_all = .mycBindFn(calico_sce_list)
gtex_all = .mycBindFn(gtex_sce_list)
# first, we want to append our best understanding of participant IDs to each of the SCEs

sces_all = .mycBindFn(list(calico_all, gtex_all))

#we would like to append metadata about library_prep to each library 
sces_all$early = sces_all$library %in% EARLY_LIBS

print("Splitting by sort and participant")

# Split by sort and then by donor. Save.
sce_nurr = sces_all[, colData(sces_all)$sort == 'nurr']
sce_dapi = sces_all[, colData(sces_all)$sort == 'dapi']

cb_sce_nurr_donor_list = .mySplitObject(sce_nurr, "participant_id")
cb_sce_nurr_donor_list = cb_sce_nurr_donor_list[!is.na(names(cb_sce_nurr_donor_list))]

cb_sce_dapi_donor_list = .mySplitObject(sce_dapi, "participant_id")
cb_sce_dapi_donor_list = cb_sce_dapi_donor_list[!is.na(names(cb_sce_dapi_donor_list))]

print("Saving SCE lists")

qsave(cb_sce_nurr_donor_list, file.path(BASE_PATH, write_basename_nurr))
qsave(cb_sce_dapi_donor_list, file.path(BASE_PATH, write_basename_dapi))

# Harmonize. Takes 90 mins for each!
print("Generating harmonized NURR Seurat")
seurat_nurr_merged = rawSceToHarmonizedSeurat(
    cb_sce_nurr_donor_list,
    split_col = "participant_id",
    harmony_group_by_vars = "participant_id",
    n_donors_hvg = NULL, # have to be in half of the participants
    n_var_features = 5000,
    n_dims_use = 50,
    res = RES,
    var_adj_pca = VAR_ADJ_PCA)
print("Saving harmonized NURR Seurat")
qsave(seurat_nurr_merged, file.path(BASE_PATH, write_basename_nurr_harmonized))
seurat_nurr_merged$initial_clusters = seurat_nurr_merged$seurat_clusters

print("Generating harmonized DAPI Seurat")
seurat_dapi_merged = rawSceToHarmonizedSeurat(
    cb_sce_dapi_donor_list,
    split_col = "participant_id",
    harmony_group_by_vars = "participant_id",
    n_donors_hvg = NULL, # have to be in half of the participants
    n_var_features = 5000,
    n_dims_use = 50,
    res = RES,
    var_adj_pca = VAR_ADJ_PCA)
seurat_dapi_merged$initial_clusters = seurat_dapi_merged$seurat_clusters

print("Saving harmonized DAPI Seurat")
qsave(seurat_dapi_merged, file.path(BASE_PATH, write_basename_dapi_harmonized))
