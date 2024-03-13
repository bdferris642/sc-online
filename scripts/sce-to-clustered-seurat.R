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

############## CONSTANTS ##############

BASE_PATH = "/mnt/accessory/seq_data"

VIREO_DONOR_IDS_BASENAME = 'vireo_outs/no_subset/donor_ids.tsv'
CB_SCE_BASENAME = "ingested_data_no_subset/cb_data_sce_FPR_0.01.rds"
SUMMARY_BASENAME = 'vireo_outs/donor_list/summary.tsv'

# Set relatively conservative thresholds for initial filtering
PCT_INTRONIC_MIN = 20
PCT_MT_MAX = 10

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

calico_libs_long = readLines("~/calico-libs-long.txt")
calico_libs_long = calico_libs_long[calico_libs_long != ""]
calico_libs = lapply(calico_libs_long, function(x) {
    split = strsplit(x, split = "_")[[1]]
    return(paste(split[2:length(split)], collapse = "_"))
})
names(calico_libs_long) = calico_libs

gtex_libs_long = readLines("~/gtex-libs-long.txt")
gtex_libs_long = gtex_libs_long[gtex_libs_long != ""]
gtex_libs = lapply(gtex_libs_long, function(x) {
    split = strsplit(x, split = "_")[[1]]
    return(paste(split[2:length(split)], collapse = "_"))
})
names(gtex_libs_long) = gtex_libs

libs_long = c(calico_libs_long, gtex_libs_long)
libs = c(calico_libs, gtex_libs)

lib_info = read.table("~/sc-online/notebook_data/pd/pd_lib_info_20240301.tsv", sep = "\t", header = TRUE )
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
calico_sce_list = loadCbSceList(calico_libs,
    base_path="/mnt/accessory/seq_data/calico",
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

gtex_sce_list = loadCbSceList(gtex_libs[!gtex_libs %in% skipped],
    base_path="/mnt/accessory/seq_data/gtex",
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

#For GTEX this is the same as their chipwell barcode / donor_id
gtex_all$participant_id = gtex_all$donor_id

# for calico, we need to map this information ourselves
# create a vector of <library>__<chipwell> and use that to map to participant_id

# load metadata
calico_chip_well_to_donor_df = read.csv("~/sc-online/notebook_data/pd/calico_kf_master_20240307.csv", sep = ",", header = TRUE)

# drop some rows that shouldn't exist
calico_chip_well_to_donor_df = calico_chip_well_to_donor_df[
    !(calico_chip_well_to_donor_df$Correct.Donor.ID == "PD0905" & calico_chip_well_to_donor_df$Chipwell.barcode == "207762960003_R02C01"),]

calico_chip_well_to_donor_df$participant_id = calico_chip_well_to_donor_df$Correct.Donor.ID
calico_chip_well_to_donor_map = setNames(
    calico_chip_well_to_donor_df$participant_id,
    paste0(calico_chip_well_to_donor_df$Chipwell.barcode, "_1")
)

# TODO fix the heck out of this
# this is super janky, but like, use the terra manifest when it's available and not in the KF manifest
calico_chip_well_to_donor_map[["207762950055_R08C02_1"]] = "2018-130"
calico_chip_well_to_donor_map[["207762960003_R01C02_1"]] = "2017-037"
calico_chip_well_to_donor_map[["206954930093_R10C01_1"]] = "PD0809"
calico_chip_well_to_donor_map[["207758480019_R01C01_1"]] = "2019-114"
calico_chip_well_to_donor_map[["206954930093_R08C02_1"]] = "PD0899"
calico_chip_well_to_donor_map[["207762950086_R12C02_1"]] = "2019-017"
calico_chip_well_to_donor_map[["207762950086_R08C02_1"]] = "2019-071"
calico_chip_well_to_donor_map[["207762950086_R02C01_1"]] = "PD0878"
calico_chip_well_to_donor_map[["207762950086_R04C01_1"]] = "PD0876"
calico_chip_well_to_donor_map[["207762950108_R10C02_1"]] = "PD0985"
calico_chip_well_to_donor_map[["207762950055_R07C02_1"]] = "PDC094"
calico_chip_well_to_donor_map[["207762950108_R07C02_1"]] = "PD0905"
calico_chip_well_to_donor_map[["206954930060_R07C01_1"]] = "PD0903"
calico_chip_well_to_donor_map[["207762960003_R12C02_1"]] = "PD0784"
calico_chip_well_to_donor_map[["206954930011_R10C01_1"]] = "PD0825"

# create the mapping in the colData
calico_all$participant_id = calico_chip_well_to_donor_map[calico_all$donor_id]

# TODO fix after we resolve the "206954930010_R05C02_1" issue
# there is one chipwell barcode that seems to map to different participants depending on the library??
calico_all$participant_id[calico_all$library == "pCalicoPDsHSrSND9id230921D9" & calico_all$donor_id == "206954930010_R05C02_1"] = "PDC091"
calico_all$participant_id[(
    calico_all$library %in% c("pCalicoPDsHSrSNSN_VTAiPoold230719K1", "pCalicoPDsHSrSNSN_VTAiPoold230719K2") 
    & calico_all$donor_id == "206954930010_R05C02_1")] = "PDC110"

sces_all = .mycBindFn(list(calico_all, gtex_all))

#we would like to append metadata about library_prep to each library 
sces_all$early = sces_all$library %in% EARLY_LIBS

print("Splitting by sort and participant")

# Split by sort and then by donor. Save.
sce_nurr = sces_all[, colData(sces_all)$sort == 'nurr']
sce_dapi = sces_all[, colData(sces_all)$sort == 'dapi']

cb_sce_nurr_donor_list = .mySplitObject(sce_nurr, "participant_id")
cb_sce_dapi_donor_list = .mySplitObject(sce_dapi, "participant_id")

print("Saving SCE lists")

qsave(cb_sce_nurr_donor_list, file.path(BASE_PATH, "pd_all/cb_sce_donor_list_nurr_no_subset_20240308.qs"))
qsave(cb_sce_dapi_donor_list, file.path(BASE_PATH, "pd_all/cb_sce_donor_list_dapi_no_subset_20240308.qs"))

# Harmonize. Takes 90 mins for each!
print("Generating harmonized NURR Seurat")
seurat_nurr_merged = rawSceToHarmonizedSeurat(
    cb_sce_nurr_donor_list,
    split_col = "participant_id",
    harmony_group_by_vars = "participant_id",
    n_donors_hvg = NULL, # have to be in half of the participants
    n_var_features = 5000,
    n_dims_use = 50,
    res = 1.0)
print("Saving harmonized NURR Seurat")
qsave(seurat_nurr_merged, file.path(BASE_PATH, "pd_all/seurat_nurr_merged_no_subset_20240308"))

print("Generating harmonized DAPI Seurat")
seurat_dapi_merged = rawSceToHarmonizedSeurat(
    cb_sce_dapi_donor_list,
    split_col = "participant_id",
    harmony_group_by_vars = "participant_id",
    n_donors_hvg = NULL, # have to be in half of the participants
    n_var_features = 5000,
    n_dims_use = 50,
    res = 1.0)
print("Saving harmonized DAPI Seurat")
qsave(seurat_dapi_merged, file.path(BASE_PATH, "pd_all/seurat_dapi_merged_no_subset_20240308"))
