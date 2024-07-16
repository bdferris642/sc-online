### This script is used to convert 10X Genomics h5 files to SingleCellExperiment objects
### It also adds metadata to the SingleCellExperiment object
### Including the number of UMIs, number of genes, and number of reads per cell, percent intronic, and percent mitochondrial
### It also adds metadata from the Vireo donor assignment file (prob_donor, prob_doublet, etc.)
### And, using this donor_id, joins to participant-derived metadata (age, sex, case control status)

## INPUTS ##
# lib-path: str, a full path to a newline delimited file of long (prefix-containing) library names
# base-path: str, a full path to the base directory, ultimately containing 10X Genomics data
# this dir should be like <base-path>/<data-source>/<library>/<h5-files>
# data-ing-dirname: str, the name of the output directory to save the SingleCellExperiment objects within lib-sepcific folders
# vireo-donor-ids-basename: str, the name of the file containing the Vireo donor assignments within lib-specific folders
# make-filtered: bool, whether to make and save the filtered SingleCellExperiment objects
# make-raw: bool, whether to make and save the raw SingleCellExperiment objects along with the associated col data
# make-mol-info: bool, whether to make and save the molecule info dataframes
# make-cellbender: bool, whether to make and save the CellBender SingleCellExperiment objects

## OUTPUTS ##
# For each library, the following SingleCellExperiments objects are saved  (depending on the input booleans) 
# to to <lib>/<data-ing-dirname> as .rds files
# molecule_info.rds: df derived from the molecule info.h5, with barcode, umi type (intronic or exonic), and read count
# data_sce.rds: the SCE of filtered 10X Genomics output
# raw_feature_bc_matrix_col_data.rds: the df of colData of the raw 10X Genomics output
# raw_feature_bc_matrix.rds: SCE of raw 10X Genomics output
# cb_data_sce_FPR_<FPR>.rds: the SCE of filtered CellBender output


###################### LIBRARIES ######################

source("~/sc-online/extraGenes.R")
source("~/sc-online/getData.R")
source("~/sc-online/myExpSetCreatorFn.R")
source("~/sc-online/utils.R")

library(caret)
library("DropletUtils")
library(dplyr)
library(getopt)
library(Matrix)
library(rhdf5)
library(rlang)
library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(xml2)


###################### PARSE ARGUMENTS ######################

spec <- matrix(c(
    'lib-path', 'l', 1, "character",
    'base-path', 'b', 1, "character",
    'data-ing-dirname', 'd', 1, "character",
    'vireo-donor-ids-basename', 'v', 1, "character",
    'make-filtered', 'mf', 1, "logical",
    'make-raw', 'mr', 1, "logical",
    'make-mol-info', 'mm', 1, "logical",
    'make-cellbender', 'mcb', 1, "logical",
    'supplementary-metadata-path', 's', 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

LIB_PATH = ifelse(
    is.null(opt[['lib-path']]), 
    "~/missing-gtex.txt", 
    opt[['lib-path']])
BASE_PATH = ifelse(
    is.null(opt[['base-path']]), 
    "/mnt/accessory/seq_data", 
    opt[['base-path']])
DATA_ING_DIRNAME = ifelse(
    is.null(opt[['data-ing-dirname']]), 
    "ingested_data_no_subset", 
    opt[['data-ing-dirname']])
VIREO_DONOR_IDS_BASENAME = ifelse(
    is.null(opt[['vireo-donor-ids-basename']]), 
    "vireo_outs/no_subset/donor_ids.tsv", 
    opt[['vireo-donor-ids-basename']])
MAKE_AND_SAVE_FILTERED_SCE_RDS = ifelse(
    is.null(opt[['make-filtered']]), 
    TRUE, 
    opt[['make-filtered']])
SAVE_RAW_RDS = ifelse(
    is.null(opt[['make-raw']]), 
    FALSE, 
    opt[['make-raw']])
SAVE_MOLECULE_INFO_RDS = ifelse(
    is.null(opt[['make-mol-info']]), 
    FALSE, 
    opt[['make-mol-info']])
SAVE_CELLBENDER_RDS = ifelse(
    is.null(opt[['make-cellbender']]), 
    TRUE, 
    opt[['make-cellbender']])
SUPPLEMENTARY_METADATA_PATH = ifelse(
    is.null(opt[['supplementary-metadata-path']]), 
    "~/sc-online/notebook_data/pd_lib_info_20240510.csv", 
    opt[['supplementary-metadata-path']]
)
    
###################### CONSTANTS ######################



ORGANISM = "human"
MIN_EXP_CELLS=0

MOLECULE_INFO_BASENAME = "cr_outs/molecule_info.h5"
RAW_COUNTS_BASENAME = "cr_outs/raw_feature_bc_matrix.h5"
FILTERED_COUNTS_BASENAME = "cr_outs/filtered_feature_bc_matrix.h5"

PARTICIPANT_METADATA_PATH = "~/sc-online/notebook_data/chip_well_barcode-to-participant_id__fp_coalesced.csv"
FPR = 0.01

CALICO_DIRS_LONG = readLines("~/calico-libs-long.txt")
GTEX_DIRS_LONG = readLines("~/gtex-libs-long.txt")

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

# MAJOR TODO: FIX WHEN THIS SAMPLE SWAP HAS BEEN DISAMBIGUATED
AMBIGUOUS_BARCODES = c("206954930011_R11C01_1", "206954930010_R11C01_1")

###################### HELPER FUNCTIONS ######################



# Extract text after the first underscore in each directory name if it begins with a date
# Function to extract text based on the pattern
extract_text <- function(dir) {
    # Regular expression to match one or more digits followed by an underscore
    pattern <- "^[0-9]+_"
  
    # If the directory name matches the pattern, extract the part after the underscore
    # Otherwise, return the whole directory name
    if (grepl(pattern, dir)) {
        return(sub(pattern, "", dir))
    } else {
        return(dir)
    }
}


createMetaDataFromDGCMatrix=function(
    expt_df,
    umi_dgc, 
    mol_dgc,
    donor_ids,
    participant_metadata){
    
    # INPUTS
    # expt_df, a row of a dataframe 
    # umi_dgc, a dgc matrix derived from a Read 10X filtered feature matrix h5 file
    # mol_dgc, a dgc matrix derived from a Read 10X molecule info h5 file
    # donor_ids, a dataframe derived from vireo's donor assignment tsv file
    # participant_metadata, a dataframe that maps donor_ids (chip_well_barcodes) to participant_ids along with associated metadata

    # OUTPUTS 
    #a dataframe, whose rownames are the colnames of umi_dgc
    
    # output containing columns:
    #   nUMI: total number of UMIs within a cell
    #   nGene: number of unique genes expressed within a cell
    #   nRead: number of total reads within a cell
    #   donor_id: str, either the ID of the donor as specified in the donor_ids.tsv file, or 'unassigned'
    #   prob_donor: float, the probability that the cell is from the donor specified in donor_id
    #   prob_doublet: probability that the droplet represents a doublet, as determined by vireo
    #   prob_singlet: 1-prob_doublet
    # and ALL COLUMNS IN `expt_df` (this should always include 'library')

    umi_dgc_rownames = rownames(umi_dgc)
    umi_dgc_colnames = colnames(umi_dgc)

    # only consider cells that are in the filtered umi data
    mol_df = mol_dgc[mol_dgc$barcode %in% colnames(umi_dgc),]
    mol_df_grouped = mol_df %>% group_by(barcode) %>% summarize(nUmi=n(), nRead=sum(count), pct_intronic=sum(umi_type==0)/nUmi)
    rownames(mol_df_grouped) = mol_df_grouped$barcode

    # match order of other matrices / dfs to barcodes of umi_dgc 
    mol_df_grouped_matched = mol_df_grouped[umi_dgc_colnames, ]

    # join donor_ids to manifest
    donor_ids = merge(
        donor_ids,
        participant_metadata,
        all.x=TRUE, # 'unassigned' and 'doublet' don't get any metadata, but we do need the rows for now
        by.x = "donor_id",
        by.y = "chip_well_barcode")        

    rownames(donor_ids) = donor_ids$cell
    donor_ids_matched = donor_ids[umi_dgc_colnames,]

    # sum within cells
    nGene = colSums(umi_dgc > 0)

    meta_df = data.frame(
        mol_df_grouped_matched$nUmi,
        nGene, 
        mol_df_grouped_matched$nRead,
        mol_df_grouped_matched$pct_intronic,
        donor_ids_matched$donor_id, # note that donor_ids_matched$donor_id is sometimes 'unassigned', even with high prob_max...
        donor_ids_matched$prob_max,
        1-donor_ids_matched$prob_doublet,
        donor_ids_matched$prob_doublet,
        donor_ids_matched$coalesced_participant_id,
        donor_ids_matched$age,
        donor_ids_matched$sex,
        donor_ids_matched$case_control
    )

    colnames(meta_df) = c(
        "nUMI", "nGene", 'nRead', 'pct_intronic',
        'donor_id', 'prob_donor', 'prob_singlet', 'prob_doublet',
        "participant_id", "age", "sex", "case_control")
    for (cname in colnames(expt_df)){
        meta_df[[cname]] = expt_df[[cname]]
    }
    rownames(meta_df) = colnames(umi_dgc)
    return(meta_df)
}     

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


###################### MAIN ######################

DIRS_TO_READ_LONG = readLines(LIB_PATH)
DIRS_TO_READ = sapply(DIRS_TO_READ_LONG, extract_text)
names(DIRS_TO_READ_LONG) = DIRS_TO_READ


# TODO: remove this "skipped" block once all gtex libraries have been processed
SKIPPED  = c()
for (d in DIRS_TO_READ){
    if (d == ""){next}
    if (!file.exists(file.path(BASE_PATH, "gtex", d, "cr_outs", "_cmdline"))){next}
    if (!file.exists(file.path(BASE_PATH, "gtex", d, "vireo_outs", "no_subset", "donor_ids.tsv"))){
        SKIPPED  = c(SKIPPED , d)
    }
}
# these libraries are bad and should probably jsut be removed from the gtex list
SKIPPED  = c(SKIPPED 
    , "pCalico_GTExsHSrSNA11iNURRd231120"
    ,"pCalico_GTExsHSrSNB11iNURRd231120"
    ,"pCalico_GTExsHSrSNC11iDAPId231120"
    ,"pCalico_GTExsHSrSND11iNURRd231120"
    ,"pCalico_GTExsHSrSNE11iNURRd231120"
    ,"pCalico_GTExsHSrSNF11iDAPId231120"
    ,""
)

filtered_dgc_list = list()
sce_list = list()
raw_dgc_list = list()
raw_cd_list = list()
mol_info_list = list()
vireo_donor_ids_list = list()

# read in participant metadata
participant_metadata_df = read.table(PARTICIPANT_METADATA_PATH, sep=",", header=TRUE)

# read in lib-level info
expt_df = read.table(SUPPLEMENTARY_METADATA_PATH, sep=",", header=TRUE)
expt_df = expt_df[,
    c("library", "library_long", "source", "sort", "region", "flowcell")
]

# note if library was sorted using early sorting protocol
expt_df$early = FALSE 
expt_df$early[expt_df$library %in% EARLY_LIBS] = TRUE

for (long_name in names(DIRS_TO_READ)){
    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    if(long_name == ''){
        next
    }

    name = DIRS_TO_READ[[long_name]]
    if (long_name %in% CALICO_DIRS_LONG){
        this_base_path = file.path(BASE_PATH, 'calico')
    } else if (long_name %in% GTEX_DIRS_LONG){
        this_base_path = file.path(BASE_PATH, 'gtex')
    }
    
    if (name %in% SKIPPED){
        next
    }
    
    print(name)

    raw_counts_path = file.path(this_base_path, name, RAW_COUNTS_BASENAME)
    filtered_counts_path = file.path(this_base_path, name, FILTERED_COUNTS_BASENAME)
    molecule_info_path = file.path(this_base_path, name, MOLECULE_INFO_BASENAME)
    vireo_donor_ids_path = file.path(this_base_path, name, VIREO_DONOR_IDS_BASENAME)

    vireo_donor_ids = read.table(vireo_donor_ids_path, sep="\t", header=TRUE)

    # MAJOR TODO: FIX WHEN THIS SAMPLE SWAP HAS BEEN DISAMBIGUATED
    vireo_donor_ids$donor_id[vireo_donor_ids$donor_id %in% AMBIGUOUS_BARCODES] = "unassigned"

    print('CREATING DATA OBJECT DIRECTORY')
    dir.create(file.path(this_base_path, name, DATA_ING_DIRNAME), showWarnings = FALSE)

    print('READING FILTERED DGC .H5')
    # always read this
    filtered_dgc = Read10X_h5(filtered_counts_path)
    filtered_dgc_list[[name]] = filtered_dgc

    if (MAKE_AND_SAVE_FILTERED_SCE_RDS){
    
        print("FETCHING MOLECULE INFO FROM .H5")
        # create a molecular data dataframe my fetching specific fields from the h5 and grouping by barcode
        mol_h5_file = file.path(molecule_info_path)
        h5fetch = function(x){return(h5read(mol_h5_file, x))}
        mol_df = data.frame(
            barcode=h5fetch("barcodes")[h5fetch("barcode_idx")+1] %>% paste0("-1"),
            feature=h5fetch("features/name")[h5fetch("feature_idx")+1],
            umi=h5fetch("umi"),
            umi_type=h5fetch("umi_type"),
            count=h5fetch("count")
        )
        
        if (SAVE_MOLECULE_INFO_RDS){
            print("SAVING MOLECULE INFO .RDS")
            saveRDS(mol_df, file=file.path(this_base_path, name, DATA_ING_DIRNAME, "molecule_info.rds"))
        }

        print("COLLECTING METADATA")
        # only consider cells that are in the filtered umi data
        mol_df = mol_df[mol_df$barcode %in% colnames(filtered_dgc),]

        vireo_donor_ids_list[[name]] = vireo_donor_ids

        expt_df_subset = expt_df[expt_df$library == name,]

        meta_df = createMetaDataFromDGCMatrix(
            expt_df=expt_df_subset,
            umi_dgc=filtered_dgc, 
            mol_dgc=mol_df,
            donor_ids=vireo_donor_ids,
            participant_metadata=participant_metadata_df
        )

        print("BUILDING FILTERED SCE")
        data_sce=.myExpSetCreatorFn(inputExpData=filtered_dgc,
                        organism=ORGANISM,
                        minExpCells=MIN_EXP_CELLS,
                        inputPdata=meta_df,
                        inputFdata=NULL,
                        addExtraAnno=T,
                        ncores=14)
        
        print("SAVING FILTERED SCE .RDS")
        sce_list[[name]] = data_sce
        saveRDS(data_sce, file.path(this_base_path, name, DATA_ING_DIRNAME, "data_sce.rds"))
    }

    if (SAVE_RAW_RDS){
        print("READING RAW DGC .H5")
        raw_counts = Read10X_h5(raw_counts_path)
        raw_dgc_list[[name]] = raw_counts

        cd = data.frame(
            row.names=colnames(raw_counts),
            nUMI=colSums(raw_counts),
            log10_nUMI = log10(colSums(raw_counts) + 1),
            nGene=colSums(raw_counts > 0),
            is_in_filtered=colnames(raw_counts) %in% colnames(filtered_dgc_list[[name]])
        )
        raw_cd_list[[name]] = cd

        print("SAVING RAW DGC COL DATA RDS")
        saveRDS(cd, file.path(this_base_path, name, DATA_ING_DIRNAME, 'raw_feature_bc_matrix_col_data.rds'))
        print("SAVING RAW DGC RDS")
        saveRDS(raw_counts, file.path(this_base_path, name, DATA_ING_DIRNAME, 'raw_feature_bc_matrix.rds'))
    }
    print(paste0("DONE WITH LIBRARY", name))
    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
}

# CELLBENDER RDS MAKER
print("MAKING CELLBENDER RDS's")
for (name in DIRS_TO_READ){
    print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    if (name %in% SKIPPED | name ==""){
        next
    }

    long_name = DIRS_TO_READ_LONG[[name]]

    if (long_name %in% CALICO_DIRS_LONG){
        this_base_path = file.path(BASE_PATH, 'calico')
    } else if (long_name %in% GTEX_DIRS_LONG){
        this_base_path = file.path(BASE_PATH, 'gtex')
    }    
    
    print(name)

    long_name = DIRS_TO_READ_LONG[[name]]

    if (long_name %in% CALICO_DIRS_LONG){
        this_base_path = file.path(BASE_PATH, 'calico')
    } else if (long_name %in% GTEX_DIRS_LONG){
        this_base_path = file.path(BASE_PATH, 'gtex')
    }

    orig_sce = readRDS(file.path(this_base_path, name, DATA_ING_DIRNAME, "data_sce.rds"))
    filtered_cb_path = file.path(this_base_path, name, "cb_outs", paste0(name, "_out_filtered.h5"))

    file_path = filtered_cb_path # the path to the *_out_filtered.h5 file
    h5_data = h5read(file_path, "/matrix")
    counts = h5_data$data
    indices = h5_data$indices
    indptr = h5_data$indptr
    barcodes = h5_data$barcodes
    num_genes = length(h5_data$features$name)  # Number of genes
    num_cells = length(barcodes)  # Number of cells

    counts_matrix = sparseMatrix(i=indices + 1,  # +1 for 1-based indexing in R
                                    p=indptr,
                                    x=counts,
                                    dims=c(num_genes, num_cells))

    rownames(counts_matrix) = h5_data$features$name 
    colnames(counts_matrix) = barcodes
    sce_filt = orig_sce[, colnames(orig_sce) %in% colnames(counts_matrix)]
    counts_matrix_filtered = counts_matrix[, colnames(counts_matrix) %in% colnames(sce_filt)]
    # for whatever reason, there can be duplicated rownames (i.e. genes) in the cellbender output, causing errors later on. Sum these together.
    counts_matrix_filtered = sum_duplicate_rownames_of_dgc_matrix(counts_matrix_filtered) 
    counts_matrix_filtered = counts_matrix_filtered[, match(colnames(sce_filt), colnames(counts_matrix_filtered))]

    cd = colData(sce_filt)
    cd = cd[rownames(cd) %in% colnames(counts_matrix_filtered), ]

    non_qc_cols <- colnames(cd)[!grepl("^QC", colnames(cd))]
    cd = cd[, non_qc_cols]

    h5fetch = function(x){return(h5read(filtered_cb_path, x))}
    cb_cd_supp = data.frame(
            barcodes = h5fetch('/matrix/barcodes'),
            cell_probability =  h5fetch('droplet_latents/cell_probability'),
            background_fraction = h5fetch('droplet_latents/background_fraction')
        )
    rownames(cb_cd_supp) = cb_cd_supp$barcodes
    cb_cd_supp = cb_cd_supp[rownames(cd), ]
    cb_cd_supp = cb_cd_supp[match(rownames(cd), rownames(cb_cd_supp)),]
    cd = cbind(cd, cb_cd_supp)

    cb_sce=.myExpSetCreatorFn(
            inputExpData=counts_matrix_filtered,
            organism=ORGANISM,
            minExpCells=MIN_EXP_CELLS,
            inputPdata=cd,
            inputFdata=NULL,
            addExtraAnno=T,
            ncores=12
        )
    
    # need to re-add important QC data that gets scrambled above for reasons I don't understand
    sce_filt = sce_filt[, match(colnames(cb_sce), colnames(sce_filt))]
    cb_sce$nUMI = sce_filt$nUMI
    cb_sce$nGene = sce_filt$nGene
    cb_sce$nRead = sce_filt$nRead
    cb_sce$pct_intronic = sce_filt$pct_intronic
    
    qc_cols = colnames(colData(sce_filt))[grepl("^QC", colnames(colData(sce_filt)))]
    for (c in qc_cols){
        cb_sce[[c]] = sce_filt[[c]]
    }
    
    print("SAVING CELLBENDER SCE AS .RDS")
    cb_rds_basename = paste0("cb_data_sce_FPR_", FPR, ".rds")
    saveRDS(cb_sce, file.path(this_base_path, name, DATA_ING_DIRNAME, cb_rds_basename))

    print(paste0("DONE WITH LIBRARY", name))
    print(":::::::::::::::::::::::::::::::::::::::::::::")
}

