################ LIBRARIES ################
library(dplyr)
library(Matrix)
library(qs)
library(Seurat)
library(SingleCellExperiment)

################ INPUTS ################

THIS_BASE_PATH = "/broad/macosko/data/libraries" # path to your libraries
MOLECULE_INFO_BASENAME = "outs/molecule_info.h5" # path to counts matrix within each library
FILTERED_COUNTS_BASENAME = "outs/filtered_feature_bc_matrix.h5" # path to molecule info matrix within each library
OUTNAME = "/broad/macosko/mehrdad/seurat_objects/merged_seurat.qs"  # output path of merged seurat object
SUPPLEMENTARY_METADATA_PATH = "" # path to a csv with metadata about each LIBRARY. Must contain column called "library"
PARTICIPANT_METADATA_PATH = "" # path to a csv with metadata about each PARTICIPANT. Must contain column called "chip_well_barcode"
VIREO_BASE_PATH = "/broad/macosko/data/multiome_demultiplex/vireo_output" # path to the vireo output for each library
VIREO_DONOR_IDS_BASENAME = "donor_ids.tsv" # path to the vireo donor assignment file within each library


################ FUNCTIONS ##################
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
        seurat_merged = Reduce(function(x, y) merge(x, y, project=project, add.cell.ids=cell_ids), 
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
    # a dataframe, whose rownames are the colnames of umi_dgc
    
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

    # MEHRDAD: change these depending on what is in your participant_metadata
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

    # MEHRDAD: change these depending on what is in your participant_metadata
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


################ MAIN #######################

# this is pseudo-code to read dcg matrices into seurat objects, assuming you have a list of library names called libs
# and that the requistite data is located at the paths specified in the code below
# will build a seurat object for every library

# the locations of these files within every library folder on your vm
libs = c() # add your list of library names here

expt_df = read.table(SUPPLEMENTARY_METADATA_PATH, sep=",", header=TRUE)
participant_metadata_df = read.table(PARTICIPANT_METADATA_PATH, sep=",", header=TRUE)

filtered_dgc_list = list()
mol_df_list = list()
sobj_list = list()

for (library in names(libs)){

    expt_df_subset = expt_df[expt_df$library == library,]

    vireo_donor_ids_path = file.path(VIREO_BASE_PATH, library, VIREO_DONOR_IDS_BASENAME)
    vireo_donor_ids = read.table(vireo_donor_ids_path, sep="\t", header=TRUE)

    filtered_counts_path = file.path(THIS_BASE_PATH, library, FILTERED_COUNTS_BASENAME)
    molecule_info_path = file.path(THIS_BASE_PATH, library, MOLECULE_INFO_BASENAME)

    print('READING FILTERED DGC .H5')
    # always read this
    filtered_dgc = Read10X_h5(filtered_counts_path)
    filtered_dgc_list[[library]] = filtered_dgc

    print("FETCHING MOLECULE INFO FROM .H5")
    
    
    # need the molecular info to get pct_introni
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

    # only consider cells that are in the filtered umi data
    mol_df = mol_df[mol_df$barcode %in% colnames(filtered_dgc),]

        
    meta_df = createMetaDataFromDGCMatrix(
        expt_df=expt_df_subset,
        umi_dgc=filtered_dgc, 
        mol_dgc=mol_df,
        donor_ids=vireo_donor_ids,
        participant_metadata=participant_metadata_df
    )

    sobj = CreateSeuratObject(counts = filtered_dgc, project = library, meta.data = md)
    sobj = RenameCells(sobj, new.names = paste(library, colnames(sobj), sep="_"))
    mt_genes = rownames(sobj)[grepl("^MT-", rownames(sobj))]
    sobj$pct_mito = 100 * colSums(sobj@assays$RNA@counts[mt_genes, ]) / colSums(sobj@assays$RNA@counts)
    sobj_list[[library]] = sobj
}

merged_sobj = mergeSeuratListWithMetadata(sobj_list)
qsave(merged_sobj, OUTNAME)
