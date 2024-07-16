
library(Seurat)
library(Matrix)
library(rhdf5)
library(dplyr)

# THIS CODE TAKES IN
# (1) A LIST OF LIBRARY NAMES
# (2) THE BASE PATH IN WHICH THOSE LIBRARIES LIVE
# (3) THE BASENAMES OF RAW, FILTERED, AND MOLECULE INFO H5's

# THIS SCRIPT THEN LOADS THE OBJECTS INTO MEMORY, STORING THEM AS LISTS FOR FURTHER INVESTIGATION

raw_dgc_list = list()
filtered_dgc_list = list()
mol_df_list = list()

this_base_path = "/PATH/TO/LIBS"
library_list = list(
    "LIBRARY",
    "NAMES"
)

RAW_COUNTS_BASENAME = "raw_feature_bc_matrix.h5"
FILTERED_COUNTS_BASENAME = "filtered_feature_bc_matrix.h5"
MOLECULE_INFO_BASENAME = "molecule_info.h5"

# loop through libraries
for (name in library_list){
    cat(name)
    cat("\n")
    
	raw_counts_path = file.path(this_base_path, name, RAW_COUNTS_BASENAME)
    filtered_counts_path = file.path(this_base_path, name, FILTERED_COUNTS_BASENAME)
    molecule_info_path = file.path(this_base_path, name, MOLECULE_INFO_BASENAME)
	
    # load the filtered counts 
    filtered_dgc = Read10X_h5(filtered_counts_path)
    filtered_dgc_list[[name]] = filtered_dgc

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


    # load the molecule data such that you can map the barcodes to what came out of CellRanger
    # load umi_type to get whether or not the umi was intronic
    h5fetch = function(x){return(h5read(molecule_info_path, x))}
    mol_df = data.frame(
                barcode=h5fetch("barcodes")[h5fetch("barcode_idx")+1] %>% paste0("-1"),
                feature=h5fetch("features/name")[h5fetch("feature_idx")+1],
                umi=h5fetch("umi"),
                umi_type=h5fetch("umi_type"),
                count=h5fetch("count")
        )
        
    # only consider cells that are in the filtered umi data
    mol_df = mol_df[mol_df$barcode %in% colnames(filtered_dgc),]

    mol_df_grouped = mol_df %>% group_by(barcode) %>% 
        summarize(
            nUmi=n(),
            nRead=sum(count), 
            pct_intronic=sum(umi_type==0)/nUmi)
    rownames(mol_df_grouped) = mol_df_grouped$barcode
    mol_df_list[[name]] = mol_df_grouped

}