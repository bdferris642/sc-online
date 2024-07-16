library(qs)
library(Seurat)
library(dplyr)

base_path = "/mnt/accessory/seq_data/pd_all/240514"

swap_participants = function(sobj, pid_pair_list){

    sobj_copy = sobj 

    for (pid_pair in pid_pair_list){
        pid1 = pid_pair[[1]]
        pid2 = pid_pair[[2]]

        sobj_copy$sex[sobj$participant_id == pid1] = sobj$sex[sobj$participant_id == pid2]
        sobj_copy$age[sobj$participant_id == pid1] = sobj$age[sobj$participant_id == pid2]
        sobj_copy$case_control[sobj$participant_id == pid1] = sobj$case_control[sobj$participant_id == pid2]

        sobj_copy$sex[sobj$participant_id == pid2] = sobj$sex[sobj$participant_id == pid1]
        sobj_copy$age[sobj$participant_id == pid2] = sobj$age[sobj$participant_id == pid1]
        sobj_copy$case_control[sobj$participant_id == pid2] = sobj$case_control[sobj$participant_id == pid1]
        
        sobj_copy$participant_id[sobj$participant_id == pid1] = pid2
        sobj_copy$participant_id[sobj$participant_id == pid2] = pid1
    }

    return(sobj_copy)
}

samples = c(
    'dapi_astro_clean.qs',
    'dapi_astro_seurat.qs',
    'dapi_da_clean.qs',
    'dapi_da_clean_hi_numi.qs',
    'dapi_endo_canonical_endo.qs',
    'dapi_endo_clean.qs',
    'dapi_endo_fibroblast.qs',
    'dapi_endo_perictye.qs',
    'dapi_endo_seurat.qs',
    'dapi_immune_clean.qs',
    'dapi_immune_mg_clean.qs',
    'dapi_immune_seurat.qs',
    'dapi_mg_clean.qs',
    'dapi_mg_clean_ctr.qs',
    'dapi_mg_clean_ctr_models.qs',
    'dapi_mg_ctr.qs',
    'dapi_mg_seurat.qs',
    'dapi_nonda_clean.qs',
    "dapi_oligo_seurat.qs",
    "dapi_opc_seurat.qs",
    "dapi_seurat.qs",
    "dapi_seurat_clean.qs",
    "dapi_seurat_clean_ctr.qs",
    "dapi_seurat_ctr.qs",
    "dapi_seurat_da.qs",
    "dapi_seurat_immune.qs",
    "nurr_da_double_clean.qs",
    "nurr_da_double_clean_ctr.qs",
    "nurr_da_double_clean_ctr_sox6_agtr1_nxph1.qs",
    "nurr_da_init.qs",
    "nurr_seurat.qs",
    "nurr_seurat_astro.qs",
    "nurr_seurat_clean.qs",
    "nurr_seurat_da.qs",
    "nurr_seurat_da_clean_clustered.qs",
    "nurr_seurat_da_clean_clustered_ctr.qs",
    'nurr_seurat_endo.qs',
    "nurr_seurat_mg.qs",
    'nurr_seurat_nonda.qs',
    'nurr_seurat_nonda_clean.qs',
    'nurr_seurat_nonda_clean_calico.qs',
    'nurr_seurat_nonda_clean_ctr.qs',
    'nurr_seurat_nonda_clean_tmp.qs',
    'nurr_seurat_oligo.qs',
    'nurr_seurat_opc.qs'
)

file_name = file.path(base_path, "swapped_samples.txt")
for (sample in samples) {
    cat("\n", "Processing", sample, "\n")
    tryCatch({
        sobj = qread(file.path(base_path, sample))
        sobj_swap = swap_participants(sobj, list(c("PDC091", "PDC110")))
        qsave(sobj_swap, file.path(base_path, sample))
        # write the sample to a new line of the file
        write(sample, file = file_name, append = TRUE)
        cat("Successfully processed", sample, "\n")
    }, error = function(e) {
        print(paste("ERROR IN", sample, "!"))
        print(e)
    })
}