library(dplyr)
library(ggplot2)
library(Seurat)
library(qs)

library(zellkonverter)
library(SingleCellExperiment)

BASE_PATH = "/mnt/accessory/seq_data/pd_all/240514"


cat("Reading data\n")
nurr_da = qread(file.path(BASE_PATH, "nurr_da_double_clean.qs"))
nurr_da_sox6 = nurr_da[,nurr_da$da_subtype %in% c("SOX6_AGTR1", "SOX6_NXPH1")]
nurr_da_sox6_ctr = nurr_da_sox6[,nurr_da_sox6$case_control == "ctr"]

barcode_df = data.frame(
    CBC=colnames(nurr_da_sox6),
    sample=nurr_da_sox6$participant_id
)

barcode_df_ctr = data.frame(
    CBC=colnames(nurr_da_sox6_ctr),
    sample=nurr_da_sox6_ctr$participant_id
)

barcode_df$output_col = paste0(barcode_df$CBC, "\t", barcode_df$sample)
barcode_df_ctr$output_col = paste0(barcode_df_ctr$CBC, "\t", barcode_df_ctr$sample)

output_vector = as.character(barcode_df$output_col)
output_vector_ctr = as.character(barcode_df_ctr$output_col)

output_vector = c("CBC\tsample", output_vector)
output_vector_ctr = c("CBC\tsample", output_vector_ctr)

cat("Writing Barcodes\n")
writeLines(output_vector, 
    file.path(BASE_PATH, "nmf/nurr_da_double_clean_sox6_agtr1_sox6_nxph1_barcodes.tsv"))
writeLines(output_vector_ctr, 
    file.path(BASE_PATH, "nmf/nurr_da_double_clean_sox6_agtr1_sox6_nxph1_ctr_barcodes.tsv"))

cat("Converting Seurat --> SCE\n")
nurr_da_sox6_ctr_sce <- as.SingleCellExperiment(nurr_da_sox6_ctr)
nurr_da_sox6_sce <- as.SingleCellExperiment(nurr_da_sox6)

cat("Writing control .h5ad\n")
zellkonverter::writeH5AD(nurr_da_sox6_ctr_sce, file.path(BASE_PATH, "nmf/nurr_da_double_clean_sox6_agtr1_sox6_nxph1_ctr.h5ad"))

cat("Writing case + control .h5ad\n")
zellkonverter::writeH5AD(nurr_da_sox6_sce, file.path(BASE_PATH, "nmf/nurr_da_double_clean_sox6_agtr1_sox6_nxph1.h5ad"))