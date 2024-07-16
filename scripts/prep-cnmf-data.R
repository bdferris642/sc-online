library(dplyr)
library(Seurat)
library(qs)
library(zellkonverter)
library(SingleCellExperiment)

BASE_PATH = "/mnt/accessory/seq_data/pd_all/240514"
SEURAT_BASENAME = "dapi_mg_clean.qs"
GCLOUD_PATH = "gs://macosko_data/ferris/h5ad_objs/"

if (!dir.exists(file.path(BASE_PATH, "nmf"))) {
    dir.create(file.path(BASE_PATH, "nmf"))
}

slogan = gsub(".qs$", "", SEURAT_BASENAME)
barcode_output_basename = paste0(slogan, ".barcodes.txt")
h5ad_output_basename = paste0(slogan, ".h5ad")

barcode_outpath = file.path(BASE_PATH, "nmf", barcode_output_basename)
h5ad_outpath = file.path(BASE_PATH, "nmf", h5ad_output_basename)

cat("Reading data\n")
sobj = qread(file.path(BASE_PATH, SEURAT_BASENAME))

barcode_df = data.frame(
    CBC=colnames(sobj),
    sample=sobj$participant_id
)

barcode_df$output_col = paste0(barcode_df$CBC, "\t", barcode_df$sample)
output_vector = as.character(barcode_df$output_col)
output_vector = c("CBC\tsample", output_vector)

cat("Writing Barcodes\n")
writeLines(output_vector, barcode_outpath)

cat("Converting Seurat --> SCE\n")
sce <- as.SingleCellExperiment(sobj)

cat("Writing control .h5ad\n")
zellkonverter::writeH5AD(sce, h5ad_outpath)

cat("Copying to GCloud\n")
system(paste("gcloud storage cp", h5ad_outpath, GCLOUD_PATH))
system(paste("gcloud storage cp", barcode_outpath, GCLOUD_PATH))

cat("Done\n")

