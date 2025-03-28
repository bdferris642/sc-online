suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))

get_fname_slogan = function(fname) {
    slogan = gsub(".rds", "", basename(fname))
    slogan = gsub("eqtl_", "", slogan)
    return(slogan)
}

spec = matrix(c(
    'base', 'b', 1, "character"
), byrow = TRUE, ncol = 4)

base = getopt(spec)[["base"]]

# list files in base ending in .rds, NOT ending in _sig.rds
files = list.files(base, pattern = "\\.rds$", full.names = T)
files = files[!grepl("_sig.rds", files)]

names(files) = sapply(files, get_fname_slogan)

rds_list = list()
for (slogan in names(files)) {
    cat(paste0("READING ", slogan, "\n"))
    rds = readRDS(files[[slogan]])
    rds$cell_class = slogan
    rds$SNP_probe = paste(rds$SNP, rds$Probe, sep="_")
    rds_list[[slogan]]
}

all_eqtls = do.call(rbind, rds_list)

# cat("READING ASTRO\n")
# astro = readRDS(file.path(base, "eqtl_astro.rds"))
# astro$SNP_probe = paste(astro$SNP, astro$Probe, sep="_")

# cat("READING DA\n")
# da = readRDS(file.path(base, "eqtl_da.rds"))
# da$SNP_probe = paste(da$SNP, da$Probe, sep="_")

# cat("READING MG\n")
# mg = readRDS(file.path(base, "eqtl_mg.rds"))
# mg$SNP_probe = paste(mg$SNP, mg$Probe, sep="_")

# cat("READING OLIGO\n")
# oligo = readRDS(file.path(base, "eqtl_oligo.rds"))
# oligo$SNP_probe = paste(oligo$SNP, oligo$Probe, sep="_")

# cat("READING NONDA\n")
# nonda = readRDS(file.path(base, "eqtl_nonda.rds"))
# nonda$SNP_probe = paste(nonda$SNP, nonda$Probe, sep="_")

# astro$cell_class = "astro"
# da$cell_class = "da"
# mg$cell_class = "mg"
# oligo$cell_class = "oligo"
# nonda$cell_class = "nonda"

# cat("SUBSETTING TO COMMON SNPS\n")
# all_eqtls = rbind(astro, da, mg, oligo, nonda)
# rm(astro)
# rm(da)
# rm(mg)
# rm(oligo)
# rm(nonda)

cat("SUBSETTING TO COMMON SNPS\n")
cell_class_counts = (all_eqtls 
    %>% group_by(SNP_probe) 
    %>% summarise(n_cell_classes=n_distinct(cell_class)) 
    %>% ungroup())
common_snps = (cell_class_counts 
    %>% filter(n_cell_classes == length(files))
    %>% select(SNP_probe))
eqtl_common_snps = all_eqtls %>% filter(SNP_probe %in% common_snps$SNP_probe)

cat("SAVING DATA\n")
saveRDS(eqtl_common_snps, file.path(base, "eqtl_present_in_all.rds"))

cat("SAVING SIG DATA\n")
common_snps_sig_only = eqtl_common_snps %>% filter(padj_snp < 0.01)
snp_probe_common_all_sig_one = unique(common_snps_sig_only$SNP_probe)
common_snps_sig_in_one = eqtl_common_snps %>% filter(SNP_probe %in% snp_probe_common_all_sig_one)
saveRDS(common_snps_sig_in_one, file.path(base, "eqtl_present_in_all_sig_in_one.rds"))

cat("DONE\n")