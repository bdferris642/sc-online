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

basenames = sapply(files, basename)
cat('\nREADING FOLLOWING FILES FROM ', base, '\n')
cat(paste0(basenames, "\n"))

rds_list = list()
common_snp_probes = NULL
for (slogan in names(files)) {
    cat(paste0("READING ", slogan, "\n"))
    rds = readRDS(files[[slogan]])
    rds$cell_class = slogan
    rds$SNP_probe = paste(rds$SNP, rds$Probe, sep="_")
    rds_list[[slogan]] = rds

    this_snp_probe_set = unique(rds$SNP_probe)
    if (is.null(common_snp_probes)) {
        common_snp_probes = this_snp_probe_set
    } else {
        common_snp_probes = intersect(common_snp_probes, this_snp_probe_set)
    }
    rm(rds)
    gc(verbose = FALSE)
}

cat("SUBSETTING TO COMMON SNPS\n")
for (slogan in names(files)) {
    cat(paste0("SUBSETTING ", slogan, "\n"))
    rds = rds_list[[slogan]]
    rds = rds %>% filter(SNP_probe %in% common_snp_probes)
    rds_list[[slogan]] = rds
    rm(rds)
    gc(verbose = FALSE)
}

eqtl_common_snps = do.call(rbind, rds_list)
rm(rds_list)


cat("SAVING DATA\n")
saveRDS(eqtl_common_snps, file.path(base, "eqtl_present_in_all.rds"))

cat("SAVING SIG DATA\n")
common_snps_sig_only = eqtl_common_snps %>% filter(padj_snp < 0.01)
snp_probe_common_all_sig_one = unique(common_snps_sig_only$SNP_probe)
common_snps_sig_in_one = eqtl_common_snps %>% filter(SNP_probe %in% snp_probe_common_all_sig_one)
saveRDS(common_snps_sig_in_one, file.path(base, "eqtl_present_in_all_sig_in_one.rds"))

cat("DONE\n")