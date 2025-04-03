# TODO: the id_cols / rownames, names_from / colnames, and values_from / values_to parameters are hard coded here. Parametrize

# Runs the mashr pipeline to estimate the relative sharing / uniqueness of effect sizes

# Input: a path to a long format eQTL file saved as RDS
# These eQTLs should be generated from the OSCA pipeline, and have a column uniquely identifying SNP-probe pairs
# every SNP-probe pair must be present in all groups on which you would like to perform the analysis

suppressMessages(suppressWarnings(library(ashr)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(mashr)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(tidyr)))

set.seed(642)

spec <- matrix(c(
    'path', 'p', 1, "character",
    'num-random', 'n', 1, "numeric",
    'padj-thresh', 't', 1, "numeric",
    'eps', 'e', 1, "numeric"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
PATH = opt[["path"]]

if (is.null(opt[["num-random"]])) {
    NUM_RANDOM = 1000000
} else {
    NUM_RANDOM = opt[["num-random"]]
}

if (is.null(opt[["padj-thresh"]])) {
    PADJ_THRESH = 0.01
} else {
    PADJ_THRESH = opt[["padj-thresh"]]
}

if (is.null(opt[["eps"]])) {
    EPS = 1e-6
} else {
    EPS = opt[["eps"]]
}

d = dirname(PATH)
slogan = gsub(".rds", "", basename(PATH))


cat("RUNNING MASHR PIPELINE\n")
cat("WITH PARAMETERS:\n")
cat("PATH: ", PATH, "\n")
cat("NUM_RANDOM: ", NUM_RANDOM, "\n")
cat("PADJ_THRESH: ", PADJ_THRESH, "\n")
cat("EPS: ", EPS, "\n")

cat("\nREADING DATA\n")
eqtl_long = readRDS(PATH)

# ANY SE VALUE LESS THAN EPS IS SET TO EPS
eqtl_long$SE = ifelse(eqtl_long$SE < EPS, EPS, eqtl_long$SE)

# POSITIVE BETAS < EPS --> EPS
eqtl_long$b <- ifelse(eqtl_long$b >= 0 & eqtl_long$b < EPS, EPS, eqtl_long$b)

# NEGATIVE BETAS > -EPS --> -EPS
eqtl_long$b <- ifelse(eqtl_long$b < 0 & eqtl_long$b > -EPS, -EPS, eqtl_long$b)


cat("PIVOTING DATA WIDE\n")
# Expect data in a LONG format (i.e. one row per SNP per cell class). Pivot Wide
eqtl_wide_b = (
    eqtl_long[,c("SNP_probe", "cell_class", "b")]
    %>% pivot_wider(id_cols = "SNP_probe", names_from = "cell_class", values_from = "b")
    %>% arrange(SNP_probe)
    %>% as.data.frame()
)

eqtl_wide_se = (
    eqtl_long[,c("SNP_probe", "cell_class", "SE")]
    %>% pivot_wider(id_cols = "SNP_probe", names_from = "cell_class", values_from = "SE")
    %>% arrange(SNP_probe)
    %>% as.data.frame()
)

eqtl_wide_padj = (
    eqtl_long[,c("SNP_probe", "cell_class", "padj_snp")]
    %>% pivot_wider(id_cols = "SNP_probe", names_from = "cell_class", values_from = "padj_snp")
    %>% arrange(SNP_probe)
    %>% as.data.frame()
)


rownames(eqtl_wide_b) = eqtl_wide_b$SNP_probe
rownames(eqtl_wide_se) = eqtl_wide_se$SNP_probe
rownames(eqtl_wide_padj) = eqtl_wide_padj$SNP_probe

eqtl_wide_b$SNP_probe = NULL
eqtl_wide_se$SNP_probe = NULL
eqtl_wide_padj$SNP_probe = NULL

cat("NUM SNP-PROBE PAIRS: ", nrow(eqtl_wide_b), "\n")

eqtl_wide_b = as.matrix(eqtl_wide_b[complete.cases(eqtl_wide_b),])
eqtl_wide_se = as.matrix(eqtl_wide_se[complete.cases(eqtl_wide_se),])
eqtl_wide_padj = as.matrix(eqtl_wide_padj[complete.cases(eqtl_wide_padj),])

cat("NUM SNP-PROBE PAIRS COMMON ACROSS ALL CELL CLASSES: ", nrow(eqtl_wide_b), "\n")

# identify a subset of strong tests. These are ones with minimum padj < PADJ_THRESH
# take min across columns for each row
padj_mins = apply(eqtl_wide_padj, 1, min)
strong_subset = rownames(eqtl_wide_padj)[padj_mins < PADJ_THRESH]

num_strong = length(strong_subset)
num_orig = nrow(eqtl_wide_b)
cat("Strong tests: ", num_strong, "\n")
cat("Total tests: ", num_orig, "\n")
cat("Fraction of Strong Tests: ", round(num_strong / num_orig, 5), "\n")

cat("ESTIMATING NULL CORRELATION STRUCTURE WITH ", NUM_RANDOM, " RANDOM SAMPLES\n")

# identify a random subset of tests. 
random_subset = sample(1:nrow(eqtl_wide_b), NUM_RANDOM)

# Estimate Correlation Structure in the null tests from random data 
data.temp = mash_set_data(
    eqtl_wide_b[random_subset,],
    eqtl_wide_se[random_subset,]
)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

cat("SUBSETTING DATA INTO 'RANDOM' & 'STRONG' SUBSETS\n")
# Now we can set up our main data objects with this correlation structure in place:
data.random = mash_set_data(
    eqtl_wide_b[random_subset,],
    eqtl_wide_se[random_subset,],
    V=Vhat
)
data.strong = mash_set_data(
    eqtl_wide_b[strong_subset,],
    eqtl_wide_se[strong_subset,], 
    V=Vhat
)

cat("RUNNING PCA\n")
# Use strong tests to set up data-driven covariances
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
U.c = cov_canonical(data.random)

cat("SETTING UP CUSTOM COVARIANCES\n")
# Add to the covariances neuron-specific and glia-specific covariances
# cell classes = c("astro", "da", "mg", "nonda", "oligo",)
astro_col_ind = which(colnames(eqtl_wide_b) == "astro")
da_col_ind = which(colnames(eqtl_wide_b) == "da")
mg_col_ind = which(colnames(eqtl_wide_b) == "mg")
nonda_col_ind = which(colnames(eqtl_wide_b) == "nonda")
oligo_col_ind = which(colnames(eqtl_wide_b) == "oligo")

neurons_only = matrix(0, nrow=5, ncol=5)
neurons_only[da_col_ind, da_col_ind] = 1
neurons_only[da_col_ind, nonda_col_ind] = 1
neurons_only[nonda_col_ind, da_col_ind] = 1
neurons_only[nonda_col_ind, nonda_col_ind] = 1

glia_only = matrix(0, nrow=5, ncol=5)
glia_only[astro_col_ind, astro_col_ind] = 1
glia_only[astro_col_ind, mg_col_ind] = 1
glia_only[astro_col_ind, oligo_col_ind] = 1
glia_only[mg_col_ind, astro_col_ind] = 1
glia_only[mg_col_ind, mg_col_ind] = 1
glia_only[mg_col_ind, oligo_col_ind] = 1
glia_only[oligo_col_ind, astro_col_ind] = 1
glia_only[oligo_col_ind, mg_col_ind] = 1
glia_only[oligo_col_ind, oligo_col_ind] = 1

U.custom = list(
    neurons_only = neurons_only,
    glia_only = glia_only
)

cat("FITTING MASH -- RANDOM DATA\n")
m = mash(data.random, Ulist = c(U.ed,U.c,U.custom))

cat("FITTING MASH -- STRONG DATA\n")
# Finally, compute posterior summaries etc for any subset of tests using the above mash fit. 
# Here we do this for the strong tests. We do this using the same mash function as above, 
# but we specify to use the fit from the previous run of mash by specifying
# `g=get_fitted_g(m), fixg=TRUE`. (In mash the parameter g is used to denote the mixture model which we learned above.)
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)


cat("SAVING DATA\n")
saveRDS(m2, file.path(d, paste0(slogan, "__mash_results_sig.rds")))
saveRDS(m, file.path(d, paste0(slogan, "__mash_results_random.rds")))
saveRDS(data.strong, file.path(d, paste0(slogan, "__mash_data_sig.rds")))
saveRDS(data.random, file.path(d, paste0(slogan, "__mash_data_random.rds")))
