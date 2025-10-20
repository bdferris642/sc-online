
# Detect script path when running via Rscript
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir <- dirname(normalizePath(script_path))
  setwd(script_dir)
  message("Working directory set to: ", script_dir)
} else {
  stop("Cannot determine script path. Are you running via Rscript?")
}

# path to a tsv file with columns "score_name" and "genes". 
# The genes column should be in the format of a comma-delimited string.

print("**************** LOADING LIBRARIES ****************")
suppressWarnings(suppressMessages({
    library(dplyr)
    library(getopt)
    library(qs)
    library(lme4)
    library(Seurat)
    library(UCell)
    source("../../utils.R")
    library(glue)
}))

g = glue::glue

spec <- matrix(c(
    'sobj-path', 'p', 1, "character",
    'contrast-col', 'cc', 1, "character",
    'covs', 'c', 2, "character",
    'subject-id-col', 'sc', 1, "character",
    'column-base', 'b', 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)

SOBJ_PATH = opt[["sobj-path"]]
if (is.null(opt[['contrast-col']])) {
    CONTRAST_COL = "case_control"   
} else {
    CONTRAST_COL = opt[["contrast-col"]]
}

if (is.null(opt[['subject-id-col']])) {
    SUBJECT_ID_COL = "donor_id"   
} else {
    SUBJECT_ID_COL = opt[["subject-id-col"]]
}

if (is.null(opt[['column-base']])) {
    COLUMN_BASE = "PC_"   
} else {
    COLUMN_BASE = opt[["column-base"]]
}

if(! is.null(opt[['covs']])) {
    COVS = strsplit(opt[['covs']], ",")[[1]]
    COV_STR = g('{paste(COVS, collapse=" + ")} ')
} else {
    COVS = c()
    COV_STR = ""
}

print( paste("Reading Seurat object from", SOBJ_PATH))
sobj = load_obj(SOBJ_PATH)

pvals = list()
input_df = sobj@meta.data
factor_cols = colnames(input_df)[ grepl("cnmf_", colnames(input_df))]
for (cname in factor_cols){
    print(cname)
    formula = as.formula(g("{cname} ~ {CONTRAST_COL} + {COV_STR} + (1|{SUBJECT_ID_COL})"))
    null_formula = as.formula(g("{cname} ~ {COV_STR} + (1|{SUBJECT_ID_COL})"))

    print("Full formula:")
    print(formula)

    print("Null formula:")
    print(null_formula)

        
    # train a linear model on the formula using the input data
    model = lmer(formula, data = input_df, REML = TRUE)
    null_model = lmer(null_formula, data = input_df, REML = TRUE)

    # calculate the p-value for the model
    anova_result = anova(model, null_model, ddf = "Satterthwait")
    p_value = anova_result$`Pr(>Chisq)`[2]
    pvals[[cname]] = p_value
}

print("P-values:")
print(pvals)

sig_factors = factor_cols[pvals < 0.05]

outfile = g("{dirname(SOBJ_PATH)}/{COLUMN_BASE}.txt")
if (! file.exists(outfile)){
    file.create(outfile)
}

if (length(sig_factors) > 0){
    sig_factor_str = paste(sig_factors, collapse=",")
    cov_str_csep = paste(COVS, collapse=",")
    command = paste(
        g("~/sc-online/scripts/run-nebula-parallel.sh --path {SOBJ_PATH}"),
        g("--id-col {SUBJECT_ID_COL} --offset-col nCount_RNA --covs {cov_str_csep},{sig_factor_str}"),
        g("--n-folds 60 --n-cores 60")
    )
    write(command, file=outfile, append=TRUE)
}
