print("**************** LOADING LIBRARIES ****************")
# Detect script path when running via Rscript
args = commandArgs(trailingOnly = FALSE)
script_path = sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir = dirname(normalizePath(script_path))
  message("Script located in directory: ", script_dir)
} else {
  stop("Cannot determine script path. Are you running via Rscript?")
}

suppressMessages(suppressWarnings({
    library(glue)
    library(dplyr)
    library(getopt)
    library(Matrix)
    library(qs)
    library(Seurat)

    source(file.path(script_dir, "../utils.R"))

    g = glue::glue
}))


spec <- matrix(c(
    'p1', 'p', 1, "character",
    'p2', 'q', 1, "character",
    'out-dir', 'o', 1, "character",
    'out-name', 'n', 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
P1 = opt[['p1']]
P2 = opt[['p2']]
if (is.null(opt[['out-dir']])) {
    BASE = dirname(P1)
    print(g("Warning: --out-dir not specified, using directory of P1: {dirname(P1)}"))
} else {
    BASE = opt[['out-dir']]
}

if (is.null(opt[['out-name']])){
    bname = basename(P2)
    SLOGAN = gsub(".qs", "", bname)
    SLOGAN = gsub("_case_", "_combined_", SLOGAN)
    SLOGAN = gsub("_control_", "_combined_", SLOGAN)
    print(g("Warning: --out-name not specified, using derived name from P2: {SLOGAN}"))
} else {
    SLOGAN = opt[['out-name']]
}

outpath = file.path(BASE, paste0(SLOGAN, ".qs"))


print(g("Reading Seurat object: from {P1}"))
s1 = qread(P1)

print(g("Reading Seurat object: from {P2}"))
s2 = qread(P2)


s1_md_cols = colnames(s1@meta.data)
s2_md_cols = colnames(s2@meta.data)

missing_in_s1 = setdiff(s2_md_cols, s1_md_cols)
missing_in_s2 = setdiff(s1_md_cols, s2_md_cols)

print("Colnames missing from first seurat object:")
print(missing_in_s1)
print("Colnames missing from second seurat object:")
print(missing_in_s2)

for (cname in missing_in_s1){
    s1@meta.data[[cname]] = NA
}

for (cname in missing_in_s2){
    s2@meta.data[[cname]] = NA
}

cnames = union(s1_md_cols, s2_md_cols)
s1@meta.data = s1@meta.data[, cnames]
s2@meta.data = s2@meta.data[, cnames]
combined = combine_sobjs(s1, s2, use_shared_features=T)

print("meta.data after alignment:")
print(head(combined@meta.data))

message("Saving Combined Seurat Object to: ", outpath)
qsave(combined, outpath)
