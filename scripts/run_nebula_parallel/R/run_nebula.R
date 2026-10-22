#!/usr/bin/env Rscript
# Copy of ~/sc-online/scripts/run-nebula.R (NEBULA LN on one gene-chunk Seurat .qs), vendored into
# run_nebula_parallel. Changes vs upstream: (1) self-activation below; (2) the source("../utils.R")
# line is removed (nothing from utils.R is used; it loads biomaRt/sva); (3) nebula::scToNeb() is
# replaced by scToNeb_seurat() (same code, GetAssayData(layer=) for SeuratObject >= 5.4).

# ---- self-activation: re-exec under the workflow env Rscript (env path derived from this file)
local({
    args_all = commandArgs(trailingOnly = FALSE)
    script = normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
    env_prefix = file.path(dirname(dirname(dirname(script))), "micromamba_root", "envs", "run-nebula-parallel-env")
    if (!dir.exists(env_prefix)) stop("workflow env not found: ", env_prefix, " (run run_nebula_parallel/setup.sh)")
    if (normalizePath(R.home()) != normalizePath(file.path(env_prefix, "lib", "R"))) {
        status = system2(file.path(env_prefix, "bin", "Rscript"), c(shQuote(script), shQuote(commandArgs(trailingOnly = TRUE))))
        quit(status = status, save = "no")
    }
})

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
    library(nebula)
    library(future)
    #future::plan(multicore, workers = 4)  

    g = glue::glue
}))

spec = matrix(c(
    'path', 'p', 1, "character",
    'covs', 'c', 2, "character",
    'id-col', 'r', 2, "character",
    'offset-col', 'f', 2, "character",
    'suffix', 's', 2, "character",
    'min-num-cells-per-id', 'm', 2, "integer",
    'pct-nz-threshold', 'z', 2, "double",
    'n-cores', 'n', 2, "integer",
    'out-dir', 'o', 2, "character",
    'clobber', 'b', 2, "logical"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
PATH = opt[['path']]

if (is.null(opt[['covs']])){
    COVS = NULL
} else {
    COVS = strsplit(opt[['covs']], ",")[[1]]
}

if (is.null(opt[['id-col']])){
    ID_COL = "donor_id"
} else {
    ID_COL = opt[['id-col']]
}

if (is.null(opt[['offset-col']])){
    OFFSET_COL = "nCount_RNA"
} else {
    OFFSET_COL = opt[['offset-col']]
}

if (is.null(opt[['pct-nz-threshold']])){
    PCT_NZ_THRESHOLD = 0.01
} else {
    PCT_NZ_THRESHOLD = opt[['pct-nz-threshold']]
}

if (is.null(opt[["min-num-cells-per-id"]])){
    MIN_NUM_CELLS_PER_ID = 5
} else {
    MIN_NUM_CELLS_PER_ID = opt[['min-num-cells-per-id']]
}

if (is.null(opt[['suffix']])){
    SUFFIX = NULL
} else {
    SUFFIX = opt[['suffix']]
}

if (is.null(opt[['n-cores']])){
    N_CORES = 1
} else {
    N_CORES = opt[['n-cores']]
}

if (is.null(opt[['clobber']])){
    CLOBBER = FALSE
} else {
    CLOBBER = opt[['clobber']]
}

geneset_dir = dirname(PATH)
split_dir = dirname(geneset_dir)
geneset = basename(geneset_dir)
slogan = gsub(".qs","__nebula_ln.qs", basename(PATH))

if (! is.null(SUFFIX)){
    slogan = gsub(".qs", g("__{SUFFIX}.qs"), slogan)
}

if (is.null(opt[['out-dir']])){
    de_dir = g("{split_dir}/de_results")
} else {
    de_dir = opt[['out-dir']]
}
if (!dir.exists(de_dir)){
    dir.create(de_dir, recursive=TRUE)
}

outfile = file.path(de_dir, slogan)

if (file.exists(outfile) && ! CLOBBER){
    print(g("**************** OUTPUT FILE: {outfile} already exists; skipping ****************"))
    quit(status=0)
}

print(g ("**************** OUTPUT FILE: {outfile} ****************"))


print(g("**************** LOADING DATA FROM {PATH} ****************"))
sobj = qread(PATH)
sobj$log10_nUMI = log10(sobj$nCount_RNA + 1)

# check if COVS and ID_COL and OFFSET_COL exist in metadata
meta_cols = sort(unique(c(COVS, ID_COL, OFFSET_COL)))
# separately consider interactions containing : or * 
meta_interaction_cols = grepl(":|\\*", meta_cols)
meta_interaction_cols_parts = c()
for (mk in meta_cols[meta_interaction_cols]) {
    parts = unlist(strsplit(mk, ":|\\*"))
    meta_interaction_cols_parts = c(meta_interaction_cols_parts, parts)
}

meta_cols = meta_cols[!meta_interaction_cols]
meta_cols = unique(c(meta_cols, meta_interaction_cols_parts))
missing_cols = c()
for (cname in meta_cols){
  if (! cname %in% colnames(sobj@meta.data)){
    missing_cols = c(missing_cols, cname)
  }
}
if (length(missing_cols) > 0){
  stop(g("ERROR: The following columns are missing from metadata: {paste(missing_cols, collapse=', ')}"))
}

# remove cells from donors with fewer than MIN_NUM_CELLS_PER_ID
cell_counts_per_id = table(sobj@meta.data[[ID_COL]])
num_cells_before = ncol(sobj)
num_donors_before = length(cell_counts_per_id)
valid_ids = names(cell_counts_per_id[cell_counts_per_id >= MIN_NUM_CELLS_PER_ID])
keep <- sobj@meta.data[[ID_COL]] %in% valid_ids
sobj <- sobj[, keep]
# warn the user how many cells and donors were removed out of original
num_removed_cells = num_cells_before - ncol(sobj)
num_removed_donors = num_donors_before - length(valid_ids)
message(g("Removed {num_removed_cells} cells from {num_removed_donors} donors with fewer than {MIN_NUM_CELLS_PER_ID} cells."))

print ("**************** PREPARING DATA FOR NEBULA ****************")
# nebula 1.5.5's scToNeb() calls GetAssayData(slot = "counts"), which is defunct in SeuratObject >= 5.4.
# scToNeb_seurat() is its Seurat branch verbatim, with layer = "counts" instead of slot = "counts".
scToNeb_seurat = function(obj, assay, id, pred, offset) {
    Seurat::DefaultAssay(obj) <- assay
    covs = colnames(obj@meta.data)
    p_df <- list()
    for (k in pred) {
        if (k %in% covs) {
            p_df[k] <- obj[[k]]
        } else {
            stop("A specified variable is not available in the meta data of the Seurat object.")
        }
    }
    p_df <- data.frame(p_df, row.names = colnames(obj))
    data_neb <- list(count = Seurat::GetAssayData(obj, layer = "counts"), pred = p_df)
    if (!(id %in% covs)) stop("The sample ID is not available in the meta data of the Seurat object.")
    data_neb$id <- as.character(obj[[id]][[id]])
    if (!(offset %in% covs)) stop("The offset is not available in the meta data of the Seurat object.")
    data_neb$offset <- as.double(obj[[offset]][[offset]])
    data_neb
}
seuratdata = scToNeb_seurat(obj = sobj, assay = "RNA", id = ID_COL, pred = COVS[!grepl(":|\\*", COVS)], offset=OFFSET_COL)

df = model.matrix(
  as.formula(g("~{paste(COVS, collapse='+')}")), 
  data=seuratdata$pred)

data_g = group_cell(
  count=seuratdata$count,
  id=seuratdata$id,
  pred=df,
  offset=seuratdata$offset)


if (is.null(data_g) || is.null(data_g$count)) {
    message("NEBULA: bypassing group_cell(); using single-cell inputs.")
    data_g = list(
        count  = seuratdata$count,
        id     = factor(seuratdata$id),     # ensure factor
        pred   = as.matrix(df),             # ensure matrix
        offset = as.numeric(seuratdata$offset)
    )
}

print ("**************** RUNNING NEBULA LN ****************")

empty_nebula_result = function() {
    list(
        summary        = data.frame(),   # 0x0; combine script bind_rows() will handle it
        overdispersion = data.frame(),
        convergence    = numeric(0),
        algorithm      = character(0),
        covariance     = NULL,
        random_effect  = NULL
    )
}

re_ln = tryCatch({
    nebula(
        count  = data_g$count,
        id     = data_g$id,
        pred   = data_g$pred,
        offset = data_g$offset,
        method = 'LN',
        ncore  = N_CORES,
        output_re = TRUE,
        mincp = round(ncol(sobj) * 0.01)
    )}, error = function(e) {
    msg = conditionMessage(e)
    # fail gracefully if no genes pass filtering
    if (grepl("No gene passed the filtering", msg, ignore.case = TRUE)) {
        message("NEBULA: No genes passed the filtering in this chunk. ",
                "Emitting empty results and exiting successfully.")
        empty_nebula_result()
    } else {
        stop(e)  # real error: keep failing fast
    }
})

print(g("**************** SAVING RESULTS to {file.path(de_dir, slogan)} ****************"))
qsave(re_ln, file=outfile)
print("**************** ALL DONE ****************")


# print ("**************** RUNNING NEBULA HL ****************")
# re_hl = nebula(
#     count=data_g$count,
#     id=data_g$id,
#     pred=data_g$pred,
#     offset=data_g$offset,
#     method='HL',
#     ncore=4)

# hlout = gsub(".qs","__nebula_hl.qs",PATH)
# print(g("**************** SAVING RESULTS to {hlout} ****************"))
# qsave(re_hl, file=hlout)




