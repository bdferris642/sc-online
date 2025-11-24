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

    source(file.path(script_dir, "../utils.R"))
    g = glue::glue
}))

spec = matrix(c(
    'path', 'p', 1, "character",
    'covs', 'c', 2, "character",
    'id-col', 'r', 2, "character",
    'offset-col', 'f', 2, "character",
    'suffix', 's', 2, "character",
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


print ("**************** PREPARING DATA FOR NEBULA ****************")
seuratdata = scToNeb(obj = sobj, assay = "RNA", id = ID_COL, pred = COVS[!grepl(":|\\*", COVS)], offset=OFFSET_COL)

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


# fail gracefully if no genes pass filtering
re_ln = tryCatch({
    nebula(
        count  = data_g$count,
        id     = data_g$id,
        pred   = data_g$pred,
        offset = data_g$offset,
        method = 'LN',
        ncore  = N_CORES,
        output_re = TRUE
    )}, error = function(e) {
    msg = conditionMessage(e)
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




