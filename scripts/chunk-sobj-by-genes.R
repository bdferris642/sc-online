#!/usr/bin/env Rscript
suppressMessages(suppressWarnings({
  library(getopt)
  library(qs)
  library(Seurat)
  library(Matrix)
  library(glue)
}))

g=glue::glue
spec = matrix(c(
    'path',       'p', 1, "character",
    'id-col',     'r', 2, "character",
    'covs',       'c', 2, "character",
    'offset-col', 'f', 2, "character",
    'n-folds',    'n', 1, "integer",
    'outdir',     'o', 1, "character",
    'new-offset-colname', 'x', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

PATH_QS    = opt[['path']]
ID_COL     = if (is.null(opt[['id-col']])) "donor_id" else opt[['id-col']]
COVS_RAW   = opt[['covs']]
OFFSET_COL = if (is.null(opt[['offset-col']])) "nCount_RNA" else opt[['offset-col']]
N_FOLDS    = opt[['n-folds']]
OUTDIR     = opt[['outdir']]
if (is.null(opt[['new-offset-colname']])) {
    NEW_OFFSET_COLNAME = "nUMI_nebula"
} else {
    NEW_OFFSET_COLNAME = opt[['new-offset-colname']]
}

if (is.null(PATH_QS) || is.null(N_FOLDS) || is.null(OUTDIR)) {
    stop("Required: --path, --n-folds, --outdir")
}

message("***** LOADING sobj: ", PATH_QS)
sobj = qread(PATH_QS)

if (!OFFSET_COL %in% colnames(sobj@meta.data)) {
    stop(g("offset-col {OFFSET_COL} not in metadata."))
}

# Create nUMI_nebula from OFFSET_COL
sobj@meta.data[[NEW_OFFSET_COLNAME]] = sobj@meta.data[[OFFSET_COL]]

covs = if (!is.null(COVS_RAW) && nchar(COVS_RAW) > 0) strsplit(COVS_RAW, ",")[[1]] else character(0)
meta_keep = unique(c(ID_COL, covs, NEW_OFFSET_COLNAME))

# separately consider interactions containing : or * 
meta_keep_interactions = grepl(":|\\*", meta_keep)
meta_keep_interaction_parts = c()
for (mk in meta_keep[meta_keep_interactions]) {
    parts = unlist(strsplit(mk, ":|\\*"))
    meta_keep_interaction_parts = c(meta_keep_interaction_parts, parts)
}
meta_keep = meta_keep[!meta_keep_interactions]
meta_keep = unique(c(meta_keep, meta_keep_interaction_parts))

missing_meta = setdiff(meta_keep, colnames(sobj@meta.data))
if (length(missing_meta) > 0) {
    stop(g("Missing metadata: {missing_meta}"))
}

counts = GetAssayData(sobj, assay="RNA", slot="counts")
g_total = nrow(counts)
if (N_FOLDS < 1) stop("n-folds must be >= 1")
if (N_FOLDS > g_total) {
  warning("n-folds > #genes; reducing to #genes")
  N_FOLDS = g_total
}

chunk_sizes = rep(floor(g_total / N_FOLDS), N_FOLDS)
remainder = g_total %% N_FOLDS
if (remainder > 0) chunk_sizes[seq_len(remainder)] = chunk_sizes[seq_len(remainder)] + 1
starts = c(1, cumsum(chunk_sizes)[-length(chunk_sizes)] + 1)
ends   = cumsum(chunk_sizes)

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
message(g("Splitting {g_total} genes into {N_FOLDS} chunks"))

for (i in seq_len(N_FOLDS)) {
    out_i = file.path(OUTDIR, sprintf("chunk_%03d_of_%03d.qs", i, N_FOLDS))
    if (file.exists(out_i)) {
        message(g("Skipping {out_i} (already exists)"))
        next
    }

    g_idx = starts[i]:ends[i]
    counts_i = counts[g_idx, , drop=FALSE]
    meta_i = sobj@meta.data[, meta_keep, drop=FALSE]
    sobj_i = CreateSeuratObject(counts=counts_i, meta.data=meta_i, assay="RNA")
    
    qsave(sobj_i, out_i)
    message(g("Wrote {out_i} (genes {starts[i]} --> {ends[i]}; nGenes={nrow(counts_i)})"))
}

message("***** SPLITTING DONE.")
