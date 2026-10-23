#!/usr/bin/env Rscript
# Build one gene-chunk Seurat object (.qs) from the raw pieces written by h5ad_to_gene_chunks.py.
#
# The pieces of a chunk are CSR rows (cells) restricted to a gene range, which is exactly the
# dgCMatrix (genes x cells) layout: p = cumsum of per-cell counts, i = local gene index, x = counts.
# The object is created exactly as chunk-sobj-by-genes.R does:
#   CreateSeuratObject(counts = <genes x all cells>, meta.data = {id, covs, nUMI_nebula}, assay = "RNA")
# with factor columns stored as character.

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

suppressMessages(suppressWarnings({
    library(getopt)
    library(Matrix)
    library(Seurat)
    library(qs)
    library(data.table)
    library(jsonlite)
}))

spec = matrix(c(
    "raw-dir", "r", 1, "character",
    "chunk",   "c", 1, "integer",
    "outdir",  "o", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
if (is.null(opt[["raw-dir"]]) || is.null(opt[["chunk"]]) || is.null(opt[["outdir"]])) {
    stop("Required: --raw-dir --chunk --outdir")
}
RAW_DIR = opt[["raw-dir"]]
CHUNK = opt[["chunk"]]

summary = fromJSON(file.path(RAW_DIR, "summary.json"))
N_FOLDS = summary$n_folds
out_qs = file.path(opt[["outdir"]], sprintf("chunk_%03d_of_%03d.qs", CHUNK, N_FOLDS))
if (file.exists(out_qs)) {
    message("Skipping ", out_qs, " (already exists)")
    quit(status = 0, save = "no")
}
chunk_dir = file.path(RAW_DIR, sprintf("chunk_%03d", CHUNK))

# ---- metadata, coerced to the R types recorded by the exporter
dtypes = fromJSON(file.path(RAW_DIR, "meta_dtypes.json"))
col_classes = c(cell = "character", unlist(dtypes))
meta = fread(file.path(RAW_DIR, "meta.tsv"), sep = "\t", colClasses = col_classes,
             na.strings = "NA", data.table = FALSE, showProgress = FALSE)
rownames(meta) = meta$cell
meta$cell = NULL
cells = rownames(meta)
n_cells = length(cells)
genes = readLines(file.path(chunk_dir, "genes.txt"))

# ---- assemble dgCMatrix from cell-block pieces (blocks are in row order)
blocks = sort(sub("\\.n\\.bin$", "", list.files(chunk_dir, pattern = "^block_[0-9]+\\.n\\.bin$")))
if (length(blocks) != summary$n_blocks) {
    stop(sprintf("chunk %d: found %d blocks, expected %d", CHUNK, length(blocks), summary$n_blocks))
}
read_piece = function(stem, ext, what, size) {
    f = file.path(chunk_dir, paste0(stem, ".", ext, ".bin"))
    readBin(f, what = what, n = file.size(f) / size, size = size)
}
n_per_cell = unlist(lapply(blocks, read_piece, ext = "n", what = "integer", size = 4L), use.names = FALSE)
i_local    = unlist(lapply(blocks, read_piece, ext = "i", what = "integer", size = 4L), use.names = FALSE)
x_counts   = unlist(lapply(blocks, read_piece, ext = "x", what = "numeric", size = 4L), use.names = FALSE)

expected_nnz = summary$chunk_nnz[CHUNK]
if (length(n_per_cell) != n_cells) stop(sprintf("per-cell counts: %d != %d cells", length(n_per_cell), n_cells))
if (length(i_local) != expected_nnz || length(x_counts) != expected_nnz) {
    stop(sprintf("chunk %d nnz mismatch: i=%d x=%d expected=%.0f", CHUNK, length(i_local), length(x_counts), expected_nnz))
}
p = c(0L, cumsum(n_per_cell))
if (anyNA(p)) stop("integer overflow building column pointers")

counts = new("dgCMatrix", i = i_local, p = p, x = x_counts,
             Dim = c(length(genes), n_cells), Dimnames = list(genes, cells))
rm(i_local, x_counts, n_per_cell, p); invisible(gc())
invisible(validObject(counts))

sobj = CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")
dir.create(opt[["outdir"]], showWarnings = FALSE, recursive = TRUE)
tmp_qs = paste0(out_qs, ".tmp")
qsave(sobj, tmp_qs)
invisible(file.rename(tmp_qs, out_qs))
message(sprintf("Wrote %s (%d genes x %d cells, nnz=%.0f)", out_qs, nrow(sobj), ncol(sobj), expected_nnz))
