# ---- thread hygiene (do this first) ----
threads <- 1L

# 1) BLAS/LAPACK (vendor-agnostic env vars)
Sys.setenv(
  OMP_NUM_THREADS        = threads,   # OpenMP default
  OPENBLAS_NUM_THREADS   = threads,   # OpenBLAS
  MKL_NUM_THREADS        = threads,   # Intel MKL
  VECLIB_MAXIMUM_THREADS = threads,   # Apple Accelerate
  BLIS_NUM_THREADS       = threads    # BLIS (if used)
)

# 2) If RhpcBLASctl is available, also tell R's BLAS
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(threads)
  RhpcBLASctl::omp_set_num_threads(threads)
}

# 3) Popular R packages with their own pools
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(threads)
}
if (requireNamespace("RcppParallel", quietly = TRUE)) {
  RcppParallel::setThreadOptions(numThreads = threads)  # controls TBB users (e.g., RcppParallel)
}
if (requireNamespace("fst", quietly = TRUE)) {
  fst::threads_fst(threads)
}
if (requireNamespace("arrow", quietly = TRUE)) {
  options(arrow.use_threads = (threads > 1), arrow.num_threads = threads)
}

# 4) R-level backends: keep sequential unless you explicitly use them
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

# 5) Package-specific knobs
# nebula: pass ncore=threads (you already do -n 1 at callsite)
# Seurat: if you ever use future.apply with Seurat, keeping plan(sequential) avoids nested parallelism


# Detect script path when running via Rscript
args = commandArgs(trailingOnly = FALSE)
script_path = sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir = dirname(normalizePath(script_path))
  message("Running script from: ", script_dir)
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
    library(RhpcBLASctl)
    library(Seurat)
    library(UCell)
    source(file.path(script_dir, "../utils.R"))
}))


spec <- matrix(c(
    'sobj-path', 'p', 1, "character",
    'score-path', 's', 1, "character",
    'assay', 'a', 1, "character",
    'num-threads', 'n', 1, "integer"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)

SOBJ_PATH = opt[["sobj-path"]]
SCORE_PATH = opt[["score-path"]]
if(is.null(opt[['assay']])) {
    ASSAY = "RNA"
} else {
    ASSAY = opt[["assay"]]
}

blas_set_num_threads(opt[['num-threads']])


print( paste0("Reading gene scores from ", SCORE_PATH))
score_list_df = read.table(SCORE_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
print(head(score_list_df))
gene_sets = lapply(score_list_df$genes, function(genes) {
    strsplit(genes, ",")[[1]]
})
names(gene_sets) = score_list_df$score_name

print( paste("Reading Seurat object from", SOBJ_PATH))
sobj = load_obj(SOBJ_PATH)

for (name in names(gene_sets)){
    gene_sets[[name]] = gene_sets[[name]][gene_sets[[name]] %in% rownames(GetAssayData(sobj, slot="counts", assay = ASSAY))]
    print(paste0("Calculating score for ", name))
    print("Genes:")
    print(paste(gene_sets[[name]][1:5], collapse=', '))
    print("...")
}
print(sapply(gene_sets, length))

missing_genes <- setdiff(unlist(gene_sets), rownames(GetAssayData(sobj, slot="counts", assay = ASSAY)))
if (length(missing_genes) > 0) {
    print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
}

print("Calculating U Scores")
u_scores = ScoreSignatures_UCell(
    as.matrix(GetAssayData(sobj, slot="counts", assay = ASSAY)),
    gene_sets, 
    assay = ASSAY)

u_scores = u_scores[match(rownames(sobj@meta.data), rownames(u_scores)),]
u_scores = as.data.frame(u_scores)

if (all(rownames(u_scores) == rownames(sobj@meta.data))){
    for (score_name in colnames(u_scores)){
        sobj@meta.data[[paste0(score_name, "__", ASSAY)]] = u_scores[[score_name]]
    }
}

print(paste0("Saving Seurat object to ", SOBJ_PATH))
save_obj(sobj, SOBJ_PATH)