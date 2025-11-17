# Takes in paths to a 10x Gene Expression H5 and molecule_info.h5 files from CellRanger, 
# an output directory
# Along with Optional CB inputs H5 file

# Outputs a CSV file with the results of the Dropsift analysis
# Along with a PDF report of the SVM Nucleus Caller
# CB inputs are optional, but if provided, will be used to calculate the contamination fraction

################################### IMPORTS ###################################

suppressMessages(library(dplyr))
suppressMessages(library(DropSift))
suppressMessages(library(Matrix))
suppressMessages(library(optparse))
suppressMessages(library(rhdf5))
suppressMessages(library(Seurat))


################################### ARGUMENTS ###################################

option_list = list(
  make_option(c("-g", "--gex_h5_path"),   type = "character",
              help = "Path to 10x Gene Expression H5 (required).", metavar = "FILE"),
  make_option(c("-m", "--mol_info_path"), type = "character",
              help = "Path to 10x molecule_info.h5 (required).", metavar = "FILE"),
  make_option(c("-c", "--cb_h5_path"),    type = "character", default = NULL,
              help = "Path to CB inputs H5 (optional).", metavar = "FILE"),
  make_option(c("-o", "--output_dir"),    type = "character",
              help = "Directory for outputs (required).", metavar = "DIR")
)

parser = OptionParser(option_list = option_list)
opt = parse_args(parser)

# Validate required args
req = c("gex_h5_path", "mol_info_path", "output_dir")
missing = req[sapply(opt[req], function(x) is.null(x) || identical(x, ""))]
if (length(missing)) {
    print_help(parser)
    stop(sprintf(
        "Missing required arguments: %s",
        paste(paste0("--", gsub("_", "-", missing)), collapse = ", ")))
}

MOL_INFO_PATH = opt$mol_info_path
GEX_H5_PATH   = opt$gex_h5_path
OUTPUT_DIR    = opt$output_dir

USE_CB_INPUTS = !is.null(opt$cb_h5_path)
CB_H5_PATH    = if (USE_CB_INPUTS) opt$cb_h5_path else NULL

dropsift_input_cols = c("cell_barcode", "num_transcripts", "pct_intronic", "pct_mt")
if (USE_CB_INPUTS){
    dropsift_input_cols = c(dropsift_input_cols, "frac_contamination")
}

# Optional: ensure output dir exists
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################### CONSTANTS ###################################

# Minimum counts threshold for a cell to be considered a cell
COUNTS_MIN = 10


################################### MAIN ###################################

# Load the raw Gene Expression counts matrix to get:
# counts matrix, along with num_transcripts and pct_mt (fraction of genes with names that start with "MT-")
# TODO: add support for ENSGs and non-human gene symbols
cr_counts = Read10X_h5(GEX_H5_PATH)

# Keep barcodes (columns) above count threshold
keep_cols = colSums(cr_counts) > COUNTS_MIN
cr_counts_filt = cr_counts[, keep_cols, drop = FALSE]
barcodes_filt  = colnames(cr_counts_filt)

num_transcripts = colSums(cr_counts_filt)
mt_rows = grepl("^MT-", rownames(cr_counts_filt), ignore.case = TRUE)
num_transcripts_mt = if (any(mt_rows)) colSums(cr_counts_filt[mt_rows, , drop = FALSE]) else rep(0, ncol(cr_counts_filt))
pct_mt = ifelse(num_transcripts > 0, num_transcripts_mt / num_transcripts, NA_real_)

counts_summary = data.frame(
    cell_barcode = colnames(cr_counts_filt),
    num_transcripts = as.numeric(num_transcripts),
    pct_mt = as.numeric(pct_mt),
    stringsAsFactors = FALSE
)

# Fraction intronic is extracted from the molecular info h5
# and is defined as the fraction of UMIs for which umi_type == 0
# (see https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-molecule-info)

mol_info_h5fetch = function(x){return(h5read(MOL_INFO_PATH, x))}
mol_df = data.frame(
    cell_barcode=mol_info_h5fetch("barcodes")[mol_info_h5fetch("barcode_idx")+1] %>% paste0("-1"),
    feature=mol_info_h5fetch("features/name")[mol_info_h5fetch("feature_idx")+1],
    umi=mol_info_h5fetch("umi"),
    umi_type=mol_info_h5fetch("umi_type"),
    count=mol_info_h5fetch("count")
)
mol_df_filt = mol_df %>% filter(cell_barcode %in% barcodes_filt)
mol_df_grouped = (mol_df_filt 
    %>% group_by(cell_barcode) 
    %>% summarize(
        n_umis=n(), 
        pct_intronic=sum(umi_type==0)/n_umis)
)

droplet_summary = merge(
    counts_summary, mol_df_grouped[,c("cell_barcode", "pct_intronic")], 
    by="cell_barcode")

if (length(droplet_summary$cell_barcode) != length(unique(droplet_summary$cell_barcode))) {
    stop("DROPSIFT FAILED: duplicate cell barcodes found in droplet_summary data frame")
}

if (any(is.na(droplet_summary))) {
    na_cols = colnames(droplet_summary)[apply(droplet_summary, 2, function(x) any(is.na(x)))]
    stop(paste("DROPSIFT FAILED: NA values found in columns:", paste(na_cols, collapse=", ")))
}

if (USE_CB_INPUTS){

    # get the frac_contamination estimate 
    # by subtracting the counts in the CB matrix from the counts in the CR matrix
    # need to construct the sparse matrix from the h5 ourselves...

    cb_h5_data = h5read(CB_H5_PATH, "/matrix")
    cb_counts = cb_h5_data$data
    cb_indices = cb_h5_data$indices
    cb_indptr = cb_h5_data$indptr
    cb_barcodes = cb_h5_data$barcodes
    cb_num_genes = length(cb_h5_data$features$name)  # Number of genes
    cb_num_cells = length(cb_barcodes)  # Number of cells

    cb_counts_matrix = sparseMatrix(i=cb_indices + 1,  # +1 for 1-based indexing in R
                                    p=cb_indptr,
                                    x=cb_counts,
                                    dims=c(cb_num_genes, cb_num_cells))
    rownames(cb_counts_matrix) = cb_h5_data$features$name 
    colnames(cb_counts_matrix) = cb_barcodes

    shared_colnames = intersect(colnames(cr_counts_filt), colnames(cb_counts_matrix))
    cr_shared = cr_counts_filt[, shared_colnames, drop= FALSE]
    cb_shared = cb_counts_matrix[, shared_colnames, drop= FALSE]

    cb_inputs = data.frame(
        cell_barcode = colnames(cb_shared),
        frac_contamination = (colSums(cr_shared) - colSums(cb_shared)) / colSums(cr_shared)
    )
    
    droplet_summary = merge(
        droplet_summary, cb_inputs, 
        by="cell_barcode")

    if (length(droplet_summary$cell_barcode) != length(unique(droplet_summary$cell_barcode))) {
        stop("DROPSIFT FAILED: duplicate cell barcodes found in droplet_summary data frame")
    }

    if (any(is.na(droplet_summary))) {
        na_cols = colnames(droplet_summary)[apply(droplet_summary, 2, function(x) any(is.na(x)))]
        stop(paste("DROPSIFT FAILED: NA values found in columns:", paste(na_cols, collapse=", ")))
    }
}

droplet_summary = droplet_summary[, dropsift_input_cols]

print("RUNNING DROPSIFT")
svmNucleusCaller = SvmNucleusCaller(
    cellFeatures=droplet_summary,
    dgeMatrix=cr_counts_filt,
    useCBRBFeatures = USE_CB_INPUTS)

dropsift_out_df = svmNucleusCaller$cell_features

print("SAVING DROPSIFT OUTPUTS")
write.csv(dropsift_out_df, file=file.path(OUTPUT_DIR,"dropsift_output.csv"), row.names=FALSE)
pdf(file.path(OUTPUT_DIR, "svmNucleusCallerReport.pdf"))
plotSvmNucleusCaller(svmNucleusCaller)
dev.off()

print("DROPSIFT DONE")