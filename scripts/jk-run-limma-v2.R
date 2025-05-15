########################################## IMPORTS ###########################################

print("Loading libraries...")

suppressWarnings(suppressMessages(source("~/code/sconline_code.R")))
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(RhpcBLASctl)))
suppressWarnings(suppressMessages(library(SingleCellExperiment)))
suppressWarnings(suppressMessages(library(tidyr)))

suppressWarnings(suppressMessages(source("/home/ferris/code/sconline_code.R"))) # <-- TODO inherit from de.R
suppressWarnings(suppressMessages(source("/home/ferris/sc-online/utils.R")))
########################################### ARGUMENTS & CONSTANTS ###########################################

spec <- matrix(c(
    'path', 'p', 1, "character",
    'num-jks', 'nj', 1, 'integer',
    'jk-frac', 'js', 1, 'integer',
    'jk-col', 'jc', 1, 'character',
    'stratification-col', 'sc', 1, 'character',
    'contrast-col', 'cnc', 1, "character",
    'cluster-col', 'clc', 1, "character",
    'cov-list', 'cl', 1, "character",
    'rand-var', 'rv', 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print(spec)
print(opt)

RAND_VAR = ifelse(
    is.null(opt[['rand-var']]), 
    "participant_id", 
    opt[['rand-var']])
message(paste('rand-var:', RAND_VAR))
if (tolower(RAND_VAR) %in% c("null", "na", "none")){
    RAND_VAR = NULL

}

if( is.null(opt[['contrast-col']])) {
    CONTRAST_COL = NULL
    CONTRAST_STR = ""
} else {
    CONTRAST_COL = opt[['contrast-col']]
    CONTRAST_STR = glue("--contrast-col={CONTRAST_COL}")
}
message(paste('contrast-col:', CONTRAST_COL))

if (! is.null(opt[['cluster-col']])){
    CLUSTER_COL = opt[['cluster-col']]
    CLUSTER_STR = glue("--cluster-col={CLUSTER_COL}")
} else {
    CLUSTER_COL = NULL
    CLUSTER_STR = ""
}

message(paste('\ncluster-col:', CLUSTER_COL))

if (is.null(opt[['cov-list']])) {
    COV_LIST = c('case_control', 'sex', 'age', 'QC_Gene_total_log', 'pct_mito', "pct_intronic")
} else {
    COV_LIST = strsplit(opt[['cov-list']], ",")[[1]]
}

if (! is.null(CONTRAST_COL) && !CONTRAST_COL %in% COV_LIST) {
    message("\nCODE ASSUMES CONTRAST COL IS IN COV LIST")
    message("\nAdding contrast col to cov list")
    COV_LIST = c(COV_LIST, CONTRAST_COL)
}

if (is.null(opt[['num-jks']])) {
    NUM_JK = 100
} else {
    NUM_JK = opt[['num-jks']]
}

if (is.null(opt[['jk-frac']])) {
    JK_FRAC = 0.5
} else {
    JK_FRAC = opt[['jk-frac']]
}

if (is.null(opt[['jk-col']])) { 
    JK_COL = "participant_id"
} else {
    JK_COL = opt[['jk-col']]
}

if (is.null(opt[['stratification-col']])) {
    STRATIFICATION_COL = NULL
} else {
    STRATIFICATION_COL = opt[['stratification-col']]
}

PATH = opt[['path']]
BASE_PATH = dirname(PATH)
PSEUDOCELL_BASENAME = basename(PATH)
SLOGAN = gsub(".qs", "", PSEUDOCELL_BASENAME)
JK_PATH = file.path(BASE_PATH, 'jk', SLOGAN)
message(glue("\nJackknife samples will be saved to {JK_PATH}"))
dir.create(JK_PATH, showWarnings = TRUE, recursive = TRUE)

########################################### MAIN ###########################################

message(glue("\nLoading data from {PATH}"))

sce = qread(PATH)

# get your list of jackknife samples
# TODO: stratify by stratification-col
jk_sample_list = list()
indiv_list = unique(as.data.frame(colData(sce))[[JK_COL]])

message(glue("\nNumber of individuals: {length(indiv_list)}"))
message(glue("\nSaving {NUM_JK} jackknife samples to {JK_PATH}"))

set.seed(642)
for (i in 1:NUM_JK) {
    message(glue("\nSaving Jackknife sample {i} of {NUM_JK}..."))
    if (is.null(STRATIFICATION_COL)){
        samples = sample(
            x = indiv_list,
            size = round(length(indiv_list) * JK_FRAC),
            replace = FALSE
        )
    } else {
        
        # Stratified sampling
        strat_df = data.frame(
            indiv = indiv_list,
            strat = colData(sce)[[STRATIFICATION_COL]][match(indiv_list, colData(sce)[[JK_COL]])]
        )
        strat_df = strat_df[!duplicated(strat_df), ]
        
        samples = unlist(lapply(split(strat_df, strat_df$strat), function(sub_df) {
            sample(
                x = sub_df$indiv,
                size = round(nrow(sub_df) * JK_FRAC),
                replace = FALSE
            )
        }))
    }
    sample_sce = sce[, colData(sce)[[JK_COL]] %in% samples]
    qsave(sample_sce, file.path(JK_PATH, paste0(SLOGAN, "_jk_", i, ".qs")))

    jk_sample_list[[i]] = paste(samples, collapse = ",")
}

# now run the jackknife samples in parallel
# must constrain number of threads (1) and number of processes (parallel jobs 36 / 40) so that the instance doesn't freeze
# the huge sed command is literally to extract the jackknife sample number from the filename
command_str = glue(
    "ls {JK_PATH} | \\
    grep {SLOGAN} | grep jk | grep qs | \\
    xargs -I @ sh -c 'f=\"@\"; s=$(echo \"$f\" | sed -n \"s/.*\\(jk_[0-9]\\+\\)\\.qs/\\1/p\"); echo Rscript ~/sc-online/scripts/run-limma-v2.R \\
    --path={JK_PATH}/$f \\
    --rand-var={RAND_VAR} \\
    --cov-list={paste(COV_LIST, collapse = \",\")} \\
    {CLUSTER_STR} {CONTRAST_STR} --num-threads=1 \\
    --suffix=$s' | parallel -j 36" 
)

writeLines(
    command_str,
    file.path(JK_PATH, paste0(SLOGAN, "_jk_command.txt"))
)

message("\nRunning jackknife samples in parallel...")
output <- system(command_str, intern = TRUE)

message("\nDone! Writing samples and output...")
lapply(jk_sample_list, write, file.path(JK_PATH, paste0(SLOGAN, "_jk_samples.txt")), append=TRUE, ncolumns=5000)


writeLines(
    output,
    file.path(JK_PATH, paste0(SLOGAN, "_jk_output.txt"))
)

# now remove the jackknife qs files in the jk folder
message("\nCleaning up jackknife files...")
unlink(
    list.files(
        path = JK_PATH,
        pattern = paste0(SLOGAN, "_jk_[0-9]+\\.qs"),
        full.names = TRUE
    )
)
