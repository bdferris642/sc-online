########################################## IMPORTS ###########################################

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


print("Loading libraries...")

suppressWarnings(suppressMessages(source("~/code/sconline_code.R")))
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(RhpcBLASctl)))
suppressWarnings(suppressMessages(library(SingleCellExperiment)))
suppressWarnings(suppressMessages(library(tidyr)))

suppressWarnings(suppressMessages(source("/home/ferris/code/sconline_code.R"))) # <-- TODO inherit from de.R
suppressWarnings(suppressMessages(source("../utils.R")))
suppressWarnings(suppressMessages(source("../de.R")))
########################################### ARGUMENTS & CONSTANTS ###########################################

spec <- matrix(c(
    'path', 'p', 1, "character",
    'bkg-gene-pct', 'bgp', 1, "numeric",
    'bkg-gene-count', 'bgc', 1, "numeric",
    'contrast-col', 'cnc', 1, "character",
    'cluster-col', 'clc', 1, "character",
    'cov-list', 'cl', 1, "character",
    'rand-var', 'rv', 1, "character",
    'de-method', 'dem', 1, "character",
    'only-clusters', 'o', 1, 'character',
    'suffix', 's', 1, 'character',
    'sva-ctr-cols', 'scs', 1, 'character',
    'n-svs', 'ns', 1, 'integer',
    'num-threads', 'n', 1, 'integer',
    'calc-purity', 'cp', 1, 'logical'
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print(spec)
print(opt)

if (! is.null(opt[['num-threads']])){
    blas_set_num_threads(opt[['num-threads']])
}

BKG_GENE_COUNT_THR = ifelse(
    is.null(opt[['bkg-gene-count']]), 
    10, 
    opt[['bkg-gene-count']])
print(paste('bkg-gene-count:', BKG_GENE_COUNT_THR))

BKG_GENE_PCT_THR = ifelse(
    is.null(opt[['bkg-gene-pct']]), 
    0.01, 
    opt[['bkg-gene-pct']])
print(paste('bkg-gene-pct:', BKG_GENE_PCT_THR))

RAND_VAR = ifelse(
    is.null(opt[['rand-var']]), 
    "participant_id", 
    opt[['rand-var']])
print(paste('rand-var:', RAND_VAR))
if (tolower(RAND_VAR) %in% c("null", "na", "none")){
    RAND_VAR = NULL
}

if( is.null(opt[['contrast-col']])) {
    CONTRAST_COL = NULL
} else {
    CONTRAST_COL = opt[['contrast-col']]
}
print(paste('contrast-col:', CONTRAST_COL))

if (! is.null(opt[['cluster-col']])){
    CLUSTER_COL = opt[['cluster-col']]
} else {
    CLUSTER_COL = NULL
}

print(paste('cluster-col:', CLUSTER_COL))

if (is.null(opt[['cov-list']])) {
    COV_LIST = c('case_control', 'sex', 'age', 'QC_Gene_total_log', 'pct_mito', "pct_intronic")
} else {
    COV_LIST = strsplit(opt[['cov-list']], ",")[[1]]
}

if (! is.null(CONTRAST_COL) && !CONTRAST_COL %in% COV_LIST) {
    print("CODE ASSUMES CONTRAST COL IS IN COV LIST")
    print("Adding contrast col to cov list")
    COV_LIST = c(COV_LIST, CONTRAST_COL)
}

DE_METHOD = ifelse(
    is.null(opt[['de-method']]),
    "Trend",
    opt[['de-method']])
    
if (is.null(opt[['only-clusters']])) {
    ONLY_CLUSTERS = NULL
} else {
    ONLY_CLUSTERS = strsplit(opt[['only-clusters']], ",")[[1]]
}

if (is.null(opt[['suffix']])){
    SUFFIX = ""
} else {
    SUFFIX = opt[['suffix']]
}
SUFFIX = paste0('__', SUFFIX, format(Sys.Date(), "__%Y%m%d"))

CALC_PURITY = ifelse(
    is.null(opt[['calc-purity']]),
    TRUE,
    opt[['calc-purity']]
)

if (!is.null(opt[['n-svs']])) {
    RUN_SVA = TRUE
    N_SVS = opt[['n-svs']]
    if (! is.null(opt[['sva-ctr-cols']])) {
        SVA_CTR_COLS = strsplit(opt[['sva-ctr-cols']], ",")[[1]]
    } else {
        SVA_CTR_COLS = NULL
    }
} else {
    RUN_SVA = FALSE
    SVA_CTR_COLS = NULL
    N_SVS = NULL
}

PATH = opt[['path']]
BASE_PATH = dirname(PATH)
PSEUDOCELL_BASENAME = basename(PATH)


# hard code these for now
QUANTILE_NORM=T
MIN_NUM_PSEUDOCELLS=length(COV_LIST)+2


########################################### MAIN ###########################################
print("Reading pseudocells and seurat object data...")
pseudocell_read_path = file.path(PATH)
pseudocells_list = qread(pseudocell_read_path)

# downstream code expects a list of pseudocell sce; 
# if the object is just sce, make it a list
if (class(pseudocells_list) == "SingleCellExperiment"){
    if (is.null(CLUSTER_COL)){
        name = "all"
    } else {
        name = as.data.frame(colData(pseudocells_list))[[CLUSTER_COL]][[1]]
    }

    pseudocells_list = list(pseudocells_list)
    names(pseudocells_list) = name
}


de_dir = file.path(BASE_PATH, "de", gsub(".qs", "", PSEUDOCELL_BASENAME))
if (!dir.exists(de_dir)) {
    print(paste("Creating directory", de_dir))
    dir.create(de_dir, recursive=TRUE)
}

# Dynamically build the contrast string based on the input parameters
# we want ctr, ctrl, cntrl, or control to be the second category
# note that if both categories are in the above list, whichever comes first will be category 2

data_sce = pseudocells_list[[1]]
if (! is.null(CONTRAST_COL)){
    categories = sort(unique(colData(data_sce)[[CONTRAST_COL]]))
    category1 = categories[[1]]
    category2 = categories[[2]]
    for (c in categories) {
        if (tolower(c) %in% c("ctr", "ctrl", "cntrl", "control", "GTEx")) {
            category2 = c
            category1 = categories[categories != c][[1]]
            break
        }
    }
    contrast_str = sprintf("%s%s - %s%s", CONTRAST_COL, category1, CONTRAST_COL, category2)
    print(paste("Contrast string:", contrast_str))
}

print(paste(
    "Performing DE analysis on clusters:", 
    paste(names(pseudocells_list), collapse=", "), 
    "from cluster column", CLUSTER_COL, "..."))
for (x_name in names(pseudocells_list)){

    if (! is.null(ONLY_CLUSTERS)) {
        print(paste("checking if", x_name, "in", ONLY_CLUSTERS))
        if (! x_name %in% ONLY_CLUSTERS) {
            print(paste("Skipping", x_name))
            next
        }
    }

    print(x_name)
    x=pseudocells_list[[x_name]]

    # convert cols to numeric if need be
    for (col in c("age", "pmi_interp")){
        if (col %in% colnames(colData(x))){
            colData(x)[[col]] = as.numeric(colData(x)[[col]])
        }
    }

    if (RUN_SVA){
        # remove pre-existing cols starting with SV
        cd = as.data.frame(colData(x))
        cd = cd[, !grepl("^SV", toupper(colnames(cd)))]

        # run SVA
        
        cd = get_df_with_svs(
            edata=as.matrix(counts(x)),
            df=cd,
            cols=COV_LIST,
            ctr_cols=SVA_CTR_COLS,
            n=N_SVS)
        

        # and we want to add the svs to the COV_LIST
        SV_COLS = colnames(cd)[grepl("^SV", toupper(colnames(cd)))]
        COV_LIST = c(COV_LIST, SV_COLS)
        
        # store the columns used to calculated SVs and the Null model
        cd$SVA_cols = paste0(c(CONTRAST_COL, COV_LIST, SV_COLS), collapse=",")
        cd$SVA_ctr_cols = paste0(SVA_CTR_COLS, collapse=",")
        colData(x) = DataFrame(cd)
    }

    if (! is.null(CONTRAST_COL)){
        if (! all(c(CONTRAST_COL, COV_LIST, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x)))){
            print("Missing columns in pseudocell data:")
            print(c(CONTRAST_COL, COV_LIST, RAND_VAR, CLUSTER_COL)[! c(CONTRAST_COL, COV_LIST, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x))])
            stop("Missing columns in pseudocell data! Exiting!")
            next
        }
        cd = colData(x)[, c(CONTRAST_COL, COV_LIST, RAND_VAR, CLUSTER_COL)]
        # todo: fix this, parametrized
        x = x[, complete.cases(cd) & cd[[CONTRAST_COL]] != ""]
    } else {
        if (! all(c(COV_LIST, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x)))){
            print("Missing columns in pseudocell data:")
            print(c(COV_LIST, RAND_VAR, CLUSTER_COL)[! c(COV_LIST, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x))])
            stop("Missing columns in pseudocell data! Exiting!")
            next
        }
        cd = colData(x)[, c(COV_LIST, RAND_VAR, CLUSTER_COL)]
        x = x[, complete.cases(cd)]
        print(dim(x))
        print(colnames(cd))
    }

    print(paste0("num pseudocells: ", ncol(x)))

    res=NULL
    res_pd=NULL

    if(ncol(x)>=MIN_NUM_PSEUDOCELLS){
        print(paste0("Calculating DE genes for ", x_name))
        print(dim(x))
        
        if(!is.null(BKG_GENE_COUNT_THR)){
            
            if (! is.null(CLUSTER_COL)){
                anno = as.data.frame(colData(x))[[CLUSTER_COL]]
                tmp_bkg_genes=counts(x)[,which(anno==x_name)]
                tmp_bkg_genes_counts=rowSums(tmp_bkg_genes>0)
                tmp_bkg_genes_frac=tmp_bkg_genes_counts/sum(anno==x_name)
                tmp_bkg_genes=row.names(x)[tmp_bkg_genes_frac>=BKG_GENE_PCT_THR&tmp_bkg_genes_counts>=BKG_GENE_COUNT_THR]
            } else {
                tmp_bkg_genes=counts(x)
                tmp_bkg_genes_counts=rowSums(tmp_bkg_genes>0)
                tmp_bkg_genes_frac=tmp_bkg_genes_counts/ncol(x)
                tmp_bkg_genes=row.names(x)[tmp_bkg_genes_frac>=BKG_GENE_PCT_THR&tmp_bkg_genes_counts>=BKG_GENE_COUNT_THR]
            }
        } else {
            tmp_bkg_genes=NULL
        }

        print(paste("Number of Background genes", length(tmp_bkg_genes)))
        tryCatch({
            res=.sconline.fitLimmaFn(
                inputExpData=x,
                DEmethod=DE_METHOD,
                covariates=COV_LIST,
                randomEffect=RAND_VAR,
                bkg_genes=tmp_bkg_genes,
                quantile.norm=QUANTILE_NORM,
                prior.count=1) #TODO: check if this is the right value for prior.count
            #check the dc object, the usual consensus.correlation that I get is in the range of ~0.2 or above if rand=T
            
            expr = res[["normData"]]$E  # ngene × nsamp
            ave_expr_vec = rowMeans(expr)
            write_path = file.path(de_dir, paste0(tolower(DE_METHOD), '__', x_name, SUFFIX, '.qs'))

            if (! is.null(CONTRAST_COL)){
                contr = makeContrasts(contrasts=contrast_str, levels=res$model)

                fit2=contrasts.fit(res$fit,contrasts=contr)
            
                # TODO: Explore setting robust=F
                fit2=eBayes(fit2, robust=T, trend=T)

                res_pd = topTable(
                    fit2,
                    number=dim(fit2)[1],
                    adjust.method="BH",
                    coef=contrast_str)

                res_pd = res_pd[order(-res_pd$logFC),]
                res_pd$gene = rownames(res_pd)
                res_pd$AveExpr = ave_expr_vec[res_pd$gene]
                if(CALC_PURITY){
                    res_pd = get_residual_corrected_purity_score(
                        res=res,
                        col=CONTRAST_COL,
                        case_col=paste0(CONTRAST_COL, category1),
                        coef_df=res_pd,
                        categorical=TRUE
                    )
                } else {
                    res_pd$rank_purity = NA
                }

                res_pd = res_pd[, c("gene", "logFC", "adj.P.Val", "rank_purity", "AveExpr", "P.Value", "t", "B")] %>% arrange(desc(logFC))
                csv_write_path = gsub(".qs", paste0("__", CONTRAST_COL, ".csv"), write_path)
                write.csv(res_pd, csv_write_path, row.names=FALSE)

            } else { 
                res_pd = "NONE"
            }

            # also make top_tables for other coefs
            fit2=eBayes(res$fit, robust=T, trend=T)

            output_list = list(
                'res' = res,
                'res_pd' = res_pd,
                'pseudocells' = pseudocell_read_path, 
                'cluster_col' = CLUSTER_COL,
                'cluster' = x_name,
                'de_method' = DE_METHOD,
                'contrast_col' = CONTRAST_COL,
                'cov_list' = COV_LIST,
                'rand_var' = RAND_VAR,
                'quantile_norm' = QUANTILE_NORM,
                'min_num_pseudocells' = MIN_NUM_PSEUDOCELLS,
                'sva_ctr_cols' = SVA_CTR_COLS,
                'n_svs' = N_SVS
            )

            for(coef in colnames(res$fit$coefficients)){
                tt = topTable(
                    fit2,
                    number=dim(fit2)[1],
                    adjust.method="BH",
                    coef=coef)

                tt = tt[order(-tt$logFC),]
                tt$gene = rownames(tt)
                tt$AveExpr = ave_expr_vec[tt$gene]

                if(CALC_PURITY){
                    tt = get_residual_corrected_purity_score(
                        res=res,
                        col=coef,
                        coef_df=tt,
                        categorical=FALSE)
                } else {
                    tt$rank_purity = NA
                }
                

                tt = tt[, c("gene", "logFC", "adj.P.Val", "rank_purity", "AveExpr", "P.Value", "t", "B")] %>% arrange(desc(logFC))
                write.csv(tt, gsub(".qs", paste0("__", coef, ".csv"), write_path))
                output_list[[paste0('tt_', coef)]] = tt 
            }

            qsave(output_list, write_path)
            

        }, error=function(e){
            print(paste("Error in DE analysis for", x_name))
            print(e)
        })
    }
}
