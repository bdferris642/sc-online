source("~/code/sconline_code.R")
########################################## IMPORTS ###########################################

print("Loading libraries...")

library(getopt)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(tidyr)

source("/home/ferris/code/sconline_code.R") # <-- TODO inherit from de.R
source("/home/ferris/sc-online/utils.R")



########################################### ARGUMENTS & CONSTANTS ###########################################

spec <- matrix(c(
    'path', 'p', 1, "character",
    'bkg-gene-pct', 'bgp', 1, "numeric",
    'bkg-gene-count', 'bgc', 1, "numeric",
    'contrast-col', 'cnc', 1, "character",
    'cluster-col', 'clc', 1, "character",
    'cov-list', 'cl', 1, "character",
    'covs-from-file', 'cff', 1, "logical",
    'rand-var', 'rv', 1, "character",
    'de-method', 'dem', 1, "character",
    'only-clusters', 'o', 1, 'character',
    'suffix', 's', 1, 'character'
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

# 'base-path', 'bp', 1, "character",
#     'pseudocell-basename', 'pbn', 1, 'character',

print(spec)
print(opt)

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

CONTRAST_COL = ifelse(
    is.null(opt[['contrast-col']]), 
    "case_control", 
    opt[['contrast-col']])
print(paste('contrast-col:', CONTRAST_COL))
if (tolower(CONTRAST_COL) %in% c("null", "na", "none")){
    CONTRAST_COL = NULL
}

CLUSTER_COL = ifelse(
    is.null(opt[['cluster-col']]), 
    "cell_class", 
    opt[['cluster-col']])
print(paste('cluster-col:', CLUSTER_COL))

if (is.null(opt[['cov-list']])) {
    COV_LIST = c('case_control', 'sex', 'age', 'QC_Gene_total_log', 'pct_mito', "pct_intronic")
} else {
    COV_LIST = strsplit(opt[['cov-list']], ",")[[1]]
}

if (is.null(opt[['covs-from-file']])) {
    COV_FROM_FILE = FALSE
} else {
    COV_FROM_FILE = opt[['covs-from-file']]
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

SUFFIX = ifelse(
    is.null(opt[['suffix']]),
    format(Sys.Date(), "__%Y%m%d"), 
    paste0('__', opt[['suffix']])
)

PATH = opt[['path']]
BASE_PATH = dirname(PATH)
PSEUDOCELL_BASENAME = basename(PATH)

# to get the seurat basename, take the pseudocell basename until the _pseudocell
# and remove everything that comes after it as well. Finally, append .qs
SEURAT_BASENAME = gsub("_pseudocell.*", ".qs", PSEUDOCELL_BASENAME)


# hard code these for now
QUANTILE_NORM=T
MIN_NUM_PSEUDOCELLS=length(COV_LIST)+2


########################################### MAIN ###########################################
print("Reading pseudocells and seurat object data...")
pseudocell_read_path = file.path(BASE_PATH, PSEUDOCELL_BASENAME)
s_obj_read_path = file.path(BASE_PATH, SEURAT_BASENAME)

pseudocells_list = qread(pseudocell_read_path)
s_obj = qread(s_obj_read_path)

# TODO! Fix for now
#s_obj = s_obj[, s_obj@meta.data[[CONTRAST_COL]] %in% c("pd", "ctr")]
if (! is.null(CONTRAST_COL)){
    print(table(s_obj@meta.data[[CONTRAST_COL]]))
    print(table(s_obj@meta.data[[CLUSTER_COL]]))
    print(table(
        s_obj@meta.data[[CONTRAST_COL]],
        s_obj@meta.data[[CLUSTER_COL]]))
}


data_sce = as.SingleCellExperiment(s_obj)
de_dir = file.path(BASE_PATH, "de", gsub(".qs", "", PSEUDOCELL_BASENAME))
if (!dir.exists(de_dir)) {
    print(paste("Creating directory", de_dir))
    dir.create(de_dir, recursive=TRUE)
}

# Dynamically build the contrast string based on the input parameters
# we want ctr, ctrl, cntrl, or control to be the second category
# note that if both categories are in the above list, whichever comes first will be category 2
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

print(paste("Performing DE analysis on clusters:", paste(names(pseudocells_list), collapse=", "), "from cluster column", CLUSTER_COL, "..."))
for (x_name in names(pseudocells_list)){

    if (! is.null(ONLY_CLUSTERS)) {
        print(paste("checking if", x_name, "in", ONLY_CLUSTERS))
        if (! x_name %in% ONLY_CLUSTERS) {
            print(paste("Skipping", x_name))
            next
        }
    }

    if (COV_FROM_FILE){
        # the file is a list of covariates to use, indexed by a cell type, with value being a comma-delimited string
        file_de_list = qread(sub(".qs$", "__covariates_used_incl_peer.qs", pseudocell_read_path))
        if (x_name %in% names(file_de_list)){
            this_cov_list = strsplit(file_de_list[[x_name]], ",")[[1]]
        } else {
            cat(paste0("No covariates found in file for cell type ", x_name, "\nFalling  back to default covariates"))
            this_cov_list = COV_LIST
        }
    } else {
        this_cov_list = COV_LIST
    }
    print(paste0('COVARIATE LIST for ct ', x_name))
    print(this_cov_list)

    print(x_name)
    x=pseudocells_list[[x_name]]

    if (! is.null(CONTRAST_COL)){
        if (! all(c(CONTRAST_COL, this_cov_list, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x)))){
            print("Missing columns in pseudocell data:")
            print(c(CONTRAST_COL, this_cov_list, RAND_VAR, CLUSTER_COL)[! c(CONTRAST_COL, this_cov_list, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x))])
            stop("Missing columns in pseudocell data! Exiting!")
            next
        }
        cd = colData(x)[, c(CONTRAST_COL, this_cov_list, RAND_VAR, CLUSTER_COL)]
        # TODO fix this
        # if (length(unique(cd[[CONTRAST_COL]] != 2))) {
        #     print("Contrast column must have exactly 2 categories! Skipping...")
        #     next
        # }
    
        # todo: fix this, parametrized
        x = x[, complete.cases(cd) & cd[[CONTRAST_COL]] != ""]
    } else {
        if (! all(c(this_cov_list, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x)))){
            print("Missing columns in pseudocell data:")
            print(c(this_cov_list, RAND_VAR, CLUSTER_COL)[! c(this_cov_list, RAND_VAR, CLUSTER_COL) %in% colnames(colData(x))])
            stop("Missing columns in pseudocell data! Exiting!")
            next
        }
        cd = colData(x)[, c(this_cov_list, RAND_VAR, CLUSTER_COL)]
        x = x[, complete.cases(cd)]
    }

    print(paste0("num pseudocells: ", ncol(x)))

    res=NULL
    res_pd=NULL


    if(ncol(x)>=MIN_NUM_PSEUDOCELLS){
        print(paste0("Calculating DE genes for ", x_name))
        print(dim(x))
        
        if(!is.null(BKG_GENE_COUNT_THR)){
            anno = as.data.frame(colData(data_sce))[[CLUSTER_COL]]
            #background genes should be selectd based on the cell level expression data
            tmp_bkg_genes=counts(data_sce)[,which(anno==x_name)]
            tmp_bkg_genes_counts=rowSums(tmp_bkg_genes>0)
            tmp_bkg_genes_frac=tmp_bkg_genes_counts/sum(anno==x_name)
            tmp_bkg_genes=row.names(data_sce)[tmp_bkg_genes_frac>=BKG_GENE_PCT_THR&tmp_bkg_genes_counts>=BKG_GENE_COUNT_THR]
        } else {
            tmp_bkg_genes=NULL
        }

        print(paste("Number of Background genes", length(tmp_bkg_genes)))
        
        #inputExpData=x;covariates=COV_LIST;randomEffect=rand_effect;bkg_genes=tmp_bkg_genes;quantile.norm=quantile_norm;prior.count=1
        tryCatch({
            res=.sconline.fitLimmaFn(
                inputExpData=x,
                DEmethod=DE_METHOD,
                covariates=this_cov_list,
                randomEffect=RAND_VAR,
                bkg_genes=tmp_bkg_genes,
                quantile.norm=QUANTILE_NORM,
                prior.count=1) #TODO: check if this is the right value for prior.count
            #check the dc object, the usual consensus.correlation that I get is in the range of ~0.2 or above if rand=T
            
            write_path = file.path(de_dir, paste0(tolower(DE_METHOD), '__', x_name, SUFFIX, '.qs'))
            # TODO: parametrize this!

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

                
                csv_write_path = gsub(".qs", paste0("__", CONTRAST_COL, ".csv"), write_path)
                write.csv(res_pd, csv_write_path, row.names=FALSE)
            } else { 
                res_pd = NULL
            }

            # also make top_tables for other coefs
            fit2=eBayes(res$fit, robust=T, trend=T)
                    
            # tt_age = topTable(
            #     fit2,
            #     number=dim(fit2)[1],
            #     adjust.method="BH",
            #     coef="age")

            # tt_age = tt_age[order(-tt_age$logFC),]
            # tt_age$gene = rownames(tt_age)

            # write.csv(tt_age, gsub(".qs", "__age.csv", write_path))

            # tt_sex = topTable(
            #     fit2,
            #     number=dim(fit2)[1],
            #     adjust.method="BH",
            #     coef="sexMale")
            # write.csv(tt_sex, gsub(".qs", "__sex.csv", write_path))

            output_list = list(
                'res' = res,
                'res_pd' = res_pd,
                # 'res_sex' = tt_sex,
                # 'res_age' = tt_age,
                'seurat_object' = s_obj_read_path,
                'pseudocells' = pseudocell_read_path, 
                'cluster_col' = CLUSTER_COL,
                'cluster' = x_name,
                'de_method' = DE_METHOD,
                'contrast_col' = CONTRAST_COL,
                'cov_list' = this_cov_list,
                'rand_var' = RAND_VAR,
                'quantile_norm' = QUANTILE_NORM,
                'min_num_pseudocells' = MIN_NUM_PSEUDOCELLS
            )

            for(coef in colnames(res$fit$coefficients)){
                tt = topTable(
                    fit2,
                    number=dim(fit2)[1],
                    adjust.method="BH",
                    coef=coef)

                tt = tt[order(-tt$logFC),]
                tt$gene = rownames(tt)
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
