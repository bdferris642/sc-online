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
    'base-path', 'bp', 1, "character",
    'pseudocell-basename', 'pbn', 1, 'character',
    'seurat-basename', 'sbn', 1, 'character',
    'bkg-gene-pct', 'bgp', 1, "numeric",
    'bkg-gene-count', 'bgc', 1, "numeric",
    'contrast-col', 'cnc', 1, "character",
    'cluster-col', 'clc', 1, "character",
    'cov-list', 'cl', 1, "character",
    'rand-var', 'rv', 1, "character",
    'de-method', 'dem', 1, "character",
    'only-clusters', 'o', 1, 'character',
    'suffix', 's', 1, 'character'
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print(spec)
print(opt)

BKG_GENE_COUNT_THR = ifelse(
    is.null(opt[['bkg-gene-count']]), 
    10, 
    opt[['bkg-gene-count']])
print(paste('bkg-gene-count:', BKG_GENE_COUNT_THR))

BKG_GENE_COUNT_THR = ifelse(
    is.null(opt[['bkg-gene-pct']]), 
    0.01, 
    opt[['bkg-gene-pct']])
print(paste('bkg-gene-pct:', BKG_GENE_COUNT_THR))

RAND_VAR = ifelse(
    is.null(opt[['rand-var']]), 
    "participant_id", 
    opt[['rand-var']])
print(paste('rand-var:', RAND_VAR))

CONTRAST_COL = ifelse(
    is.null(opt[['contrast-col']]), 
    "case_control", 
    opt[['contrast-col']])
print(paste('contrast-col:', CONTRAST_COL))

CLUSTER_COL = ifelse(
    is.null(opt[['cluster-col']]), 
    "cell_class", 
    opt[['cluster-col']])
print(paste('cluster-col:', CLUSTER_COL))

if (is.null(opt[['cov-list']])) {
    COV_LIST = c('case_control', 'sex', 'age', 'QC_Gene_total_log', 'QC_MT.pct')
} else {
    COV_LIST = strsplit(opt[['cov-list']], ",")[[1]]
}
print('cov-list:')
print(COV_LIST)

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


if (is.null(opt[['base-path']]) || is.null(opt[['pseudocell-basename']]) || is.null(opt[['seurat-basename']])) {
    stop("Please provide a base-path, a pseudocell-basename, and a seurat-basename")
}

BASE_PATH = opt[['base-path']]
PSEUDOCELL_BASENAME = opt[['pseudocell-basename']]
SEURAT_BASENAME = opt[['seurat-basename']]


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
print(table(s_obj@meta.data[[CONTRAST_COL]]))
s_obj = s_obj[, s_obj@meta.data[[CONTRAST_COL]] %in% c("pd", "ctr")]
print(table(s_obj@meta.data[[CONTRAST_COL]]))


data_sce = as.SingleCellExperiment(s_obj)


de_dir = file.path(BASE_PATH, "de", gsub(".qs", "", PSEUDOCELL_BASENAME))
if (!dir.exists(de_dir)) {
    print(paste("Creating directory", de_dir))
    dir.create(de_dir, recursive=TRUE)
}

# Dynamically build the contrast string based on the input parameters
# we want ctr, ctrl, cntrl, or control to be the second category
# note that if both categories are in the above list, whichever comes first will be category 2
categories = sort(unique(colData(data_sce)[[CONTRAST_COL]]))
category1 = categories[[1]]
category2 = categories[[2]]
for (c in categories) {
    if (tolower(c) %in% c("ctr", "ctrl", "cntrl", "control")) {
        category2 = c
        category1 = categories[categories != c][[1]]
        break
    }
}
contrast_str = sprintf("%s%s - %s%s", CONTRAST_COL, category1, CONTRAST_COL, category2)
print(paste("Contrast string:", contrast_str))

print(paste("Performing DE analysis on clusters:", paste(names(pseudocells_list), collapse=", "), "from cluster column", CLUSTER_COL, "..."))
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
    cd = colData(x)[, c(CONTRAST_COL, COV_LIST, RAND_VAR, CLUSTER_COL)]
    
    # todo: fix this, parametrized
    x = x[, complete.cases(cd) & cd$case_control != ""]

    # if (length(unique(cd[[CONTRAST_COL]] != 2))) {
    #     print("Contrast column must have exactly 2 categories! Skipping...")
    #     next
    # }

    res=NULL
    res_pd=NULL

    if(ncol(x)>=MIN_NUM_PSEUDOCELLS){
        print(paste0("Calculating DE genes for ", x_name))
        
        if(!is.null(BKG_GENE_COUNT_THR)){
            anno = as.data.frame(colData(data_sce))[[CLUSTER_COL]]
            #background genes should be selectd based on the cell level expression data
            tmp_bkg_genes=counts(data_sce)[,which(anno==x_name)]
            tmp_bkg_genes_counts=rowSums(tmp_bkg_genes>0)
            tmp_bkg_genes_frac=tmp_bkg_genes_counts/sum(anno==x_name)
            tmp_bkg_genes=row.names(data_sce)[tmp_bkg_genes_frac>=BKG_GENE_COUNT_THR&tmp_bkg_genes_counts>=BKG_GENE_COUNT_THR]
        } else {
            tmp_bkg_genes=NULL
        }
        
        #inputExpData=x;covariates=COV_LIST;randomEffect=rand_effect;bkg_genes=tmp_bkg_genes;quantile.norm=quantile_norm;prior.count=1
        res=.sconline.fitLimmaFn(
            inputExpData=x,
            DEmethod=DE_METHOD,
            covariates=COV_LIST,
            randomEffect=RAND_VAR,
            bkg_genes=tmp_bkg_genes,
            quantile.norm=QUANTILE_NORM,
            prior.count=1) #TODO: check if this is the right value for prior.count
        #check the dc object, the usual consensus.correlation that I get is in the range of ~0.2 or above if rand=T
        
        # TODO: parametrize this!
        #contr = makeContrasts(contrast_str, levels=res$model)
        contr = makeContrasts(case_controlpd - case_controlctr, levels = res$model)
        fit2=contrasts.fit(res$fit,contrasts=contr)
        
        # TODO: Explore setting robust=F
        fit2=eBayes(fit2, robust=T, trend=T)

        res_pd = topTable(
            fit2,
            number=dim(fit2)[1],
            adjust.method="BH",
            coef=contrast_str)

        res_pd = res_pd[order(-res_pd$adj.P.Val),]
        res_pd$gene = rownames(res_pd)
        res_pd = res_pd[, c('gene', 'logFC', 'adj.P.Val', 'AveExpr')]

        output_list = list(
            'res' = res,
            'res_pd' = res_pd,
            'seurat_object' = s_obj_read_path,
            'pseudocells' = pseudocell_read_path, 
            'cluster_col' = CLUSTER_COL,
            'cluster' = x_name,
            'de_method' = DE_METHOD,
            'contrast_col' = CONTRAST_COL,
            'cov_list' = COV_LIST,
            'rand_var' = RAND_VAR,
            'quantile_norm' = QUANTILE_NORM,
            'min_num_pseudocells' = MIN_NUM_PSEUDOCELLS
        )

        write_path = file.path(de_dir, paste0(tolower(DE_METHOD), '__', x_name, SUFFIX, '.qs'))
        qsave(output_list, write_path)

        csv_write_path = gsub(".qs", ".csv", write_path)
        write.csv(res_pd, csv_write_path, row.names=FALSE)
    }
}
