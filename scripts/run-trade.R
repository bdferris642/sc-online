suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(devtools)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(ashr)))
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(qs)))
suppressWarnings(suppressMessages(library(RhpcBLASctl)))
suppressWarnings(suppressMessages(library(TRADEtools)))

suppressWarnings(suppressMessages(source("~/sc-online/utils.R")))

# path is absolute path to qs or rds file of limma results from run-limma-v2.R
# contrast-col, optional: the name of the coefficient to be used as a contrast
# contrast-str, optional: the string to be used as a contrast, e.g. "conditionA - conditionB"

spec <- matrix(c(
    'path', 'p', 1, "character",
    'contrast-col', 'cc', 1, "character",
    'contrast-str', 'cs', 1, "character",
    'num-threads', 'n', 1, 'integer'
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

if (! is.null(opt[['num-threads']])){
    blas_set_num_threads(opt[['num-threads']])
}

PATH = opt[['path']]
CONTRAST_COL = opt[['contrast-col']]
CONTRAST_STR = opt[['contrast-str']]

trade_dir = file.path(dirname(PATH), "trade")
if (!dir.exists(trade_dir)) {
    dir.create(trade_dir, recursive=TRUE)
}
slogan = gsub(".qs$", "", basename(PATH))
slogan = gsub(".rds$", "", slogan)

res = load_obj(PATH)[["res"]]
fit = res$fit 
coefs = colnames(res$fit)

trade_input_dfs = list()

# 1. prepare TRADE inputs: SE's, logFCs, t stats, un-adjusted P values, all pre-shrinking
# first specify inputs for the contrast (should it be specified), then for all other coefs
cat("\nPreparing TRADE inputs...\n")
if (!is.null(CONTRAST_COL) & !is.null(CONTRAST_STR)){
    cat("Preparing Inputs for contrast column:", CONTRAST_COL, "\n")
    coefs = coefs[!grepl(CONTRAST_COL, coefs)]

    contr = makeContrasts(CONTRAST_STR, levels = res$model)
    fit2 = contrasts.fit(res$fit, contrasts=contr)
    # Calculate standard errors for the contrast coefficients (pre-shrunk)
    se_before_eBayes = as.data.frame(fit2$stdev.unscaled * fit2$sigma)
    se_before_eBayes = se_before_eBayes[rownames(fit2),]
    logFCs = fit2$coef
    logFCs = logFCs[rownames(fit2),]

    ordinary_tstat = fit2$coef / (fit2$stdev.unscaled * fit2$sigma)

    colnames(ordinary_tstat) = c("tstat")
    ordinary_tstat = ordinary_tstat[rownames(fit2),]

    #Two sided t-test p-value
    p_value = as.data.frame(2 * pt(-abs(ordinary_tstat), fit2$df.residual))
    p_value = p_value[rownames(fit2),]

    res_no_ebayes = as.data.frame(cbind(se_before_eBayes, logFCs, ordinary_tstat, p_value))
    colnames(res_no_ebayes) = c("lfcSE", "log2FoldChange", "tstat", "pvalue")
    res_no_ebayes$coef = CONTRAST_COL
    trade_input_dfs[[CONTRAST_COL]] = res_no_ebayes

} else {
    warning("WARNING:\nNo contrast col and/or contrast string specified.\nPerforming TRADE on all coefficients without contrasts.")
}
# now prepare TRADE inputs for all non-contrast coefficients
cat("Preparing Inputs for all other coefficients...\n")
for (coef in coefs){
    se_before_eBayes = setNames(as.data.frame(fit$stdev.unscaled)[[coef]] * fit$sigma, rownames(fit))
    logFCs = setNames(as.data.frame(fit$coef)[[coef]], rownames(fit))
    ordinary_tstat = as.data.frame(fit$coef)[[coef]] / (as.data.frame(fit$stdev.unscaled)[[coef]] * fit$sigma)
    p_value = setNames(2 * pt(-abs(ordinary_tstat), fit$df.residual), rownames(fit)) #Two sided t-test p-value
    res_no_ebayes = data.frame(
        lfcSE = se_before_eBayes, 
        log2FoldChange = logFCs, 
        tstat = ordinary_tstat, 
        pvalue = p_value,
        coef = coef)
    trade_input_dfs[[coef]] = res_no_ebayes
}

# 2. run TRADE for all coefs
trade_outputs = list()
for (coef in names(trade_input_dfs)){
    cat("\nRunning TRADE for coefficient:", coef, "\n")

    input_df = trade_input_dfs[[coef]]
    # OJO: exclude MALAT1 which can have anomalously high p values for low logFCs
    input_df = input_df[!rownames(res_no_ebayes) %in% c("MALAT1"),]
    
    tryCatch({
        trade = TRADE(
            mode="univariate",
            results1=input_df
        )
        trade$coef = coef
        trade_outputs[[coef]] = trade
    }, error = function(e) {
        message(paste("Error in TRADE for coefficient", coef, ":", e$message))})
}

# summarize TRADE outputs
summary_list = list()
cat("\nSummarizing TRADE outputs...\n")
for (coef in names(trade_input_dfs)){
    distsum = trade_outputs[[coef]][["distribution_summary"]]
    distsum$coef = coef
    summary_list[[coef]] = distsum
}

input_df = do.call(rbind, trade_input_dfs)
distsum_df = do.call(rbind, summary_list)

cat("\nSaving TRADE inputs, outputs, and summary...\n")
write.csv(input_df, file.path(trade_dir, paste0(slogan, "__trade_input.csv")), row.names=TRUE)
write.csv(distsum_df, file.path(trade_dir, paste0(slogan, "__trade_distribution_summary.csv")), row.names=TRUE)
qsave(trade_outputs, file.path(trade_dir, paste0(slogan, "__trade_outputs.qs")))
cat("\nTRADE analysis completed and saved in directory:", trade_dir, "\n")
