library(Seurat)
library(devtools)
library(Matrix)
library(qs)
library(TRADE)
library(ashr)
library(limma)
library(ggplot2)
library(getopt)
library(Matrix)

# `path` is the absolute path to a psuedocell directory (should be in `de` directory)
# this pseudocell directory should contain .qs files of DE results

# TODO: make generalizable to whatever contrast is
# TODO: make sexMale vs sexFemale generalizable

spec <- matrix(c(
    'path', 'pp', 1, "character"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

print(spec)
print(opt)

path = opt[['path']]
base_path = dirname(path)

pseudocell_basename = basename(path)
pseudocell_slogan = gsub("\\.qs$", "", pseudocell_basename)

print("DE FILES:")
list.files(path, pattern = ".*.qs$", full.names = FALSE)
de_basenames = list.files(path, pattern = ".*.qs$", full.names = FALSE)

trade_dir = file.path(path, "trade")
print("")
print(paste("WRITING OUTPUT TO:", trade_dir))
if (!dir.exists(trade_dir)) {
    print(paste("Creating directory", trade_dir))
    dir.create(trade_dir, recursive=TRUE)
    dir.create(file.path(trade_dir, "case_control"), recursive=TRUE)
    dir.create(file.path(trade_dir, "age"), recursive=TRUE)
    dir.create(file.path(trade_dir, "sex"), recursive=TRUE)
}

trade_dfs_pd = list()
trade_outputs_pd = list()
summary_list_pd = list()

trade_dfs_age = list()
trade_outputs_age = list()
summary_list_age = list()

trade_dfs_sex = list()
trade_outputs_sex = list()
summary_list_sex = list()

print("")
print("LOADING TRADE INPUTS:")
for (de_basename in de_basenames) {

    slogan = sub(".qs$", "", de_basename)
    if (length(strsplit(slogan, "__")[[1]]) <= 3){
        ct = strsplit(slogan, "__")[[1]][[2]]
    } else {
        ct = paste0(strsplit(slogan, "__")[[1]][[2]], '__', strsplit(slogan, "__")[[1]][[3]])
    }

    print(paste("Loading", de_basename, "and storing as", ct))
    de_path = file.path(path, de_basename)
    de_obj = qread(de_path)
    res = de_obj[['res']]
    
    print("Case/Control")
    # PD vs Control
    contr = makeContrasts(case_controlpd - case_controlctr, levels = res$model)
    fit2 = contrasts.fit(res$fit, contrasts=contr)

    # PD
    # Calculate standard errors for the contrast coefficients (pre-shrunk)
    se_before_eBayes_pd = as.data.frame(fit2$stdev.unscaled * fit2$sigma)
    se_before_eBayes_pd = se_before_eBayes_pd[rownames(fit2),]

    logFCs_pd = fit2$coef
    logFCs_pd = logFCs_pd[rownames(fit2),]
    
    ordinary_tstat_pd = fit2$coef / (fit2$stdev.unscaled * fit2$sigma)
    colnames(ordinary_tstat_pd) = c("tstat")
    ordinary_tstat_pd = ordinary_tstat_pd[rownames(fit2),]
    
    #Two sided t-test p-value
    p_value_pd = as.data.frame(2 * pt(-abs(ordinary_tstat_pd), fit2$df.residual))
    p_value_pd = p_value_pd[rownames(fit2),]
    
    res_no_ebayes_pd = as.data.frame(cbind(se_before_eBayes_pd, logFCs_pd, ordinary_tstat_pd, p_value_pd))
    colnames(res_no_ebayes_pd) = c("lfcSE", "log2FoldChange", "tstat", "pvalue")
    trade_dfs_pd[[ct]] = res_no_ebayes_pd

    print("Age")
    fit = res$fit 
    # Age
    se_before_eBayes_age = setNames(as.data.frame(fit$stdev.unscaled)$age * fit$sigma, rownames(fit))
    logFCs_age = setNames(as.data.frame(fit$coef)$age, rownames(fit))
    ordinary_tstat_age = as.data.frame(fit$coef)$age / (as.data.frame(fit$stdev.unscaled)$age * fit$sigma)
    p_value_age = setNames(2 * pt(-abs(ordinary_tstat_age), fit$df.residual), rownames(fit)) #Two sided t-test p-value
    res_no_ebayes_age = data.frame(
        lfcSE = se_before_eBayes_age, 
        log2FoldChange = logFCs_age, 
        tstat = ordinary_tstat_age, 
        pvalue = p_value_age)
    trade_dfs_age[[ct]] = res_no_ebayes_age

    # Sex
    print("Sex")
    se_before_eBayes_sex = setNames(as.data.frame(fit$stdev.unscaled)$sexMale * fit$sigma, rownames(fit))
    logFCs_sex = setNames(as.data.frame(fit$coef)$sexMale, rownames(fit))
    ordinary_tstat_sex = as.data.frame(fit$coef)$sexMale / (as.data.frame(fit$stdev.unscaled)$sexMale * fit$sigma)
    p_value_sex = setNames(2 * pt(-abs(ordinary_tstat_sex), fit$df.residual), rownames(fit)) #Two sided t-test p-value
    res_no_ebayes_sex = data.frame(
        lfcSE = se_before_eBayes_sex, 
        log2FoldChange = logFCs_sex, 
        tstat = ordinary_tstat_sex, 
        pvalue = p_value_sex)
    trade_dfs_sex[[ct]] = res_no_ebayes_sex

}
print("")
print("")
print("RUNNING TRADE: Case/Control")
for (cc in names(trade_dfs_pd)) {
    print(paste("Running TRADE for", cc))
    res = trade_dfs_pd[[cc]]
    trade = TRADE(
        mode="univariate",
        results1=res
    )
    trade$cell_class = cc
    trade_outputs_pd[[cc]] = trade
}
print("")
print("RUNNING TRADE: Age")
for (cc in names(trade_dfs_age)) {
    print(paste("Running TRADE for", cc))
    res = trade_dfs_age[[cc]]
    trade = TRADE(
        mode="univariate",
        results1=res
    )
    trade$cell_class = cc
    trade_outputs_age[[cc]] = trade
}
print("")
print("RUNNING TRADE: Sex")
for (cc in names(trade_dfs_sex)) {
    print(paste("Running TRADE for", cc))
    res = trade_dfs_sex[[cc]]
    trade = TRADE(
        mode="univariate",
        results1=res
    )
    trade$cell_class = cc
    trade_outputs_sex[[cc]] = trade
}

print("")
print("")
print("SUMMARIZING TRADE OUTPUTS")
print("Case/Control")
for (cc in names(trade_outputs_pd)) {
    res = trade_outputs_pd[[cc]][["distribution_summary"]]
    res$cell_class = cc
    summary_list_pd[[cc]] = res
}
summary_df_pd = do.call(rbind, summary_list_pd)

print("Age")
for (cc in names(trade_outputs_age)) {
    res = trade_outputs_age[[cc]][["distribution_summary"]]
    res$cell_class = cc
    summary_list_age[[cc]] = res
}
summary_df_age = do.call(rbind, summary_list_age)

print("Sex")
for (cc in names(trade_outputs_sex)) {
    res = trade_outputs_sex[[cc]][["distribution_summary"]]
    res$cell_class = cc
    summary_list_sex[[cc]] = res
}
summary_df_sex = do.call(rbind, summary_list_sex)

print("")
print("")
print("WRITING DATA")
print("Case/Control")
write.csv(summary_df_pd, file.path(trade_dir, "case_control/trade_summary_df__case_control.csv"))
qsave(trade_dfs_pd, file.path(trade_dir, "case_control/trade_inputs__case_control.qs"))
qsave(trade_outputs_pd, file.path(trade_dir, "case_control/trade_outputs__case_control.qs"))

print("Age")
write.csv(summary_df_age, file.path(trade_dir, "age/trade_summary_df__age.csv"))
qsave(trade_dfs_age, file.path(trade_dir, "age/trade_inputs__age.qs"))
qsave(trade_outputs_age, file.path(trade_dir, "age/trade_outputs__age.qs"))

print("Sex")
write.csv(summary_df_sex, file.path(trade_dir, "sex/trade_summary_df__sex.csv"))
qsave(trade_dfs_sex, file.path(trade_dir, "sex/trade_inputs__sex.qs"))
qsave(trade_outputs_sex, file.path(trade_dir, "sex/trade_outputs__sex.qs"))