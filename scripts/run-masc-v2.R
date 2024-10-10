################## LIBRARIES #################
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(lme4)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(qs)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(tidyverse)))

suppressWarnings(suppressMessages(source("~/sc-online/masc.R")))


################## ARGUMENTS #################
spec = matrix(c(
    'path', 'p', 1, "character",
    'covariates', 'c', 1, 'character',
    'contrast-col', 'cc', 1, 'character',
    'cluster-col', 'cl', 1, 'character',
    'rand-col', 'rc', 1, 'character',
    'suffix', 's', 1, 'character',
    'continuous', 'C', 1, 'logical',
    'filter-cluster-n', 'fc', 1, 'numeric',
    'filter-rand-n', 'fr', 1, 'numeric',
    'jk-samples', 'js', 1, 'logical',
    'leave-out', 'lo', 1, 'character'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
COVARIATES = strsplit(opt[['covariates']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
CLUSTER_COL = opt[['cluster-col']]
RAND_COL = opt[['rand-col']]
CONTINUOUS = if(is.null(opt[['continuous']])){ FALSE }else{ opt[['continuous']] }
JK_SAMPLES = if(is.null(opt[['jk-samples']])){ FALSE }else{ opt[['jk']] }
LEAVE_OUT = if(is.null(opt[['leave-out']])){ FALSE }else{ opt[['leave-out']] }

# today = format(Sys.Date(), "%y_%m_%d")
suffix = if(!is.null(opt[['suffix']])){
    suffix = paste(CLUSTER_COL, CONTRAST_COL, opt[['suffix']], sep="__") 
} else {
    suffix = paste(CLUSTER_COL, CONTRAST_COL, sep="__")
}
if (! is.null(LEAVE_OUT)){
    suffix = paste0(suffix, "__leave_out_", LEAVE_OUT)
}

FILTER_CLUSTER_N = if(is.null(opt[['filter-cluster-n']])){ NULL }else{ opt[['filter-cluster-n']] }
FILTER_RAND_N = if(is.null(opt[['filter-rand-n']])){ NULL }else{ opt[['filter-rand-n']] }

print(paste("PATH:", PATH))
print(paste("CONTRAST_COL:", CONTRAST_COL))
print(paste("COVARIATES:", COVARIATES))
print(paste("RAND_COL:", RAND_COL))
print(paste("CLUSTER_COL:", CLUSTER_COL))
print(paste("CONTINUOUS:", CONTINUOUS))
print(paste("FILTER_CLUSTER_N:", FILTER_CLUSTER_N))
print(paste("FILTER_RAND_N:", FILTER_RAND_N))
print(paste("JK_SAMPLES:", JK_SAMPLES))
print(paste("LEAVE_OUT:", LEAVE_OUT))

################## FUNCTIONS #################

str_to_title = function(s){
    s_list = strsplit(s, " ")[[1]]
    s_out = paste(
        lapply(s_list, function(x) {
            paste(toupper(substring(x, 1, 1)), substring(x, 2), sep="")
        }), collapse=" ")
    return(s_out)
}

MASC = function(dataset, cluster_col, contrast_col, random_effects = NULL, fixed_effects = NULL, verbose = TRUE) {

    # Convert cluster assignments to string
    cluster = dataset[[cluster_col]]
    cluster = as.character(cluster)
    contrast = dataset[[contrast_col]]

    # Check inputs
    if (is.factor(contrast) == FALSE) {
        stop("Specified contrast term is not coded as a factor in dataset!")
    }

    if (is.null(fixed_effects)){
        stop("No fixed effects specified!")
    }

    # For now, do not allow models without mixed effects terms
    if (is.null(random_effects)){
        stop("No random effects specified!")
    }

    # Prepend design matrix generated from cluster assignments
    designmat = model.matrix(~ cluster + 0, data.frame(cluster = cluster))
    dataset = cbind(designmat, dataset)

    # Create output list to hold results
    res = vector(mode = "list", length = length(unique(cluster)))
    names(res) = attributes(designmat)$dimnames[[2]]
    
    # Create model formulas

    model_rhs = paste0(c(
        paste0(fixed_effects, collapse = " + "),
        paste0("(1|", random_effects, ")", collapse = " + ")),
    collapse = " + ")

    if (verbose == TRUE) {
        message(paste("Using null model:", "cluster ~", model_rhs))
    }

    # Run nested mixed-effects models for each cluster
    for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
        test_cluster = attributes(designmat)$dimnames[[2]][i]
        if (verbose == TRUE) {
            message(paste("Creating logistic mixed models for", test_cluster))
        }
        null_fm = as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                        model_rhs), collapse = ""))
        full_fm = as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                        model_rhs), collapse = ""))
        # Run null and full mixed-effects models
        null_model = lme4::glmer(formula = null_fm, data = dataset,
                                    family = binomial, nAGQ = 1, verbose = 0,
                                    control = glmerControl(optimizer = "bobyqa"))
        full_model = lme4::glmer(formula = full_fm, data = dataset,
                                    family = binomial, nAGQ = 1, verbose = 0,
                                    control = glmerControl(optimizer = "bobyqa"))
        model_lrt = anova(null_model, full_model)
        # calculate confidence intervals for contrast term beta
        contrast_lvl2 = paste0(contrast, levels(contrast)[2])
        contrast_ci = confint.merMod(full_model, method = "Wald",
                                        parm = contrast_lvl2)
        # Save model objects to list
        cluster_models[[i]]$null_model = null_model
        cluster_models[[i]]$full_model = full_model
        cluster_models[[i]]$model_lrt = model_lrt
        cluster_models[[i]]$confint = contrast_ci
    }
    # Organize results into output dataframe
    output = data.frame(cluster = attributes(designmat)$dimnames[[2]],
                        size = colSums(designmat))
    output$model.pvalue = sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
    output[[paste(contrast_lvl2, "OR", sep = ".")]] = sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
    output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] = sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
    output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] = sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))

    return(output)

}


################## MAIN #################

model_cols = c(CONTRAST_COL, COVARIATES, CLUSTER_COL, RAND_COL)
base_path_list = strsplit(PATH, "/")[[1]]
base_path = paste(base_path_list[1:(length(base_path_list)-1)], collapse="/")
basename = base_path_list[[length(base_path_list)]]
slogan = gsub(".qs", "", basename)
out_slogan = paste0(slogan, "__masc", suffix)

print(paste("Output Slogan:", out_slogan))

if (!dir.exists(file.path(base_path, "masc"))) {
    dir.create(file.path(base_path, "masc"))
}

sobj = qread(PATH)
print("Seurat Object Dimensions")
print(dim(sobj))
pd = sobj@meta.data[,model_cols]
print("Seurat Object Meta Data Columns")
print(colnames(pd))

if (!is.null(LEAVE_OUT)) {
    cluster_vec = pd[[CLUSTER_COL]]
    print(head(cluster_vec))
    print(length(cluster_vec))
    print(sum(cluster_vec == LEAVE_OUT))
    print(dim(pd))
    print(head(pd))
    pd = pd[cluster_vec != LEAVE_OUT, ]
}

if (!is.null(FILTER_CLUSTER_N)) {
    pd = pd %>% group_by(pd[[CLUSTER_COL]]) %>% filter(n() >= FILTER_CLUSTER_N) 
}

if (!is.null(FILTER_RAND_N) & !is.null(RAND_COL)) {
    pd = pd %>% group_by(pd[[RAND_COL]]) %>% filter(n() >= FILTER_RAND_N)
}

pd = pd[complete.cases(pd),]

if ("source" %in% model_cols) {
    pd$source = factor(pd$source, levels=c("Calico", "GTEx"))
}
if ("sex" %in% model_cols) {
    if ("F" %in% unique(pd$sex)) {
        pd$sex = factor(pd$sex, levels=c("M", "F"))
    }
    if ("Female" %in% unique(pd$sex)) {
        pd$sex = factor(pd$sex, levels=c("Male", "Female"))
    }
}


# replace dashes with underscores in actual column data
# ojo: update here if some rude person puts +, &, |, etc. in the column names
for (col in model_cols) {
    pd[[col]] = gsub("-", "_", pd[[col]])
}

# note the strange beahvior of the case_control column
# by ordering asc, the column will be called case_controlPD or case_controlXDP 
# but when we grab the case_name a second, we need to reverse the order for some reason.
if ("case_control" %in% model_cols) {
    pd$case_control = factor(pd$case_control, levels=sort(unique(as.character(pd$case_control)), decreasing=F))
}
case_name = sort(unique(as.character(pd[[CONTRAST_COL]])), decreasing=T)[[1]]
ctr_name = sort(unique(as.character(pd[[CONTRAST_COL]])), decreasing=T)[[2]]

print(paste("Case Name:", case_name))
print(paste("Control Name:", ctr_name))

or_colname = paste0(CONTRAST_COL, case_name, ".OR")
ci_low_colname = paste0(CONTRAST_COL, case_name, ".OR.95pct.ci.lower")
ci_high_colname = paste0(CONTRAST_COL, case_name, ".OR.95pct.ci.upper")

print(paste("OR Column Name:", or_colname))
print(paste("CI High Column Name:", ci_high_colname))
print(paste("CI Low Column Name:", ci_low_colname))

if (CONTINUOUS) {
    print("foo")
} else {
    masc_df = MASC(
        pd, 
        cluster_col=CLUSTER_COL, 
        contrast_col=CONTRAST_COL, 
        random_effects=RAND_COL, 
        fixed_effects=COVARIATES,
        verbose=TRUE)

    covariate_str = paste(COVARIATES, collapse=", ")
    masc_df$covariates = covariate_str
    masc_df$random_effect = RAND_COL
    masc_df$contrast = CONTRAST_COL
    masc_df$cluster_by = CLUSTER_COL
    masc_df$filter_cluster_n = FILTER_CLUSTER_N
    masc_df$filter_rand_n = FILTER_RAND_N
    masc_df$cluster_name = sub("cluster", "", masc_df$cluster)
    write.csv(masc_df, file.path(base_path, "masc", paste0(out_slogan, ".csv")))
    
    masc_df$log2_or = log2(masc_df[[or_colname]])
    masc_df$log2_or_ci_low = log2(masc_df[[ci_low_colname]])
    masc_df$log2_or_ci_high = log2(masc_df[[ci_high_colname]])
    masc_df = masc_df[order(-masc_df$log2_or), ]
    masc_df$cluster_name = factor(masc_df$cluster_name, levels = masc_df$cluster_name[order(masc_df$log2_or)])
    masc_df$rank = rank(masc_df$log2_or)
    write.csv(masc_df, file.path(base_path, "masc", paste0(out_slogan, ".csv")))
} 



if (CONTINUOUS){
    print("bar")
} else {
    # Create the forest plot with ticks on error bars, axis lines with ticks, RdBu color map, and opaque white circles on top
    p = ggplot(masc_df, aes(y = cluster_name, x = log2_or)) +
        ggtitle(paste0(str_to_title(gsub("_", " ", CLUSTER_COL)), " Enrichment in ", toupper(case_name), " vs. ", str_to_title(tolower(ctr_name)))) +  # Title
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +  # Add dotted vertical line at x=0
        geom_segment(aes(x = log2_or_ci_low, xend = log2_or_ci_high, y = cluster_name, yend = cluster_name, color = rank), size = 1) +  # Add horizontal error bars
        geom_point(size = 3, aes(color = log2_or), shape = 1) +  # Add points for effect sizes
        geom_point(size = 3, shape = 21, fill = "white") +  # Add opaque white circle on top of the error bar line
        scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdBu")) +  # Use RdBu color map
        theme_minimal() +  # Minimal theme
        labs(x = "log2(OR)", y = str_to_title(gsub("_", " ", CLUSTER_COL))) +  # Axis labels
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = "none",
            plot.title = element_text(size=16),
            axis.line = element_line(color = "black"),  # Add axis lines
            axis.ticks = element_line(color = "black"),  # Add axis ticks
            axis.text = element_text(size = 14),  # Increase tick label font size
            axis.title = element_text(size = 15)  # Increase axis label font size
        )

    ggsave(file.path(base_path, "masc", paste0(out_slogan, ".png")), plot = p, width = 11, height = 7, bg="white", dpi=400)

}