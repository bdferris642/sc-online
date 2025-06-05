# A script that, given a Seurat Object with metadata,
# Trains a mixed effect model to infer statistically significant differences in cluster proportions
# Saving output as a csv as well as figures

# path: path to a seurat .qs
# covariates: comma-separated list of covariates to include in the model
# contrast-col: the column in the metadata that contains the contrast to test
# cluster-col: the column in the metadata that contains the cluster annotations to iterate over (e.g. cell_class)
# rand-col: the column in the metadata that contains the random effects annotations (e.g. donor_id)
# suffix (optional): a suffix to add to the output file name
# continuous (optional): whether the contrast is continuous or categorical (default: FALSE)
# filter-cluster-n (optional): minimum number of cells in a cluster to include it in the analysis
    # (cell classes with very few cells can have extremely high error estimates)
# filter-rand-n (optional): minimum number of cells in a random effect split to include in the analysis 
    # (don't want to be biased by proportions coming from donors with very few cells)
# leave-out (optional): a comma-separated list of clusters to leave out of the analysis
# num-threads (optional): number of threads to use for the analysis

################## LIBRARIES #################
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(lme4)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(qs)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(RhpcBLASctl)))

suppressWarnings(suppressMessages(source("~/sc-online/masc.R")))


################## ARGUMENTS #################
spec <- matrix(c(
    'path', 'p', 1, "character",
    'covariates', 'c', 1, 'character',
    'contrast-col', 'cc', 1, 'character',
    'cluster-col', 'cl', 1, 'character',
    'rand-col', 'rc', 1, 'character',
    'suffix', 's', 1, 'character',
    'continuous', 'C', 1, 'logical',
    'filter-cluster-n', 'fc', 1, 'numeric',
    'filter-rand-n', 'fr', 1, 'numeric',
    'leave-out', 'lo', 1, 'character',
    'num-threads', 'n', 1, 'integer'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
COVARIATES = strsplit(opt[['covariates']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
CLUSTER_COL = opt[['cluster-col']]
RAND_COL = opt[['rand-col']]
CONTINUOUS = if(is.null(opt[['continuous']])){ FALSE }else{ opt[['continuous']] }
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

if (! is.null(opt[['num-threads']])){
    blas_set_num_threads(opt[['num-threads']])
}

# TODO: parametrize
JK_SAMPLES = FALSE

message(paste("PATH:", PATH))
message(paste("CONTRAST_COL:", CONTRAST_COL))
message(paste("COVARIATES:", COVARIATES))
message(paste("RAND_COL:", RAND_COL))
message(paste("CLUSTER_COL:", CLUSTER_COL))
message(paste("CONTINUOUS:", CONTINUOUS))
message(paste("FILTER_CLUSTER_N:", FILTER_CLUSTER_N))
message(paste("FILTER_RAND_N:", FILTER_RAND_N))
message(paste("JK_SAMPLES:", JK_SAMPLES))
message(paste("LEAVE_OUT:", LEAVE_OUT))
message(paste("SUFFIX:", suffix))
message(paste("NUM_THREADS:", opt[['num-threads']]))

################## FUNCTIONS #################

str_to_title = function(s){
    s_list = strsplit(s, " ")[[1]]
    s_out = paste(
        lapply(s_list, function(x) {
            paste(toupper(substring(x, 1, 1)), substring(x, 2), sep="")
        }), collapse=" ")
    return(s_out)
}


################## MAIN #################

model_cols = c(CONTRAST_COL, COVARIATES, CLUSTER_COL, RAND_COL)
base_path_list = strsplit(PATH, "/")[[1]]
base_path = paste(base_path_list[1:(length(base_path_list)-1)], collapse="/")
basename = base_path_list[[length(base_path_list)]]
slogan = gsub(".qs", "", basename)
out_slogan = paste0(slogan, "__masc_", suffix)

message(paste("Output Slogan:", out_slogan))

if (!dir.exists(file.path(base_path, "masc"))) {
    dir.create(file.path(base_path, "masc"))
}

sobj = qread(PATH)
message("Seurat Object Dimensions")
message(dim(sobj))
pd = sobj@meta.data[,model_cols]
pd[[CLUSTER_COL]] = factor(
    pd[[CLUSTER_COL]], 
    levels= sort(unique(as.character(pd[[CLUSTER_COL]]))))
message("Seurat Object Meta Data Columns")
message(colnames(pd))


if (!is.null(LEAVE_OUT)) {
    cluster_vec = pd[[CLUSTER_COL]]
    message("Cluster Vector Length and Head")
    message(length(cluster_vec))
    message(head(cluster_vec))
    message("Leave Out Vector Length and Head")
    message(sum(cluster_vec == LEAVE_OUT))
    message("Phenotype Data Dimensions")
    message(dim(pd))
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
for (col in model_cols) {
    pd[[col]] = gsub("-", "_", pd[[col]])
    pd[[col]] = gsub(" ", "_", pd[[col]])
    pd[[col]] = gsub("\\.", "", pd[[col]])
    pd[[col]] = gsub("\\,", "", pd[[col]])
    pd[[col]] = gsub("\\?", "", pd[[col]])
    pd[[col]] = gsub("\\+", "", pd[[col]])
    pd[[col]] = gsub("\\&", "", pd[[col]])
    pd[[col]] = gsub("\\!", "", pd[[col]])
    pd[[col]] = gsub("\\|", "", pd[[col]])
    pd[[col]] = gsub("\\/", "", pd[[col]])
}

# note the strange beahvior of the case_control column
# by ordering asc, the column will be called case_controlPD or case_controlXDP 
# but when we grab the case_name a second, we need to reverse the order for some reason.
if ("case_control" %in% model_cols) {
    pd$case_control = factor(pd$case_control, levels=sort(unique(as.character(pd$case_control)), decreasing=F))
}
case_name = sort(unique(as.character(pd[[CONTRAST_COL]])), decreasing=T)[[1]]
ctr_name = sort(unique(as.character(pd[[CONTRAST_COL]])), decreasing=T)[[2]]

message(paste("Case Name:", case_name))
message(paste("Control Name:", ctr_name))

or_colname = paste0(CONTRAST_COL, case_name, ".OR")
ci_low_colname = paste0(CONTRAST_COL, case_name, ".OR.95pct.ci.lower")
ci_high_colname = paste0(CONTRAST_COL, case_name, ".OR.95pct.ci.upper")

message(paste("OR Column Name:", or_colname))
message(paste("CI High Column Name:", ci_high_colname))
message(paste("CI Low Column Name:", ci_low_colname))

if (CONTINUOUS) {
    masc_df = masc_fn_continuous(
        dataset=pd,
        cluster_colname=CLUSTER_COL,
        test_colname=CONTRAST_COL,
        fixed_effects=COVARIATES, 
        random_effects=RAND_COL
    )
    
} else {
    message("Running MASC for categorical data")
    message("Dimension of input data:")
    message(dim(pd))
    message("Number of Complete Cases")
    message(sum(complete.cases(pd)))

    # print summary of all columns, including class
    message("Summary of input data:")
    print(summary(pd))
    print(head(pd))

    pd[[RAND_COL]] = as.factor(pd[[RAND_COL]])
    
    masc_df = .sconline.MASCfn(
        dataset=pd,
        cluster=pd[[CLUSTER_COL]], # cluster annotations
        contrast=CONTRAST_COL, # name of contrast annotations (what you want to run the test for)
        random_effects=RAND_COL, # name of random effects annotations
        fixed_effects=COVARIATES,
        jackknife = JK_SAMPLES
    )

    masc_df$cluster_name = sub("cluster", "", masc_df$cluster)
    masc_df$log2_or = log2(masc_df[[or_colname]])
    masc_df$log2_or_ci_low = log2(masc_df[[ci_low_colname]])
    masc_df$log2_or_ci_high = log2(masc_df[[ci_high_colname]])
    masc_df = masc_df[order(-masc_df$log2_or), ]
    masc_df$cluster_name = factor(masc_df$cluster_name, levels = masc_df$cluster_name[order(masc_df$log2_or)])
    masc_df$rank = rank(masc_df$log2_or)
}

covariate_str = paste(COVARIATES, collapse=", ")
masc_df$covariates = covariate_str
masc_df$random_effect = RAND_COL
masc_df$contrast = CONTRAST_COL
masc_df$cluster_by = CLUSTER_COL
masc_df$filter_cluster_n = FILTER_CLUSTER_N
masc_df$filter_rand_n = FILTER_RAND_N

write.csv(masc_df, file.path(base_path, "masc", paste0(out_slogan, ".csv")))

if (CONTINUOUS){

    # Create the forest plot with ticks on error bars, axis lines with ticks, RdBu color map, and opaque white circles on top
    p = ggplot(masc_df, aes(y = cluster_name, x = or)) +
        ggtitle(paste(
            "Change in", str_to_title(gsub("_", " ", CLUSTER_COL)), 
            "Proportion per Additional Unit of", str_to_title(CONTRAST_COL))) +
        xlab("Odds Ratio") +
        ylab(str_to_title(gsub("_", " ", CLUSTER_COL))) +
        geom_vline(xintercept = 1, linetype = "dotted", color = "gray") +  # Add dotted vertical line at x=0
        geom_segment(aes(
            x = ci_lower, xend = ci_upper, 
            y = cluster_name, yend = cluster_name, color = rank), size = 1) +  # Add horizontal error bars
        geom_point(size = 3, aes(color = or), shape = 1) +  # Add points for effect sizes
        geom_point(size = 3, shape = 21, fill = "white") +  # Add opaque white circle on top of the error bar line
        scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "RdBu")) +  # Use RdBu color map
        theme_minimal()  + # Apply a minimal theme
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

} else {

    # Create the forest plot with ticks on error bars, axis lines with ticks, RdBu color map, and opaque white circles on top
    p = ggplot(masc_df, aes(y = cluster_name, x = log2_or)) +
        ggtitle(paste0(str_to_title(
            gsub("_", " ", CLUSTER_COL)), 
            " Enrichment in ", toupper(case_name), " vs. ", 
            str_to_title(tolower(ctr_name)))) +  # Title
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray") +  # Add dotted vertical line at x=0
        geom_segment(aes(
            x = log2_or_ci_low, xend = log2_or_ci_high, 
            y = cluster_name, yend = cluster_name, color = rank), size = 1) +  # Add horizontal error bars
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

    ggsave(file.path(base_path, "masc", paste0(out_slogan, ".png")), 
        plot = p, width = 11, height = 7, bg="white", dpi=400)

}

