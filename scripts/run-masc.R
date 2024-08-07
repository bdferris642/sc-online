################## LIBRARIES #################
library(dplyr)
library(getopt)
library(ggplot2)
library(lme4)
library(Matrix)
library(qs)
library(Seurat)
library(RColorBrewer)

source("~/sc-online/masc.R")


################## ARGUMENTS #################
spec <- matrix(c(
    'path', 'p', 1, "character",
    'covariates', 'c', 1, 'character',
    'contrast-col', 'cc', 1, 'character',
    'cluster-col', 'cl', 1, 'character',
    'rand-col', 'rc', 1, 'character',
    'suffix', 's', 1, 'character',
    'continuous', 'C', 1, 'logical'
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
COVARIATES = strsplit(opt[['covariates']], ",")[[1]]
CONTRAST_COL = opt[['contrast-col']]
CLUSTER_COL = opt[['cluster-col']]
RAND_COL = opt[['rand-col']]
CONTINUOUS = if(is.null(opt[['continuous']])){ FALSE }else{ opt[['continuous']] }
# today = format(Sys.Date(), "%y_%m_%d")
suffix = if(!is.null(opt[['suffix']])){
    suffix = paste(CONTRAST_COL, opt[['suffix']], sep="_") 
} else {
    suffix = CONTRAST_COL
}

print(paste("PATH:", PATH))
print(paste("CONTRAST_COL:", CONTRAST_COL))
print(paste("COVARIATES:", COVARIATES))
print(paste("RAND_COL:", RAND_COL))
print(paste("CLUSTER_COL:", CLUSTER_COL))
print(paste("CONTINUOUS:", CONTINUOUS))

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
out_slogan = paste0(slogan, "__masc__", suffix)

print(paste("Output Slogan:", out_slogan))

if (!dir.exists(file.path(base_path, "masc"))) {
    dir.create(file.path(base_path, "masc"))
}

sobj = qread(PATH)

pd = sobj@meta.data[,model_cols]
pd = pd[complete.cases(pd),]

if ("source" %in% model_cols) {
    pd$source = factor(pd$source, levels=c("Calico", "GTEx"))
}
if ("sex" %in% model_cols) {
    pd$sex = factor(pd$sex, levels=c("Male", "Female"))
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
    masc_df = masc_fn_continuous(
        dataset=pd,
        cluster_colname=CLUSTER_COL,
        test_colname=CONTRAST_COL,
        fixed_effects=COVARIATES, 
        random_effects=RAND_COL
    )
    
} else {
    masc_df = .sconline.MASCfn(
        dataset=pd,
        cluster=pd[[CLUSTER_COL]], # cluster annotations
        contrast=CONTRAST_COL, # name of contrast annotations (what you want to run the test for)
        random_effects=RAND_COL, # name of random effects annotations
        fixed_effects=COVARIATES, 
    )

write.csv(masc_df, file.path(base_path, "masc", paste0(out_slogan, ".csv")))

    masc_df$cluster_name = sub("cluster", "", masc_df$cluster)
    masc_df$log2_or = log2(masc_df[[or_colname]])
    masc_df$log2_or_ci_low = log2(masc_df[[ci_low_colname]])
    masc_df$log2_or_ci_high = log2(masc_df[[ci_high_colname]])
    masc_df = masc_df[order(-masc_df$log2_or), ]
    masc_df$cluster_name = factor(masc_df$cluster_name, levels = masc_df$cluster_name[order(masc_df$log2_or)])
    masc_df$rank = rank(masc_df$log2_or)
}

write.csv(masc_df, file.path(base_path, "masc", paste0(out_slogan, ".csv")))

if (CONTINUOUS){

    # Create the forest plot with ticks on error bars, axis lines with ticks, RdBu color map, and opaque white circles on top
    p = ggplot(masc_df, aes(y = cluster_name, x = or)) +
        ggtitle(paste("Change in", str_to_title(gsub("_", " ", CLUSTER_COL)), "Proportion per Additional Unit of", str_to_title(CONTRAST_COL))) +
        xlab("Odds Ratio") +
        ylab(str_to_title(gsub("_", " ", CLUSTER_COL))) +
        geom_vline(xintercept = 1, linetype = "dotted", color = "gray") +  # Add dotted vertical line at x=0
        geom_segment(aes(x = ci_lower, xend = ci_upper, y = cluster_name, yend = cluster_name, color = rank), size = 1) +  # Add horizontal error bars
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




