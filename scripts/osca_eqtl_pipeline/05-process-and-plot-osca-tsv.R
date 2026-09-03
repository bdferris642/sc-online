# 05-process-and-plot-osca-tsv.R — Apply two-stage FDR to OSCA output; save RDS; generate plots.
#
# Arguments:
#   --path / -p   (required) Path to eqtl_{cell_class}.tsv produced by step 3.
#                            File must be named eqtl_{cell_class}.tsv; cell_class is extracted
#                            from the filename as everything after the leading "eqtl_" prefix.
#   --out-dir / -o (optional) Output directory. Defaults to dirname(--path).
#
# Two-stage FDR procedure:
#   1. Within each gene (Probe): Bonferroni across all cis-SNPs → p_bonf
#   2. Minimum p_bonf per gene → BH across genes → padj_gene
#   3. eGenes: padj_gene < FDR_THRESH
#   4. Within eGenes: BH across all SNPs → padj_snp
#   5. Significant eSNPs: padj_gene < FDR_THRESH AND padj_snp < FDR_THRESH
#
# Outputs (in out-dir):
#   eqtl_{cell_class}.rds       — full data frame with FDR columns
#   eqtl_{cell_class}_sig.rds   — significant rows only
#   plots/eqtl_{cell_class}_manhattan.png
#   plots/eqtl_{cell_class}_min_p_gene_hist.png
#   plots/eqtl_{cell_class}_padj_gene_hist.png
#   plots/eqtl_{cell_class}_pval_hist.png
#   plots/eqtl_{cell_class}_padj_snp_hist.png

print("**************** LOADING LIBRARIES ****************")
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(ggrepel)))


print("**************** PARSING ARGUMENTS ****************")
spec = matrix(c(
    'path', 'p', 1, "character",
    'out-dir', 'o', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
if (is.null(opt[["out-dir"]])) {
    OUT_DIR = dirname(PATH)
} else {
    OUT_DIR = opt[["out-dir"]]
}

base = dirname(PATH)
# Extract cell_class from filename: input must be named eqtl_{cell_class}.tsv.
# sub("^eqtl_", ...) handles underscores within cell class names (e.g. "da_neuron").
slogan = tools::file_path_sans_ext(basename(PATH))   # "eqtl_da_neuron"
cell_class = sub("^eqtl_", "", slogan)                # "da_neuron"
print(slogan)
print(cell_class)
plot_dir = file.path(OUT_DIR, "plots")
dir.create(OUT_DIR, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

# HARD_CAP: maximum -log10(p) displayed in plots; values above this are capped to
#           prevent extreme outliers from compressing the scale.
HARD_CAP = 40
# FDR_THRESH: BH adjusted p-value threshold for declaring eGenes and eSNPs significant.
FDR_THRESH = 0.05

cat(paste0("\n READING DataFrame from", PATH, "\n"))
df = read.table(PATH, header = T, sep = "\t")
df$BP = as.numeric(df$BP)
df$negative_log10_p = -log10(df$p)
df$negative_log10_p[df$negative_log10_p > HARD_CAP] = HARD_CAP

# Two-stage FDR: within-gene Bonferroni → BH across genes (padj_gene) → BH within eGenes (padj_snp).
# The commented-out alternative approaches were explored during development; active implementation below.

cat(paste0("\n PERFORMING TWO-STAGE FDR\n"))
cat(paste0("\n Step 1: Find lead SNPs per gene\n"))

# Step 1: Bonferroni correction within each gene
df_probe = df %>%
    group_by(Probe) %>%
    mutate(p_bonf = p * n()) %>%
    mutate(p_bonf = ifelse(p_bonf > 1, 1, p_bonf)) %>%
    ungroup()

# Step 2: Get the minimum p_bonf per gene
df_min_p_probe = df_probe %>%
    group_by(Probe) %>%
    summarise(min_p_gene = min(p_bonf)) %>%
    ungroup()

# Step 3: BH correct these min_p_gene values across genes
df_min_p_probe = df_min_p_probe %>%
    mutate(padj_gene = p.adjust(min_p_gene, method = "BH"))

# Step 4: Identify eGenes (FDR < threshold)
eGenes = df_min_p_probe %>%
    filter(padj_gene < FDR_THRESH) %>%
    pull(Probe)

# Step 5: Apply BH FDR to SNPs within eGenes
df_eGenes = df %>%
    filter(Probe %in% eGenes) %>%
    mutate(padj_snp = p.adjust(p, method = "BH"))

# Step 6: Merge the min_p_gene and padj_gene back to the df
df = df %>%
    left_join(df_min_p_probe[, c("Probe", "min_p_gene", "padj_gene")], by = "Probe") %>%
    left_join(df_eGenes[, c("SNP", "Probe", "padj_snp")], by = c("SNP", "Probe")) %>%
    mutate(
        padj_snp = ifelse(is.na(padj_snp), 1, padj_snp),
        is_significant_gene = padj_gene < FDR_THRESH,
        is_significant_snp = padj_snp < FDR_THRESH
    )

df$negative_log10_padj_snp = -log10(df$padj_snp)
df$negative_log10_padj_snp[df$negative_log10_padj_snp > HARD_CAP] = HARD_CAP

df$negative_log10_padj_gene = -log10(df$padj_gene)
df$negative_log10_padj_gene[df$negative_log10_padj_gene > HARD_CAP] = HARD_CAP


print(dim(df))
print(colnames(df))
cat(paste0("\nSAVING DATA\n"))
saveRDS(df, file.path(OUT_DIR, paste0(slogan, ".rds")))

df_sig = df[df$is_significant_gene & df$is_significant_snp,]
saveRDS(df_sig, file.path(OUT_DIR, paste0(slogan, "_sig.rds")))

print(colnames(df))

cat("\n PLOTTING DATA\n")
cat("\n GENE LEVEL Min(P) AND Padj HISTOGRAMS\n")

p_gene_hist = (
    ggplot(df_min_p_probe, aes(x = min_p_gene))
    + geom_histogram(binwidth=0.02, fill = "blue", color = "black", alpha = 0.5)
    + coord_cartesian(xlim = c(0, 1))
    + ggtitle(paste0("eQTL Gene Minimum P-value Histogram:\n", cell_class))
    + theme(
            plot.title = element_text(size = 22), # title font size
            axis.text = element_text(size = 20),  # Increase tick label font size
            axis.title = element_text(size = 20),  # Increase axis label font size
    )
)
print(p_gene_hist)
ggsave(
    file.path(plot_dir, paste0(slogan, "_min_p_gene_hist.png")),
    plot=p_gene_hist, width=8, height=6, dpi=600)

padj_gene_hist = (
    ggplot(df_min_p_probe, aes(x = padj_gene))
    + geom_histogram(binwidth=0.02, fill = "blue", color = "black", alpha = 0.5)
    + coord_cartesian(xlim = c(0, 1))
    + ggtitle(paste0("eQTL Gene Adj. P-value Histogram:\n", cell_class))
    + theme(
            plot.title = element_text(size = 22), # title font size
            axis.text = element_text(size = 20),  # Increase tick label font size
            axis.title = element_text(size = 20),  # Increase axis label font size
    )
)
print(padj_gene_hist)
ggsave(
    file.path(plot_dir, paste0(slogan, "_padj_gene_hist.png")),
    plot=padj_gene_hist, width=8, height=6, dpi=600)


cat(paste0("\n RAW SNP PVAL HISTOGRAM\n"))
pval_hist = (
    ggplot(df, aes(x = p))
    + geom_histogram(binwidth=0.02, fill = "blue", color = "black", alpha = 0.5)
    + coord_cartesian(xlim = c(0, 1))
    + ggtitle(paste0("eQTL P-value Histogram:\n", cell_class))
    + theme(
            plot.title = element_text(size = 22), # title font size
            axis.text = element_text(size = 20),  # Increase tick label font size
            axis.title = element_text(size = 20),  # Increase axis label font size
    )
)

cat(paste0("\n SNP ADJ PVAL HISTOGRAM\n"))
print(pval_hist)
ggsave(
    file.path(plot_dir, paste0(slogan, "_pval_hist.png")),
    plot=pval_hist, width=8, height=6, dpi=600)

padj_hist = (
    ggplot(df, aes(x = padj_snp))
    + geom_histogram(binwidth=0.02, fill = "blue", color = "black", alpha = 0.5)
    + coord_cartesian(xlim = c(0, 1))
    + ggtitle(paste0("eQTL SNP Adj. P-value Histogram:\n", cell_class))
    + theme(
            plot.title = element_text(size = 22), # title font size
            axis.text = element_text(size = 20),  # Increase tick label font size
            axis.title = element_text(size = 20),  # Increase axis label font size
    )
)
print(padj_hist)
ggsave(
    file.path(plot_dir, paste0(slogan, "_padj_snp_hist.png")),
    plot=padj_hist, width=8, height=6, dpi=600)


options(repr.plot.width = 20, repr.plot.height = 12)

cat(paste0("\n PREPARING DATA FOR MANHATTAN PLOT\n"))
# Prepare df by calculating chromosome midpoints and cumulative base pair positions.
# Only chromosomes 1-22 are included; sex chromosomes and MT are dropped with a warning.
chr_order <- as.character(1:22)

unsupported_chrs = setdiff(unique(df$Chr), chr_order)
if (length(unsupported_chrs) > 0)
    warning(paste("Dropping", length(unsupported_chrs), "chromosomes not in 1:22:",
                  paste(unsupported_chrs, collapse=", ")))

df = df %>%
    mutate(Chr = factor(Chr, levels=chr_order)) %>%
    arrange(Chr, BP)

df_sig = df[df$is_significant_gene & df$is_significant_snp,]
df_insig = df[!df$is_significant_gene | !df$is_significant_snp,]

# Subsample 1/50th of insignificant points for plotting performance.
# set.seed ensures reproducibility of the random subsample.
set.seed(42)
df_insig_subset = if (nrow(df_insig) > 0) {
    df_insig[sample(nrow(df_insig), max(1L, floor(nrow(df_insig) / 50L))), ]
} else {
    df_insig
}
df_plot = rbind(df_sig, df_insig_subset)
df_plot = df_plot %>% arrange(Chr, BP)

# Get cumulative positions
chr_info = df_plot %>%
    group_by(Chr) %>%
    summarize(chr_len = max(BP)) %>%
    mutate(tot = cumsum(chr_len) - chr_len)

df_plot = df_plot %>%
    left_join(chr_info, by = "Chr") %>%
    mutate(BP_cum = BP + tot)

# Get chromosome center for x-axis label positions
axis_df = chr_info %>%
    mutate(center = tot + chr_len / 2)

# Pick top distinct significant genes by p-value (higher y)
top_genes = df_plot %>%
    filter(negative_log10_padj_snp > -log10(FDR_THRESH)) %>%
    arrange(desc(negative_log10_padj_snp)) %>%
    distinct(Probe, .keep_all = TRUE) %>%
    top_n(25, negative_log10_padj_snp)

yint = min(df_plot$negative_log10_padj_snp[df_plot$is_significant_snp])

cat(paste0("\n MAKING MANHATTAN PLOT\n"))
# Manhattan plot
m = (
    ggplot(df_plot, aes(x = BP_cum, y = negative_log10_padj_snp)) +
    ggtitle(paste0("eQTL Manhattan Plot:\n", cell_class)) +

    # Alternating background rectangles for chromosomes
    geom_rect(data = axis_df,
            aes(xmin = tot, xmax = tot + chr_len, ymin = -Inf, ymax = Inf, fill = as.factor(as.numeric(Chr) %% 2)),
            alpha = 0.1, inherit.aes = FALSE) +

    # Points
    geom_point(aes(color = is_significant_snp), alpha = 0.75, size = 1.2) +

    # Horizontal threshold line
    geom_hline(yintercept = yint, linetype = "dotted", color = "red") +

    # Labels for top 20 genes
    geom_label_repel(
        data = top_genes,
        aes(label = Probe),
        size = 6,
        box.padding = 0.5,
        nudge_y = 1,
        segment.color = 'grey50',
        max.overlaps = Inf
    ) +

    # Customize scales and theme
    scale_x_continuous(labels = axis_df$Chr, breaks = axis_df$center) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("white", "grey75")) +
    scale_color_manual(values = rep(c("black", "steelblue"), length.out = length(unique(df_plot$Chr)))) +

    labs(x = "Chromosome", y = expression(-log[10]("BH adjusted p-value"))) +

    theme_minimal(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, size=14),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=18),
        plot.title = element_text(size = 22)
    )
)

print(m)
ggsave(
    file.path(plot_dir, paste0(slogan, "_manhattan.png")),
    plot=m, width=20, height=12, dpi=800)
