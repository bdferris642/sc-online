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
slogan = strsplit(basename(PATH), "\\.")[[1]][[1]]
cell_class = strsplit(slogan, "_")[[1]][[2]]
print(slogan)
print(cell_class)
plot_dir = file.path(OUT_DIR, "plots")
dir.create(OUT_DIR, showWarnings = F)
dir.create(plot_dir, showWarnings = F)

HARD_CAP = 40
FDR_THRESH = 0.05

cat(paste0("\n READING DataFrame from", PATH, "\n"))
df = read.table(PATH, header = T, sep = "\t")
df$BP = as.numeric(df$BP)
df$negative_log10_p = -log10(df$p)
df$negative_log10_p[df$negative_log10_p > HARD_CAP] = HARD_CAP

# perform two stage FDR
# (1) get min p-value for each gene and conduct FDR at gene level
# (2) Within each gene, perform FDR correction at the SNP level

cat(paste0("\n PERFORMING TWO-STAGE FDR\n"))
cat(paste0("\n Step 1: Find lead SNPs per gene\n"))
# Step 1: Find lead SNPs per gene
# gene_leads = (df 
#     %>% group_by(Probe)
#     %>% slice_min(order_by = p, n = 1, with_ties = FALSE)
#     %>% ungroup())

# # Step 2: Gene-level FDR on lead SNP p-values
# gene_leads = gene_leads %>% mutate(padj_gene = p.adjust(p, method = "bonferroni"))

# # Step 3: Merge gene-level FDR back into main df
# df = df %>% left_join(gene_leads %>% select(Probe, padj_gene), by = "Probe")

# # Step 4: Keep only SNPs from genes passing Stage 1 FDR threshold
# df_sig_genes = df %>% filter(padj_gene <= FDR_THRESH)


# cat(paste0("\n", nrow(df_sig_genes), " significant genes out of", nrow(gene_leads), "\n"))

# # Step 5: For those significant genes, apply SNP-level FDR correction within each gene
# df_sig_genes = (df_sig_genes 
#     %>% group_by(Probe) 
#     %>% mutate(padj_snp = p.adjust(p, method = "BH")) 
#     %>% ungroup())

# # Step 6: Merge SNP-level FDR back into main df
# df = df %>% left_join(df_sig_genes %>% select(SNP, Probe, padj_snp), by = c("SNP", "Probe"))

# # STEP 1: Get lead SNP (smallest nominal p) per gene
# gene_leads = df %>%
#     group_by(Probe) %>%
#     slice_min(order_by = p, n = 1, with_ties = FALSE) %>%
#     ungroup()

# # STEP 2: Apply gene-wise Bonferroni correction to lead SNPs
# gene_leads = gene_leads %>%
#     rowwise() %>%
#     mutate(min_p_bonf_within_gene = min(p * sum(df$Probe == Probe), 1)) %>%
#     ungroup()

# # STEP 3: Apply BH FDR to the gene-level Bonferroni p-values
# gene_leads = gene_leads %>%
#     mutate(padj_gene = p.adjust(min_p_bonf_within_gene, method = "BH"))

# # STEP 4: Identify significant eGenes
# signif_genes = gene_leads %>% filter(padj_gene < FDR_THRESH)

# # STEP 5: Apply BH correction across *all* gene–SNP pairs
# df = df %>% mutate(padj_snp = p.adjust(p, method = "BH"))

# # STEP 6: Determine nominal P threshold from significant gene–SNP pairs
# signif_pairs = df %>% filter(padj_snp < FDR_THRESH)
# nominal_p_thresh = max(signif_pairs$p)

# # STEP 7: Declare all SNPs with nominal P ≤ threshold as significant eSNPs
# signif_esnps = df %>%
#     filter(p <= nominal_p_thresh)

# # STEP 8: Set NA-adjusted p-values to 1 (optional, for plotting)
# df$padj_snp[is.na(df$padj_snp)] = 1
# df$negative_log10_padj_snp = -log10(df$padj_snp)
# df$negative_log10_padj_snp[df$negative_log10_padj_snp > HARD_CAP] = HARD_CAP

# # Add gene-level FDR values to df by joining from gene_leads
# df = df %>% left_join(gene_leads %>% select(Probe, padj_gene), by = "Probe")

# df$negative_log10_padj_gene = -log10(df$padj_gene)
# df$negative_log10_padj_gene[df$negative_log10_padj_gene > HARD_CAP] = HARD_CAP

# df$is_significant_gene = df$padj_gene < FDR_THRESH
# df$is_significant_snp = df$padj_snp < FDR_THRESH

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
    + ggtitle(paste0("GTEX EQTL Gene Minimum P-value Histogram:\n", cell_class))
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
    + ggtitle(paste0("GTEX EQTL Gene Adj. P-value Histogram:\n", cell_class))
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
    + ggtitle(paste0("GTEX EQTL P-value Histogram:\n", cell_class))
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
    + ggtitle(paste0("GTEX EQTL SNP Adj. P-value Histogram:\n", cell_class))
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
# Prepare df by calculating chromosome midpoints and cumulative base pair positions
chr_order <- as.character(1:22)
df = df %>%
    mutate(Chr = factor(Chr, levels=chr_order)) %>%
    arrange(Chr, BP)

df_sig = df[df$is_significant_gene & df$is_significant_snp,]
df_insig = df[!df$is_significant_gene | !df$is_significant_snp,]
# take 1/50th of the insignificant data
df_insig_subset = df_insig[sample(nrow(df_insig), nrow(df_insig) / 50), ]
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
    ggtitle(paste0("GTEX EQTL Manhattan Plot:\n", cell_class)) +
  
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