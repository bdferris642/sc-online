# 05-process-and-plot-osca-tsv.R — Apply two-stage FDR to OSCA output; save RDS; generate plots.
#
# Arguments:
#   --path / -p         (required) Path to eqtl_{cell_class}.tsv produced by step 3.
#                                  File must be named eqtl_{cell_class}.tsv; cell_class is extracted
#                                  from the filename as everything after the leading "eqtl_" prefix.
#   --out-dir / -o      (optional) Output directory. Defaults to dirname(--path).
#   --ensg-to-symbol /  (optional) ENSG→symbol CSV (ensembl_gene_id,hgnc_symbol);
#   -E                             default: /mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv
#
# Two-stage FDR procedure:
#   1. Within each gene (Probe/ensg_id): Bonferroni across all cis-SNPs → p_bonf
#   2. Minimum p_bonf per gene → BH across genes → padj_gene
#   3. eGenes: padj_gene < FDR_THRESH
#   4. Within eGenes: BH across all SNPs → padj_snp
#   5. Significant eSNPs: padj_gene < FDR_THRESH AND padj_snp < FDR_THRESH
#
# v2 RDS column schema:
#   SNP, Chr, BP, A1, A2, Freq — unchanged
#   ensg_id     — Ensembl gene ID (was Probe)
#   gene_symbol — HGNC symbol or NA (new; best-effort join from ensg_to_symbol.csv)
#   Probe_Chr, Probe_bp — gene chr + TSS (unchanged names)
#   Orientation — strand (unchanged)
#   b, SE, p, padj_gene, padj_snp — unchanged
#   is_significant_gene, is_significant_snp — unchanged
#
# Outputs (in out-dir):
#   eqtl_{cell_class}.rds       — full data frame with FDR columns
#   eqtl_{cell_class}_sig.rds   — significant rows only
#   plots/eqtl_{cell_class}_manhattan_capped.png   — y-axis capped at HARD_CAP
#   plots/eqtl_{cell_class}_manhattan_capped.svg
#   plots/eqtl_{cell_class}_manhattan_uncapped.png  — y-axis uncapped
#   plots/eqtl_{cell_class}_manhattan_uncapped.svg
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
    'path',           'p', 1, "character",
    'out-dir',        'o', 1, "character",
    'ensg-to-symbol', 'E', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
if (is.null(opt[["out-dir"]])) {
    OUT_DIR = dirname(PATH)
} else {
    OUT_DIR = opt[["out-dir"]]
}
ENSG_TO_SYMBOL_PATH = if (!is.null(opt[["ensg-to-symbol"]])) opt[["ensg-to-symbol"]] else
    "/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv"

# Load ENSG→symbol mapping (best-effort; non-fatal if absent)
ensg_to_sym_map = setNames(character(0), character(0))
if (file.exists(ENSG_TO_SYMBOL_PATH)) {
    ensg_to_sym = read.csv(ENSG_TO_SYMBOL_PATH) %>%
        filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%
        distinct(ensembl_gene_id, .keep_all = TRUE)
    ensg_to_sym_map = setNames(ensg_to_sym$hgnc_symbol, ensg_to_sym$ensembl_gene_id)
    cat(sprintf("Loaded %d ENSG→symbol mappings\n", length(ensg_to_sym_map)))
} else {
    cat("Warning: ensg_to_symbol CSV not found; gene_symbol column will be NA for all genes\n")
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

# v2: Probe column = ENSG ID; Gene column = ENSG (duplicate of Probe from .opi probe field)
# Rename Probe → ensg_id; drop Gene (redundant); add gene_symbol via join.
df = df %>%
    rename(ensg_id = Probe) %>%
    select(-any_of("Gene")) %>%
    mutate(gene_symbol = ifelse(
        ensg_id %in% names(ensg_to_sym_map),
        ensg_to_sym_map[ensg_id],
        NA_character_))
cat(sprintf("Mapped gene_symbol for %d / %d unique genes\n",
    sum(!is.na(unique(df$gene_symbol[df$ensg_id %in% names(ensg_to_sym_map)]))),
    length(unique(df$ensg_id))))

# Two-stage FDR: within-gene Bonferroni → BH across genes (padj_gene) → BH within eGenes (padj_snp).

cat(paste0("\n PERFORMING TWO-STAGE FDR\n"))
cat(paste0("\n Step 1: Find lead SNPs per gene\n"))

# Step 1: Bonferroni correction within each gene
df_probe = df %>%
    group_by(ensg_id) %>%
    mutate(p_bonf = p * n()) %>%
    mutate(p_bonf = ifelse(p_bonf > 1, 1, p_bonf)) %>%
    ungroup()

# Step 2: Get the minimum p_bonf per gene
df_min_p_probe = df_probe %>%
    group_by(ensg_id) %>%
    summarise(min_p_gene = min(p_bonf)) %>%
    ungroup()

# Step 3: BH correct these min_p_gene values across genes
df_min_p_probe = df_min_p_probe %>%
    mutate(padj_gene = p.adjust(min_p_gene, method = "BH"))

# Step 4: Identify eGenes (FDR < threshold)
eGenes = df_min_p_probe %>%
    filter(padj_gene < FDR_THRESH) %>%
    pull(ensg_id)

# Step 5: Apply BH FDR to SNPs within eGenes
df_eGenes = df %>%
    filter(ensg_id %in% eGenes) %>%
    mutate(padj_snp = p.adjust(p, method = "BH"))

# Step 6: Merge the min_p_gene and padj_gene back to the df
df = df %>%
    left_join(df_min_p_probe[, c("ensg_id", "min_p_gene", "padj_gene")], by = "ensg_id") %>%
    left_join(df_eGenes[, c("SNP", "ensg_id", "padj_snp")], by = c("SNP", "ensg_id")) %>%
    mutate(
        padj_snp = ifelse(is.na(padj_snp), 1, padj_snp),
        is_significant_gene = padj_gene < FDR_THRESH,
        is_significant_snp = padj_snp < FDR_THRESH
    )

df$negative_log10_p_raw        = -log10(df$p)            # raw SNP p; never capped; used for ranking

df$negative_log10_padj_snp_raw = -log10(df$padj_snp)     # uncapped; used for uncapped plot
df$negative_log10_padj_snp     = pmin(df$negative_log10_padj_snp_raw, HARD_CAP)  # capped plot

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

# Pick top 25 annotatable genes for labeling.
# Rules:
#   1. Must be a significant eSNP (padj_snp < FDR_THRESH).
#   2. Exclude genes without a resolved symbol (ENSG fallback), LINCs, antisense (-AS),
#      divergent transcripts (-DT), and readthrough (-IT/-OT) loci — these clutter plots
#      without conveying biological meaning.
#   3. Rank by raw uncapped -log10(p) (not the capped padj_snp column) to break ties
#      at the hard cap and always surface the most extreme hits.
#   4. One representative SNP per gene (lowest raw p), top 25.
JUNK_PATTERN = "^LINC|^ENSG|^AC[0-9]|^AL[0-9]|^AP[0-9]|-AS[0-9]*$|-DT$|-IT[0-9]*$|-OT$"
top_genes = df_plot %>%
    filter(is_significant_snp) %>%
    filter(!is.na(gene_symbol)) %>%
    filter(!grepl(JUNK_PATTERN, gene_symbol)) %>%
    arrange(desc(negative_log10_p_raw)) %>%
    distinct(ensg_id, .keep_all = TRUE) %>%
    slice_head(n = 25) %>%
    mutate(plot_label = gene_symbol)

cat(sprintf("  Annotating %d genes in Manhattan plot\n", nrow(top_genes)))

sig_vals = df_plot$negative_log10_padj_snp[df_plot$is_significant_snp]
yint = if (length(sig_vals) > 0) min(sig_vals) else NA_real_

# Shared theme and layer factory for the two Manhattan plots
manhattan_layers = function(df_plot, y_col, y_label, title_suffix, top_genes_y, yint_val) {
    # y-headroom: 35% above the highest label point so repel has space to work
    y_max_data  = max(df_plot[[y_col]], na.rm = TRUE)
    y_max_label = if (nrow(top_genes_y) > 0) max(top_genes_y[[y_col]], na.rm = TRUE) else y_max_data
    y_ceiling   = max(y_max_data, y_max_label) * 1.35

    ggplot(df_plot, aes_string(x = "BP_cum", y = y_col)) +
    ggtitle(paste0("eQTL Manhattan Plot", title_suffix, ":\n", cell_class)) +

    geom_rect(data = axis_df,
        aes(xmin = tot, xmax = tot + chr_len, ymin = -Inf, ymax = Inf,
            fill = as.factor(as.numeric(Chr) %% 2)),
        alpha = 0.1, inherit.aes = FALSE) +

    geom_point(aes(color = is_significant_snp), alpha = 0.75, size = 1.2) +

    { if (!is.na(yint_val))
        geom_hline(yintercept = yint_val, linetype = "dotted", color = "red")
      else geom_blank() } +

    geom_text_repel(
        data = top_genes_y,
        aes_string(label = "plot_label"),
        size          = 3.8,
        box.padding   = 0.7,
        point.padding = 0.4,
        force         = 4,
        force_pull    = 0.3,
        min.segment.length = 0,
        segment.color = "grey40",
        segment.size  = 0.35,
        max.overlaps  = Inf,
        max.iter      = 20000,
        seed          = 42
    ) +

    scale_x_continuous(labels = axis_df$Chr, breaks = axis_df$center) +
    scale_y_continuous(limits = c(0, y_ceiling),
                       expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = c("0" = "white", "1" = "grey75")) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "steelblue")) +

    labs(x = "Chromosome", y = y_label) +

    theme_minimal(base_size = 14) +
    theme(
        legend.position    = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x  = element_text(angle = 0, vjust = 0.5, size = 14),
        axis.text.y  = element_text(size = 14),
        axis.title   = element_text(size = 18),
        plot.title   = element_text(size = 22)
    )
}

cat(paste0("\n MAKING MANHATTAN PLOTS (capped and uncapped)\n"))

# --- Plot 1: capped at HARD_CAP ---
# top_genes y-position for the capped plot uses the capped column
top_genes_capped = top_genes %>%
    mutate(!!sym("negative_log10_padj_snp") := pmin(negative_log10_padj_snp_raw, HARD_CAP))

m_capped = manhattan_layers(
    df_plot      = df_plot,
    y_col        = "negative_log10_padj_snp",
    y_label      = bquote(-log[10]("BH adj. p") ~ "(capped at" ~ .(HARD_CAP) * ")"),
    title_suffix = paste0(" (capped -log10p ≤ ", HARD_CAP, ")"),
    top_genes_y  = top_genes_capped,
    yint_val     = yint
)
print(m_capped)
ggsave(
    file.path(plot_dir, paste0(slogan, "_manhattan_capped.png")),
    plot = m_capped, width = 20, height = 12, dpi = 800)
ggsave(
    file.path(plot_dir, paste0(slogan, "_manhattan_capped.svg")),
    plot = m_capped, width = 20, height = 12)

# --- Plot 2: uncapped ---
yint_uncapped = if (length(sig_vals) > 0) {
    min(-log10(df_plot$padj_snp[df_plot$is_significant_snp]), na.rm = TRUE)
} else NA_real_

m_uncapped = manhattan_layers(
    df_plot      = df_plot %>% mutate(negative_log10_padj_snp_raw = -log10(padj_snp)),
    y_col        = "negative_log10_padj_snp_raw",
    y_label      = expression(-log[10]("BH adj. p")),
    title_suffix = " (uncapped)",
    top_genes_y  = top_genes,   # already has negative_log10_padj_snp_raw
    yint_val     = yint_uncapped
)
print(m_uncapped)
ggsave(
    file.path(plot_dir, paste0(slogan, "_manhattan_uncapped.png")),
    plot = m_uncapped, width = 20, height = 12, dpi = 800)
ggsave(
    file.path(plot_dir, paste0(slogan, "_manhattan_uncapped.svg")),
    plot = m_uncapped, width = 20, height = 12)
