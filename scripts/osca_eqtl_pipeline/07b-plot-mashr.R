#!/usr/bin/env Rscript
# 07b-plot-mashr.R — Visualization of mashr eQTL results
# Produces 10 plots (PNG + SVG + companion CSV each) from mashr output.
#
# Usage (self-activating — no manual env activation needed):
#   Rscript 07b-plot-mashr.R --mashr-sig <path> --eqtl-dir <dir> \
#     --eqtl-long <path> --out-dir <dir> [--lfsr-thresh 0.05] \
#     [--nom-thresh 0.05] [--gene-loc <path>]

# ── Self-activation block ─────────────────────────────────────────────────────
local({
  args <- commandArgs(trailingOnly = FALSE)
  script_flag <- "--file="
  script_path <- sub(script_flag, "", args[grepl(script_flag, args)])
  if (length(script_path) == 0) return()  # running interactively

  script_path <- normalizePath(script_path, mustWork = FALSE)
  sandbox     <- dirname(dirname(script_path))  # cc-sandbox root
  env_rscript <- file.path(sandbox, "micromamba_root", "envs", "osca-venv",
                           "bin", "Rscript")
  env_r_home  <- file.path(sandbox, "micromamba_root", "envs", "osca-venv",
                           "lib", "R")

  if (!identical(normalizePath(R.home()), normalizePath(env_r_home)) &&
      file.exists(env_rscript)) {
    message("[07b-plot-mashr] Re-launching under osca-venv Rscript …")
    # Pass all original arguments through
    orig_args <- commandArgs(trailingOnly = TRUE)
    cmd <- c(env_rscript, "--no-save", "--no-restore",
             "--no-site-file", "--no-init-file",
             script_path, orig_args)
    code <- system2(cmd[1], shQuote(cmd[-1]))
    quit(save = "no", status = code)
  }
})

suppressPackageStartupMessages({
  library(optparse)
  library(mashr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(pheatmap)
  library(UpSetR)
  library(grid)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-s", "--mashr-sig"),  type = "character", default = NULL,
              help = "Path to *__mash_results_sig.rds [required]"),
  make_option(c("-d", "--eqtl-dir"),   type = "character", default = NULL,
              help = "Dir with eqtl_{celltype}.rds + _sig.rds [required]"),
  make_option(c("-l", "--eqtl-long"),  type = "character", default = NULL,
              help = "Path to eqtl_present_in_all.rds [required]"),
  make_option(c("-o", "--out-dir"),    type = "character", default = NULL,
              help = "Output directory [required]"),
  make_option(c("-t", "--lfsr-thresh"),type = "double",    default = 0.05,
              help = "LFSR significance threshold [default %default]"),
  make_option(c("-n", "--nom-thresh"), type = "double",    default = 0.05,
              help = "Nominal p threshold for TSS plot [default %default]"),
  make_option(c("-g", "--gene-loc"),   type = "character", default = NULL,
              help = "Path to gene_loc_v2.txt for volcano plot labels (optional)")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("mashr-sig", "eqtl-dir", "eqtl-long", "out-dir")
missing_args <- required[!required %in% names(opt) | sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing_args) > 0)
  stop("Missing required argument(s): ", paste("--", missing_args, sep = "", collapse = ", "))

mashr_sig_path <- opt[["mashr-sig"]]
eqtl_dir       <- opt[["eqtl-dir"]]
eqtl_long_path <- opt[["eqtl-long"]]
out_dir        <- opt[["out-dir"]]
lfsr_thresh    <- opt[["lfsr-thresh"]]
nom_thresh     <- opt[["nom-thresh"]]
gene_loc_path  <- opt[["gene-loc"]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
message("[07b-plot-mashr] Loading mashr results: ", mashr_sig_path)
m <- readRDS(mashr_sig_path)

lfsr_mat <- get_lfsr(m)      # rows = SNP_probe, cols = cell_class
pm_mat   <- get_pm(m)
cell_classes <- colnames(lfsr_mat)
snp_probes   <- rownames(lfsr_mat)

message("[07b-plot-mashr] ", nrow(lfsr_mat), " SNP-probe pairs × ",
        length(cell_classes), " cell classes: ",
        paste(cell_classes, collapse = ", "))

# ── Helper: save PNG + SVG + CSV ─────────────────────────────────────────────
save_gg <- function(p, stem, width = 8, height = 6) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width,
         height = height, dpi = 150)
  ggsave(file.path(out_dir, paste0(stem, ".svg")), p, width = width,
         height = height)
  message("  Saved: ", stem, ".{png,svg}")
}

save_csv <- function(df, stem) {
  write.csv(df, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)
  message("  Saved: ", stem, ".csv")
}

save_pheatmap <- function(p, stem, width = 7, height = 6) {
  png(file.path(out_dir, paste0(stem, ".png")),
      width = width * 150, height = height * 150, res = 150)
  grid.newpage(); grid.draw(p$gtable)
  dev.off()
  svg(file.path(out_dir, paste0(stem, ".svg")), width = width, height = height)
  grid.newpage(); grid.draw(p$gtable)
  dev.off()
  message("  Saved: ", stem, ".{png,svg}")
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot 01 — UpSet: LFSR-significant SNP-probes by cell class
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 01] UpSet: LFSR-sig SNP-probes by cell class")
sig_mat <- lfsr_mat < lfsr_thresh
membership_list <- lapply(cell_classes, function(cc) snp_probes[sig_mat[, cc]])
names(membership_list) <- cell_classes

# CSV: binary membership table
csv01 <- as.data.frame(sig_mat * 1)
csv01$SNP_probe <- snp_probes
save_csv(csv01, "01_upset_lfsr_sig")

# PNG only — UpSetR doesn't support SVG reliably
message("  Note: UpSetR does not support SVG; emitting PNG only for Plot 01.")
png(file.path(out_dir, "01_upset_lfsr_sig.png"), width = 1800, height = 1200, res = 150)
upset(fromList(membership_list),
      nsets       = length(membership_list),
      order.by    = "freq",
      text.scale  = 1.4,
      mainbar.y.label = "# SNP-probe pairs",
      sets.x.label    = "Total sig (per class)")
dev.off()
message("  Saved: 01_upset_lfsr_sig.png")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 02 — Histogram: # cell classes per significant eQTL
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 02] Histogram: # cell classes per sig eQTL")
n_sig <- rowSums(sig_mat)
df02  <- data.frame(
  n_classes = factor(n_sig, levels = 0:length(cell_classes)),
  SNP_probe = snp_probes
)
counts02 <- df02 %>% count(n_classes, name = "count")

p02 <- ggplot(counts02, aes(x = n_classes, y = count)) +
  geom_col(fill = "#4682B4") +
  labs(title = "Number of cell classes per significant eQTL",
       x = "# cell classes (LFSR < threshold)",
       y = "# SNP-probe pairs") +
  theme_classic(base_size = 14)

save_gg(p02, "02_n_classes_per_eqtl")
save_csv(counts02, "02_n_classes_per_eqtl")

# ─────────────────────────────────────────────────────────────────────────────
# Plots 03 & 04 — eSNPs and eGenes per cell class (total vs unique)
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 03/04] eSNP and eGene counts per cell class")

sig_files <- list.files(eqtl_dir, pattern = "eqtl_.+_sig\\.rds$", full.names = TRUE)
if (length(sig_files) == 0)
  stop("No eqtl_*_sig.rds files found in: ", eqtl_dir)

sig_data <- lapply(sig_files, readRDS)
# derive cell class from filename: eqtl_<celltype>_sig.rds
names(sig_data) <- sub("eqtl_(.+)_sig\\.rds$", "\\1", basename(sig_files))

# helper to compute total vs unique counts for a given column
make_count_df <- function(col) {
  all_vals  <- lapply(sig_data, function(d) unique(d[[col]]))
  # unique to a class = not in any other class
  count_df <- do.call(rbind, lapply(names(all_vals), function(cc) {
    others <- unique(unlist(all_vals[names(all_vals) != cc]))
    total  <- length(all_vals[[cc]])
    unique_n <- sum(!all_vals[[cc]] %in% others)
    data.frame(cell_class = cc, total = total, unique_count = unique_n,
               shared = total - unique_n, stringsAsFactors = FALSE)
  }))
  count_df
}

snp_counts  <- make_count_df("SNP")
gene_counts <- make_count_df("ensg_id")

make_stack_plot <- function(df, title, ylab) {
  plot_df <- df %>%
    select(cell_class, unique_count, shared) %>%
    pivot_longer(c(unique_count, shared), names_to = "type", values_to = "count") %>%
    mutate(type = factor(type, levels = c("unique_count", "shared"),
                         labels = c("Cell-class unique", "Shared")))
  ggplot(plot_df, aes(x = cell_class, y = count, fill = type)) +
    geom_col() +
    scale_fill_manual(values = c("Cell-class unique" = "#E07B54",
                                 "Shared"            = "#4682B4")) +
    labs(title = title, x = NULL, y = ylab, fill = NULL) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

p03 <- make_stack_plot(snp_counts,  "eSNPs per cell class", "# unique eSNPs")
p04 <- make_stack_plot(gene_counts, "eGenes per cell class", "# unique eGenes")

save_gg(p03, "03_esnp_per_cellclass")
save_csv(snp_counts, "03_esnp_per_cellclass")

save_gg(p04, "04_egene_per_cellclass")
save_csv(gene_counts, "04_egene_per_cellclass")

# ─────────────────────────────────────────────────────────────────────────────
# Plots 05a/05b — Pairwise sharing heatmaps
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 05a] Pairwise sharing by sign")
ps_sign <- get_pairwise_sharing(m, factor = 0)
p05a <- pheatmap(ps_sign,
                 cluster_rows   = TRUE, cluster_cols = TRUE,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = "Pairwise sharing (same sign)",
                 color = colorRampPalette(c("white", "#4682B4"))(50),
                 fontsize = 12, silent = TRUE)
save_pheatmap(p05a, "05a_pairwise_sharing_sign")
save_csv(as.data.frame(ps_sign), "05a_pairwise_sharing_sign")

message("[Plot 05b] Pairwise sharing by magnitude")
ps_mag <- get_pairwise_sharing(m, FUN = abs)
p05b <- pheatmap(ps_mag,
                 cluster_rows   = TRUE, cluster_cols = TRUE,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = "Pairwise sharing (magnitude)",
                 color = colorRampPalette(c("white", "#E07B54"))(50),
                 fontsize = 12, silent = TRUE)
save_pheatmap(p05b, "05b_pairwise_sharing_mag")
save_csv(as.data.frame(ps_mag), "05b_pairwise_sharing_mag")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 06 — TSS distance density by significance category
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 06] TSS distance density by significance category")
eqtl_long <- readRDS(eqtl_long_path)

# Pivot wide: one row per SNP_probe, columns = p value per cell class
p_wide <- eqtl_long %>%
  select(SNP_probe, cell_class, p, BP, Probe_bp) %>%
  distinct() %>%
  pivot_wider(names_from = cell_class, values_from = p,
              names_prefix = "p_", values_fn = min)

# For TSS dist use BP and Probe_bp from the first occurrence per SNP_probe
pos_df <- eqtl_long %>%
  select(SNP_probe, BP, Probe_bp) %>%
  distinct(SNP_probe, .keep_all = TRUE)

p_wide <- p_wide %>%
  select(-starts_with("BP"), -starts_with("Probe_bp")) %>%
  left_join(pos_df, by = "SNP_probe")

p_cols  <- paste0("p_", cell_classes)
lfsr_df <- as.data.frame(lfsr_mat)
lfsr_df$SNP_probe <- rownames(lfsr_mat)

df06 <- p_wide %>%
  left_join(lfsr_df, by = "SNP_probe") %>%
  mutate(
    lfsr_sig_any = if (length(cell_classes) > 1)
      rowSums(across(all_of(cell_classes), ~ . < lfsr_thresh)) > 0
    else
      .data[[cell_classes[1]]] < lfsr_thresh,
    nom_sig_any = if (length(p_cols) > 1)
      rowSums(across(all_of(p_cols), ~ !is.na(.) & . < nom_thresh)) > 0
    else
      !is.na(.data[[p_cols[1]]]) & .data[[p_cols[1]]] < nom_thresh,
    category = case_when(
       lfsr_sig_any &  nom_sig_any ~ "mashr sig & nom sig",
       lfsr_sig_any & !nom_sig_any ~ "mashr sig only",
      !lfsr_sig_any &  nom_sig_any ~ "nom sig only",
      TRUE                         ~ "not significant"
    ),
    cis_dist = abs(BP - Probe_bp) + 1  # +1 to avoid log(0)
  ) %>%
  filter(!is.na(cis_dist))

csv06 <- df06 %>% select(SNP_probe, cis_dist, category)
save_csv(csv06, "06_tss_distance")

cat_colors <- c(
  "mashr sig & nom sig" = "#D62728",
  "mashr sig only"      = "#FF7F0E",
  "nom sig only"        = "#1F77B4",
  "not significant"     = "#AAAAAA"
)

p06 <- ggplot(df06, aes(x = cis_dist, color = category)) +
  geom_density(linewidth = 0.9) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = cat_colors) +
  labs(title = "TSS distance by significance category",
       x = "Cis distance (bp, log10)", y = "Density", color = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom")

save_gg(p06, "06_tss_distance", width = 9, height = 6)

# ─────────────────────────────────────────────────────────────────────────────
# Plot 07 — ED_PCA_1 covariance heatmap
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 07] ED_PCA_1 covariance heatmap")
fitted_g  <- get_fitted_g(m)
ed_pca1   <- fitted_g$Ulist[["ED_PCA_1"]]

if (is.null(ed_pca1)) {
  message("  WARNING: ED_PCA_1 not found in fitted_g$Ulist — skipping Plot 07.")
  message("  Available components: ", paste(names(fitted_g$Ulist), collapse = ", "))
} else {
  cor_mat <- cov2cor(ed_pca1)
  rownames(cor_mat) <- cell_classes
  colnames(cor_mat) <- cell_classes
  p07 <- pheatmap(cor_mat,
                  cluster_rows    = TRUE, cluster_cols = TRUE,
                  display_numbers = TRUE, number_format = "%.2f",
                  main    = "ED_PCA_1 correlation matrix",
                  color   = colorRampPalette(c("#D62728", "white", "#4682B4"))(51),
                  breaks  = seq(-1, 1, length.out = 52),
                  fontsize = 14, silent = TRUE)
  save_pheatmap(p07, "07_ed_pca1_heatmap")
  save_csv(as.data.frame(cor_mat), "07_ed_pca1_heatmap")
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot 08 — Mixture proportions bar chart
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 08] Mixture proportions")
pi_vals <- get_estimated_pi(m)
df08 <- data.frame(
  component = names(pi_vals),
  est_pi    = as.numeric(pi_vals),
  stringsAsFactors = FALSE
) %>% arrange(desc(est_pi))
df08$component <- factor(df08$component, levels = df08$component)

p08 <- ggplot(df08, aes(x = component, y = est_pi)) +
  geom_col(fill = "#4682B4") +
  labs(title = "Mashr mixture proportions",
       x = "Component", y = "Estimated proportion") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_gg(p08, "08_mixture_proportions", width = max(6, length(pi_vals) * 0.7), height = 5)
save_csv(df08, "08_mixture_proportions")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 09 — Pairwise PM scatter plots
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 09] Pairwise PM scatter plots")
pairs <- combn(cell_classes, 2, simplify = FALSE)

if (length(pairs) == 0) {
  message("  Only 1 cell class — skipping pairwise scatter.")
} else {
  pm_df <- as.data.frame(pm_mat)
  pm_df$SNP_probe <- snp_probes

  sig_df <- as.data.frame(sig_mat)
  sig_df$SNP_probe <- snp_probes

  plot_list <- lapply(pairs, function(pr) {
    cA <- pr[1]; cB <- pr[2]
    dat <- pm_df %>%
      select(SNP_probe, x = all_of(cA), y = all_of(cB)) %>%
      left_join(sig_df %>% select(SNP_probe,
                                   sigA = all_of(cA),
                                   sigB = all_of(cB)),
                by = "SNP_probe") %>%
      mutate(status = case_when(
         sigA &  sigB ~ "Both",
         sigA & !sigB ~ paste("Sig in", cA),
        !sigA &  sigB ~ paste("Sig in", cB),
        TRUE           ~ "Neither"
      ))
    ggplot(dat, aes(x = x, y = y, color = status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c(
        "Both"    = "#D62728", "Neither" = "#AAAAAA",
        setNames(c("#FF7F0E", "#1F77B4"), paste("Sig in", c(cA, cB)))
      ), drop = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      labs(title = paste("Posterior mean:", cA, "vs", cB),
           x = paste("PM (", cA, ")"), y = paste("PM (", cB, ")"),
           color = NULL) +
      theme_classic(base_size = 13) +
      theme(legend.position = "bottom")
  })

  # For 2 classes: one panel; for more: facet grid would need reshaping, just save multi-page
  n_pairs <- length(plot_list)
  ncols   <- min(3, n_pairs)
  nrows   <- ceiling(n_pairs / ncols)

  # assemble with cowplot / patchwork if available, else save individual pages
  combined <- tryCatch({
    if (requireNamespace("patchwork", quietly = TRUE)) {
      Reduce(`+`, plot_list) + patchwork::plot_layout(ncol = ncols)
    } else {
      NULL
    }
  }, error = function(e) NULL)

  # CSV: long format of all pairs
  csv09 <- do.call(rbind, lapply(pairs, function(pr) {
    cA <- pr[1]; cB <- pr[2]
    pm_df %>%
      select(SNP_probe, x = all_of(cA), y = all_of(cB)) %>%
      left_join(sig_df %>% select(SNP_probe,
                                   sigA = all_of(cA),
                                   sigB = all_of(cB)),
                by = "SNP_probe") %>%
      mutate(pair = paste(cA, "vs", cB),
             classA = cA, classB = cB) %>%
      rename(PM_A = x, PM_B = y, sig_A = sigA, sig_B = sigB)
  }))
  save_csv(csv09, "09_pm_pairwise_scatter")

  if (!is.null(combined)) {
    save_gg(combined, "09_pm_pairwise_scatter",
            width = ncols * 5, height = nrows * 5)
  } else {
    # Save PNG page-by-page using base graphics multi-page PNG
    png(file.path(out_dir, "09_pm_pairwise_scatter.png"),
        width = 900, height = 700, res = 150)
    for (pl in plot_list) print(pl)
    dev.off()
    svg(file.path(out_dir, "09_pm_pairwise_scatter.svg"),
        width = 6, height = 5)
    for (pl in plot_list) print(pl)
    dev.off()
    message("  Saved: 09_pm_pairwise_scatter.{png,svg} (multi-page)")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot 10 — Volcano plot per cell class
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 10] Volcano plots per cell class")

# Load gene symbols if --gene-loc provided
gene_sym <- NULL
if (!is.null(gene_loc_path) && file.exists(gene_loc_path)) {
  gene_loc <- read.table(gene_loc_path, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
  # expected columns: ensg_id, gene_symbol
  if (all(c("ensg_id", "gene_symbol") %in% colnames(gene_loc))) {
    gene_sym <- gene_loc %>%
      select(ensg_id, gene_symbol) %>%
      distinct(ensg_id, .keep_all = TRUE)
    message("  Gene symbols loaded: ", nrow(gene_sym), " entries")
  } else {
    message("  WARNING: gene_loc_v2.txt missing ensg_id/gene_symbol columns; ",
            "labels will use ensg_id.")
  }
} else if (!is.null(gene_loc_path)) {
  message("  WARNING: --gene-loc file not found: ", gene_loc_path)
}

# Build data frame: SNP_probe, cell_class, PM, LFSR, ensg_id, label
vol_base <- data.frame(
  SNP_probe  = snp_probes,
  stringsAsFactors = FALSE
)
# extract ensg_id from SNP_probe (format: rs*_ENSG*)
vol_base$ensg_id <- sub("^[^_]+_", "", snp_probes)

if (!is.null(gene_sym)) {
  vol_base <- left_join(vol_base, gene_sym, by = "ensg_id")
  vol_base$label <- ifelse(is.na(vol_base$gene_symbol) | vol_base$gene_symbol == "",
                           vol_base$ensg_id, vol_base$gene_symbol)
} else {
  vol_base$label <- vol_base$ensg_id
}

volcano_df <- do.call(rbind, lapply(cell_classes, function(cc) {
  data.frame(
    SNP_probe  = snp_probes,
    cell_class = cc,
    PM         = pm_mat[, cc],
    LFSR       = lfsr_mat[, cc],
    ensg_id    = vol_base$ensg_id,
    label      = vol_base$label,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(neg_log10_lfsr = -log10(pmax(LFSR, 1e-10)),
         sig            = LFSR < lfsr_thresh)

# Top 10 per class by smallest LFSR for labels
top_labels <- volcano_df %>%
  filter(sig) %>%
  group_by(cell_class) %>%
  slice_min(order_by = LFSR, n = 10) %>%
  ungroup()

p10 <- ggplot(volcano_df, aes(x = PM, y = neg_log10_lfsr, color = sig)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_text_repel(data = top_labels,
                  aes(label = label),
                  size = 3, max.overlaps = 20, color = "#D62728") +
  scale_color_manual(values = c("FALSE" = "#AAAAAA", "TRUE" = "#D62728"),
                     labels = c("Not sig", "LFSR sig")) +
  facet_wrap(~cell_class, scales = "free") +
  labs(title  = "Volcano plot per cell class",
       x      = "Posterior mean effect",
       y      = expression(-log[10](LFSR)),
       color  = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

n_cc    <- length(cell_classes)
v_ncols <- min(3, n_cc)
v_nrows <- ceiling(n_cc / v_ncols)
save_gg(p10, "10_volcano_per_cellclass",
        width = v_ncols * 5, height = v_nrows * 5)
save_csv(volcano_df %>% select(SNP_probe, cell_class, PM, LFSR, label, sig),
         "10_volcano_per_cellclass")

# ─────────────────────────────────────────────────────────────────────────────
message("\n[07b-plot-mashr] All done. Outputs in: ", out_dir)
