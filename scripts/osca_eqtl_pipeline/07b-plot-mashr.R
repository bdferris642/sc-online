#!/usr/bin/env Rscript
# 07b-plot-mashr.R — Visualization of mashr eQTL results.
#
# Usage (self-activating — no manual env activation needed):
#   Rscript 07b-plot-mashr.R --mashr-sig <path> --eqtl-dir <dir> \
#     --eqtl-long <path> --out-dir <dir> [--mashr-random <path>] \
#     [--lfsr-thresh 0.05] [--nom-thresh 0.05] [--gene-loc <path>]

# ── Self-activation block ─────────────────────────────────────────────────────
local({
  args <- commandArgs(trailingOnly = FALSE)
  script_flag <- "--file="
  script_path <- sub(script_flag, "", args[grepl(script_flag, args)])
  if (length(script_path) == 0) return()

  script_path <- normalizePath(script_path, mustWork = FALSE)
  sandbox     <- dirname(dirname(script_path))
  env_rscript <- file.path(sandbox, "micromamba_root", "envs", "osca-venv", "bin", "Rscript")
  env_r_home  <- file.path(sandbox, "micromamba_root", "envs", "osca-venv", "lib", "R")

  if (!identical(normalizePath(R.home()), normalizePath(env_r_home)) &&
      file.exists(env_rscript)) {
    message("[07b-plot-mashr] Re-launching under osca-venv Rscript …")
    orig_args <- commandArgs(trailingOnly = TRUE)
    cmd <- c(env_rscript, "--no-save", "--no-restore",
             "--no-site-file", "--no-init-file", script_path, orig_args)
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
  library(tibble)
  library(ggrepel)
  library(pheatmap)
  library(UpSetR)
  library(grid)
  library(scales)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-s", "--mashr-sig"),    type = "character", default = NULL,
              help = "Path to *__mash_results_sig.rds [required]"),
  make_option(c("-r", "--mashr-random"), type = "character", default = NULL,
              help = "Path to *__mash_results_random.rds [auto-derived from --mashr-sig if omitted]"),
  make_option(c("-d", "--eqtl-dir"),     type = "character", default = NULL,
              help = "Dir with eqtl_{celltype}.rds + _sig.rds [required]"),
  make_option(c("-l", "--eqtl-long"),    type = "character", default = NULL,
              help = "Path to eqtl_present_in_all.rds [required]"),
  make_option(c("-o", "--out-dir"),      type = "character", default = NULL,
              help = "Output directory [required]"),
  make_option(c("-t", "--lfsr-thresh"),  type = "double",    default = 0.05,
              help = "LFSR significance threshold [default %default]"),
  make_option(c("-n", "--nom-thresh"),   type = "double",    default = 0.05,
              help = "Nominal padj_snp threshold for TSS and LFSR density plots [default %default]"),
  make_option(c("-g", "--gene-loc"),     type = "character", default = NULL,
              help = "Path to gene_loc_v2.txt for volcano plot labels (optional)")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("mashr-sig", "eqtl-dir", "eqtl-long", "out-dir")
missing_args <- required[sapply(required, function(x) is.null(opt[[x]]))]
if (length(missing_args) > 0)
  stop("Missing required argument(s): ", paste("--", missing_args, sep = "", collapse = ", "))

mashr_sig_path    <- opt[["mashr-sig"]]
mashr_random_path <- opt[["mashr-random"]]
if (is.null(mashr_random_path))
  mashr_random_path <- sub("_sig\\.rds$", "_random.rds", mashr_sig_path)
eqtl_dir       <- opt[["eqtl-dir"]]
eqtl_long_path <- opt[["eqtl-long"]]
out_dir        <- opt[["out-dir"]]
lfsr_thresh    <- opt[["lfsr-thresh"]]
nom_thresh     <- opt[["nom-thresh"]]
gene_loc_path  <- opt[["gene-loc"]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load mashr sig results ────────────────────────────────────────────────────
message("[07b-plot-mashr] Loading mashr sig results: ", mashr_sig_path)
m          <- readRDS(mashr_sig_path)
lfsr_mat   <- get_lfsr(m)
pm_mat     <- get_pm(m)
cell_classes <- colnames(lfsr_mat)
snp_probes   <- rownames(lfsr_mat)
sig_mat      <- lfsr_mat < lfsr_thresh

message("[07b-plot-mashr] ", nrow(lfsr_mat), " SNP-probe pairs × ",
        length(cell_classes), " cell classes: ", paste(cell_classes, collapse = ", "))

# ── Load mashr random results (for LFSR density plot) ────────────────────────
mash_random <- NULL
if (file.exists(mashr_random_path)) {
  message("[07b-plot-mashr] Loading mashr random results: ", mashr_random_path)
  mash_random <- readRDS(mashr_random_path)
} else {
  message("  WARNING: mashr random RDS not found: ", mashr_random_path,
          " — LFSR random density plot will be skipped.")
}

# ── Load eqtl_present_in_all_sig_in_one (for LFSR density plots) ─────────────
sig_in_one_path <- file.path(eqtl_dir, "eqtl_present_in_all_sig_in_one.rds")
sig_in_one <- NULL
if (file.exists(sig_in_one_path)) {
  message("[07b-plot-mashr] Loading sig-in-one: ", sig_in_one_path)
  sig_in_one <- readRDS(sig_in_one_path)
} else {
  message("  WARNING: eqtl_present_in_all_sig_in_one.rds not found in eqtl-dir — ",
          "LFSR density plots will be skipped.")
}

# ── Per-cell-class color palette (consistent across all plots) ────────────────
color_vec <- setNames(hue_pal()(length(cell_classes)), cell_classes)

# ── Helper: save PNG + CSV ────────────────────────────────────────────────────
save_gg <- function(p, stem, width = 8, height = 6) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = width,
         height = height, dpi = 150)
  message("  Saved: ", stem, ".png")
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
  message("  Saved: ", stem, ".png")
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot 01 — UpSet: LFSR-significant SNP-probes by cell class
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 01] UpSet: LFSR-sig SNP-probes by cell class")
membership_list <- lapply(cell_classes, function(cc) snp_probes[sig_mat[, cc]])
names(membership_list) <- cell_classes

csv01 <- as.data.frame(sig_mat * 1)
csv01$SNP_probe <- snp_probes
save_csv(csv01, "01_upset_lfsr_sig")

png(file.path(out_dir, "01_upset_lfsr_sig.png"), width = 1800, height = 1200, res = 150)
upset(fromList(membership_list),
      nsets           = length(membership_list),
      order.by        = "freq",
      text.scale      = 1.4,
      show.numbers    = "no",          # suppress counts on bar tops
      mainbar.y.label = "# SNP-probe pairs",
      sets.x.label    = "Total sig (per class)")
dev.off()
message("  Saved: 01_upset_lfsr_sig.png")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 02 — Bar: # cell classes with LFSR-sig per SNP-probe (strong-signal set)
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 02] Bar: # cell classes LFSR-sig per SNP-probe")
n_sig    <- rowSums(sig_mat)
counts02 <- data.frame(n_classes = factor(n_sig, levels = 0:length(cell_classes))) %>%
  count(n_classes, name = "count")

p02 <- ggplot(counts02, aes(x = n_classes, y = count)) +
  geom_col(fill = "#4682B4") +
  labs(
    title = sprintf(
      "LFSR significance breadth across cell classes\n(strong-signal SNP-probe set; LFSR < %.2f)",
      lfsr_thresh),
    x = "# cell classes with LFSR < threshold",
    y = "# SNP-probe pairs (strong-signal set)") +
  theme_classic(base_size = 14)

save_gg(p02, "02_n_classes_per_eqtl")
save_csv(counts02, "02_n_classes_per_eqtl")

# ─────────────────────────────────────────────────────────────────────────────
# Plots 03 & 04 — eSNPs and eGenes per cell class
# Excludes present_in_all* files — per-cell-class _sig.rds only.
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 03/04] eSNP and eGene counts per cell class")

sig_files <- list.files(eqtl_dir, pattern = "eqtl_.+_sig\\.rds$", full.names = TRUE)
sig_files <- sig_files[!grepl("present_in_all", sig_files)]
if (length(sig_files) == 0)
  stop("No per-cell-class eqtl_*_sig.rds files found in: ", eqtl_dir)

sig_data <- lapply(sig_files, readRDS)
names(sig_data) <- sub("eqtl_(.+)_sig\\.rds$", "\\1", basename(sig_files))

make_count_df <- function(col) {
  all_vals <- lapply(sig_data, function(d) unique(d[[col]]))
  do.call(rbind, lapply(names(all_vals), function(cc) {
    others   <- unique(unlist(all_vals[names(all_vals) != cc]))
    total    <- length(all_vals[[cc]])
    unique_n <- sum(!all_vals[[cc]] %in% others)
    data.frame(cell_class = cc, total = total, unique_count = unique_n,
               shared = total - unique_n, stringsAsFactors = FALSE)
  }))
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
    scale_fill_manual(values = c("Cell-class unique" = "#E07B54", "Shared" = "#4682B4")) +
    labs(title = title, x = NULL, y = ylab, fill = NULL) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

p03 <- make_stack_plot(snp_counts,  "eSNPs per cell class",  "# unique eSNPs")
p04 <- make_stack_plot(gene_counts, "eGenes per cell class", "# unique eGenes")
save_gg(p03, "03_esnp_per_cellclass");  save_csv(snp_counts,  "03_esnp_per_cellclass")
save_gg(p04, "04_egene_per_cellclass"); save_csv(gene_counts, "04_egene_per_cellclass")

# ─────────────────────────────────────────────────────────────────────────────
# Plots 05a/05b — Pairwise sharing heatmaps
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 05a] Pairwise sharing by sign")
ps_sign <- get_pairwise_sharing(m, factor = 0)
p05a <- pheatmap(ps_sign, cluster_rows = TRUE, cluster_cols = TRUE,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = "Pairwise sharing (same sign)",
                 color = colorRampPalette(c("white", "#4682B4"))(50),
                 fontsize = 12, silent = TRUE)
save_pheatmap(p05a, "05a_pairwise_sharing_sign")
save_csv(as.data.frame(ps_sign), "05a_pairwise_sharing_sign")

message("[Plot 05b] Pairwise sharing by magnitude")
ps_mag <- get_pairwise_sharing(m, FUN = abs)
p05b <- pheatmap(ps_mag, cluster_rows = TRUE, cluster_cols = TRUE,
                 display_numbers = TRUE, number_format = "%.2f",
                 main = "Pairwise sharing (magnitude)",
                 color = colorRampPalette(c("white", "#E07B54"))(50),
                 fontsize = 12, silent = TRUE)
save_pheatmap(p05b, "05b_pairwise_sharing_mag")
save_csv(as.data.frame(ps_mag), "05b_pairwise_sharing_mag")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 06 — TSS distance density by significance category
# Uses per-SNP_probe summaries (no wide pivot) to keep memory manageable.
# Categories:
#   mashr sig & nom sig  — in lfsr_mat with LFSR-sig in ≥1 class AND raw p-sig in eqtl_long
#   mashr sig only       — in lfsr_mat with LFSR-sig in ≥1 class AND NOT raw p-sig
#   nom sig only         — NOT in lfsr_mat LFSR-sig AND raw p-sig in eqtl_long
#   not significant      — neither
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 06] TSS distance density by significance category")
eqtl_long <- readRDS(eqtl_long_path)

# Summarise eqtl_long per SNP_probe: nominal sig in any class, and cis positions.
message("  Summarising eqtl_long per SNP_probe ...")
pos_nom <- eqtl_long %>%
  group_by(SNP_probe) %>%
  summarise(
    BP       = first(BP),
    Probe_bp = first(Probe_bp),
    nom_sig_any = any(p < nom_thresh, na.rm = TRUE),
    .groups = "drop"
  )
rm(eqtl_long); gc()

# mashr LFSR sig: any cell class significant
lfsr_sig_df <- data.frame(
  SNP_probe    = snp_probes,
  lfsr_sig_any = rowSums(sig_mat) > 0
)

df06 <- pos_nom %>%
  left_join(lfsr_sig_df, by = "SNP_probe") %>%
  mutate(
    lfsr_sig_any = !is.na(lfsr_sig_any) & lfsr_sig_any,
    category = case_when(
       lfsr_sig_any &  nom_sig_any ~ "mashr sig & nom sig",
       lfsr_sig_any & !nom_sig_any ~ "mashr sig only",
      !lfsr_sig_any &  nom_sig_any ~ "nom sig only",
      TRUE                         ~ "not significant"
    ),
    cis_dist = abs(BP - Probe_bp) + 1
  ) %>%
  filter(!is.na(cis_dist))

message("  Category counts:\n",
        paste(capture.output(print(table(df06$category))), collapse = "\n"))

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
  scale_x_log10(labels = comma) +
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
fitted_g <- get_fitted_g(m)

make_cov_heatmap <- function(mat_name, stem, plot_title) {
  u <- fitted_g$Ulist[[mat_name]]
  if (is.null(u)) {
    message("  WARNING: ", mat_name, " not found — skipping.")
    message("  Available: ", paste(names(fitted_g$Ulist), collapse = ", "))
    return(invisible(NULL))
  }
  cor_mat <- cov2cor(u)
  rownames(cor_mat) <- cell_classes; colnames(cor_mat) <- cell_classes
  ph <- pheatmap(cor_mat, cluster_rows = TRUE, cluster_cols = TRUE,
                 display_numbers = TRUE, number_format = "%.2f",
                 main    = plot_title,
                 color   = colorRampPalette(c("#D62728", "white", "#4682B4"))(51),
                 breaks  = seq(-1, 1, length.out = 52),
                 fontsize = 14, silent = TRUE)
  save_pheatmap(ph, stem)
  save_csv(as.data.frame(cor_mat), stem)
}

make_cov_heatmap("ED_PCA_1",  "07a_ed_pca1_heatmap",  "ED_PCA_1 correlation matrix")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 07b — ED_tPCA correlation heatmap
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 07b] ED_tPCA covariance heatmap")
# ED_tPCA may be named "ED_tPCA_1" or similar; find any tPCA component.
tpca_name <- grep("tPCA", names(fitted_g$Ulist), value = TRUE)[1]
if (!is.na(tpca_name)) {
  make_cov_heatmap(tpca_name, "07b_ed_tpca_heatmap",
                   paste0(tpca_name, " correlation matrix"))
} else {
  message("  WARNING: no tPCA component found in fitted_g$Ulist — skipping 07b.")
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot 08 — Mixture proportions bar chart
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 08] Mixture proportions")
pi_vals <- get_estimated_pi(m)
df08 <- data.frame(component = names(pi_vals), est_pi = as.numeric(pi_vals),
                   stringsAsFactors = FALSE) %>% arrange(desc(est_pi))
df08$component <- factor(df08$component, levels = df08$component)

p08 <- ggplot(df08, aes(x = component, y = est_pi)) +
  geom_col(fill = "#4682B4") +
  labs(title = "Mashr mixture proportions", x = "Component", y = "Estimated proportion") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_gg(p08, "08_mixture_proportions", width = max(6, length(pi_vals) * 0.7), height = 5)
save_csv(df08, "08_mixture_proportions")

# ─────────────────────────────────────────────────────────────────────────────
# Plot 09 — Pairwise PM scatter plots (one file per pair)
# Title includes sign consistency among union-significant SNP-probes.
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 09] Pairwise PM scatter plots")
pairs <- combn(cell_classes, 2, simplify = FALSE)

if (length(pairs) == 0) {
  message("  Only 1 cell class — skipping pairwise scatter.")
} else {
  pm_df  <- as.data.frame(pm_mat);  pm_df$SNP_probe  <- snp_probes
  sig_df <- as.data.frame(sig_mat); sig_df$SNP_probe <- snp_probes
  csv09_rows <- list()

  for (pr in pairs) {
    cA <- pr[1]; cB <- pr[2]
    safe_cA <- gsub("[^A-Za-z0-9]", "_", cA)
    safe_cB <- gsub("[^A-Za-z0-9]", "_", cB)
    stem <- paste0("09_pm_", safe_cA, "_vs_", safe_cB, "_scatter")

    dat <- pm_df %>%
      select(SNP_probe, x = all_of(cA), y = all_of(cB)) %>%
      left_join(sig_df %>% select(SNP_probe, sigA = all_of(cA), sigB = all_of(cB)),
                by = "SNP_probe") %>%
      mutate(status = case_when(
         sigA &  sigB ~ "Both",
         sigA & !sigB ~ paste("Sig in", cA),
        !sigA &  sigB ~ paste("Sig in", cB),
        TRUE           ~ "Neither"
      ))

    # Sign consistency among union-significant SNP-probes
    union_sig  <- dat$sigA | dat$sigB
    n_union    <- sum(union_sig, na.rm = TRUE)
    if (n_union > 0) {
      sign_agree  <- sign(dat$x[union_sig]) == sign(dat$y[union_sig])
      consistency <- mean(sign_agree, na.rm = TRUE)
      consist_str <- sprintf("\n%.1f%% sign-consistent among union-sig (n=%d)", consistency * 100, n_union)
    } else {
      consist_str <- "\n(no union-significant pairs)"
    }

    pl <- ggplot(dat, aes(x = x, y = y, color = status)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c(
        "Both" = "#D62728", "Neither" = "#AAAAAA",
        setNames(c("#FF7F0E", "#1F77B4"), paste("Sig in", c(cA, cB)))
      ), drop = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      labs(title = paste0("Posterior mean: ", cA, " vs ", cB, consist_str),
           x = paste0("PM (", cA, ")"), y = paste0("PM (", cB, ")"),
           color = NULL) +
      theme_classic(base_size = 13) +
      theme(legend.position = "bottom")

    ggsave(file.path(out_dir, paste0(stem, ".png")), pl, width = 6, height = 5, dpi = 150)
    message("  Saved: ", stem, ".png")

    csv09_rows[[length(csv09_rows) + 1]] <- dat %>%
      mutate(pair = paste(cA, "vs", cB), classA = cA, classB = cB) %>%
      rename(PM_A = x, PM_B = y, sig_A = sigA, sig_B = sigB)
  }
  save_csv(do.call(rbind, csv09_rows), "09_pm_pairwise_scatter")
}

# ─────────────────────────────────────────────────────────────────────────────
# Plot 10 — Volcano plot per cell class
# ─────────────────────────────────────────────────────────────────────────────
message("[Plot 10] Volcano plots per cell class")

gene_sym <- NULL
if (!is.null(gene_loc_path) && file.exists(gene_loc_path)) {
  gene_loc <- read.table(gene_loc_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (all(c("ensg_id", "gene_symbol") %in% colnames(gene_loc))) {
    gene_sym <- gene_loc %>% select(ensg_id, gene_symbol) %>% distinct(ensg_id, .keep_all = TRUE)
    message("  Gene symbols loaded: ", nrow(gene_sym), " entries")
  }
} else if (!is.null(gene_loc_path)) {
  message("  WARNING: --gene-loc file not found: ", gene_loc_path)
}

vol_base <- data.frame(SNP_probe = snp_probes, stringsAsFactors = FALSE)
vol_base$ensg_id <- sub("^[^_]+_", "", snp_probes)
if (!is.null(gene_sym)) {
  vol_base <- left_join(vol_base, gene_sym, by = "ensg_id")
  vol_base$label <- ifelse(is.na(vol_base$gene_symbol) | vol_base$gene_symbol == "",
                           vol_base$ensg_id, vol_base$gene_symbol)
} else {
  vol_base$label <- vol_base$ensg_id
}

volcano_df <- do.call(rbind, lapply(cell_classes, function(cc) {
  data.frame(SNP_probe = snp_probes, cell_class = cc,
             PM = pm_mat[, cc], LFSR = lfsr_mat[, cc],
             ensg_id = vol_base$ensg_id, label = vol_base$label,
             stringsAsFactors = FALSE)
})) %>%
  mutate(neg_log10_lfsr = -log10(pmax(LFSR, 1e-10)), sig = LFSR < lfsr_thresh)

top_labels <- volcano_df %>% filter(sig) %>%
  group_by(cell_class) %>% slice_min(order_by = LFSR, n = 10) %>% ungroup()

p10 <- ggplot(volcano_df, aes(x = PM, y = neg_log10_lfsr, color = sig)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_text_repel(data = top_labels, aes(label = label),
                  size = 3, max.overlaps = 20, color = "#D62728") +
  scale_color_manual(values = c("FALSE" = "#AAAAAA", "TRUE" = "#D62728"),
                     labels = c("Not sig", "LFSR sig")) +
  facet_wrap(~cell_class, scales = "free") +
  labs(title = "Volcano plot per cell class",
       x = "Posterior mean effect", y = expression(-log[10](LFSR)), color = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

n_cc <- length(cell_classes)
save_gg(p10, "10_volcano_per_cellclass",
        width = min(3, n_cc) * 5, height = ceiling(n_cc / min(3, n_cc)) * 5)
save_csv(volcano_df %>% select(SNP_probe, cell_class, PM, LFSR, label, sig),
         "10_volcano_per_cellclass")

# ─────────────────────────────────────────────────────────────────────────────
# Plots 11 / 12 / 13 — LFSR density plots
#
# Uses the LFSR values from the strong-signal mashr model (lfsr_mat) merged with
# eqtl_present_in_all_sig_in_one.rds to obtain the original padj_snp values.
# Split by nominal significance (padj_snp) to show how well mashr recovers signal:
#   11 — nominally significant (padj_snp < nom_thresh)   → expect low LFSR
#   12 — nominally non-significant (padj_snp ≥ nom_thresh) → expect high LFSR
#   13 — random SNP-probes from mash_random               → should be flat / near 1
# ─────────────────────────────────────────────────────────────────────────────
if (!is.null(sig_in_one)) {
  message("[Plot 11/12] LFSR density: nom-sig vs nom-insig SNP-probes")

  # Long-format LFSR from the strong-signal model
  lfsr_sig_long <- lfsr_mat %>%
    as.data.frame() %>%
    rownames_to_column("SNP_probe") %>%
    pivot_longer(all_of(cell_classes), names_to = "cell_class", values_to = "lfsr")

  # Merge with sig_in_one to get padj_snp per (SNP_probe, cell_class)
  merged <- inner_join(lfsr_sig_long, sig_in_one[, c("SNP_probe", "cell_class", "padj_snp")],
                       by = c("SNP_probe", "cell_class"))
  message("  Merged rows: ", nrow(merged),
          " (lfsr_sig_long: ", nrow(lfsr_sig_long),
          ", sig_in_one: ", nrow(sig_in_one), ")")

  lfsr_nom_sig   <- merged[merged$padj_snp <  nom_thresh, ]
  lfsr_nom_insig <- merged[merged$padj_snp >= nom_thresh, ]
  message("  nom-sig rows: ", nrow(lfsr_nom_sig),
          "  nom-insig rows: ", nrow(lfsr_nom_insig))

  lfsr_theme <- theme(
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    plot.title   = element_text(size = 18),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 15)
  )

  if (nrow(lfsr_nom_sig) > 0) {
    p11 <- ggplot(lfsr_nom_sig, aes(x = lfsr, color = cell_class)) +
      geom_density(linewidth = 1) +
      scale_color_manual(values = color_vec) +
      labs(title = sprintf(
             "LFSR density — nominally significant SNP-probes\n(padj_snp < %.2f)", nom_thresh),
           x = "LFSR", y = "Density", color = "Cell class") +
      theme_classic(base_size = 14) + lfsr_theme
    save_gg(p11, "11_lfsr_density_nom_sig")
    save_csv(lfsr_nom_sig[, c("SNP_probe", "cell_class", "lfsr", "padj_snp")],
             "11_lfsr_density_nom_sig")
  } else {
    message("  No nom-sig rows — skipping Plot 11.")
  }

  if (nrow(lfsr_nom_insig) > 0) {
    p12 <- ggplot(lfsr_nom_insig, aes(x = lfsr, color = cell_class)) +
      geom_density(linewidth = 1) +
      scale_color_manual(values = color_vec) +
      labs(title = sprintf(
             "LFSR density — nominally non-significant SNP-probes\n(padj_snp ≥ %.2f)", nom_thresh),
           x = "LFSR", y = "Density", color = "Cell class") +
      theme_classic(base_size = 14) + lfsr_theme
    save_gg(p12, "12_lfsr_density_nom_insig")
    save_csv(lfsr_nom_insig[, c("SNP_probe", "cell_class", "lfsr", "padj_snp")],
             "12_lfsr_density_nom_insig")
  } else {
    message("  No nom-insig rows — skipping Plot 12.")
  }
} else {
  message("[Plot 11/12] Skipped — eqtl_present_in_all_sig_in_one.rds not available.")
}

if (!is.null(mash_random)) {
  message("[Plot 13] LFSR density: random SNP-probes")
  lfsr_rand_mat  <- get_lfsr(mash_random)
  rand_probes    <- rownames(lfsr_rand_mat)
  rand_classes   <- colnames(lfsr_rand_mat)

  lfsr_rand_long <- lfsr_rand_mat %>%
    as.data.frame() %>%
    rownames_to_column("SNP_probe") %>%
    pivot_longer(all_of(rand_classes), names_to = "cell_class", values_to = "lfsr")

  rand_color_vec <- color_vec[rand_classes]
  rand_color_vec[is.na(rand_color_vec)] <- hue_pal()(sum(is.na(rand_color_vec)))

  p13 <- ggplot(lfsr_rand_long, aes(x = lfsr, color = cell_class)) +
    geom_density(linewidth = 1) +
    scale_color_manual(values = rand_color_vec) +
    labs(title = "LFSR density — random SNP-probes",
         x = "LFSR", y = "Density", color = "Cell class") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x  = element_text(size = 13),
      axis.text.y  = element_text(size = 13),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title   = element_text(size = 18),
      legend.text  = element_text(size = 13),
      legend.title = element_text(size = 14)
    )
  save_gg(p13, "13_lfsr_density_random")
  save_csv(lfsr_rand_long, "13_lfsr_density_random")
} else {
  message("[Plot 13] Skipped — mashr random RDS not available.")
}

# ─────────────────────────────────────────────────────────────────────────────
message("\n[07b-plot-mashr] All done. Outputs in: ", out_dir)
