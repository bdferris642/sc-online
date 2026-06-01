#!/usr/bin/env Rscript
# run_depletion_analysis.R
#
# Correlates per-gene DE statistics (and per-pathway GSEA NES) with cell-type
# depletion estimates from crumblr, producing scatter plots, volcano plots,
# and heatmaps — then calls run-gsea.R on the gene r-values.
#
# Self-activates de-depletion-env — no manual env activation needed.
#
# Usage:
#   Rscript de-depletion-analysis/run_depletion_analysis.R \
#     --crumblr  /path/to/crumblr_results.csv \
#     --manifest /path/to/manifest.csv \
#     --output-dir /path/to/output/

# ══════════════════════════════════════════════════════════════════════════════
# Section 0 — Self-activation
# ══════════════════════════════════════════════════════════════════════════════
invisible(local({
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg  <- grep("^--file=", args_full, value = TRUE)
  if (!length(file_arg)) return()

  script      <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
  sandbox     <- normalizePath(file.path(dirname(script), ".."))
  env_rscript <- file.path(sandbox, "micromamba_root", "envs",
                           "de-depletion-env", "bin", "Rscript")
  env_r_home  <- file.path(sandbox, "micromamba_root", "envs",
                           "de-depletion-env", "lib", "R")

  if (file.exists(env_r_home) &&
      normalizePath(R.home()) == normalizePath(env_r_home)) return()

  if (!file.exists(env_rscript)) {
    message("ERROR: de-depletion-env not found at: ", env_rscript)
    message("Run: bash ", file.path(dirname(script), "setup.sh"))
    quit(status = 1, save = "no")
  }

  env_r_lib <- file.path(sandbox, "micromamba_root", "envs",
                         "de-depletion-env", "lib", "R", "library")
  Sys.unsetenv("R_HOME")
  Sys.setenv(R_LIBS = env_r_lib, R_LIBS_SITE = "", R_LIBS_USER = "")

  status <- system2(env_rscript,
                    c("--no-save", "--no-restore", "--no-site-file", "--no-init-file",
                      shQuote(script),
                      sapply(commandArgs(trailingOnly = TRUE), shQuote)))
  quit(status = status, save = "no")
}))

# ══════════════════════════════════════════════════════════════════════════════
# Section 1 — CLI Arguments
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--crumblr",       type = "character", default = NULL,
              help = "CSV with columns: cell_type, logFC, adj.P.Val, coefficient [required]"),
  make_option("--manifest",      type = "character", default = NULL,
              help = "CSV with columns: ct_name, de_csv_path, gsea_csv_path, de_rank_col, gsea_rank_col [required]"),
  make_option("--output-dir",    type = "character", default = NULL, dest = "output_dir",
              help = "Output directory (created if missing) [required]"),
  make_option("--min-de-sets",   type = "integer",   default = 5L,   dest = "min_de_sets",
              help = "Gene must appear in >= N DE files [default: %default]"),
  make_option("--top-n-genes",   type = "integer",   default = 50L,  dest = "top_n_genes",
              help = "Top genes by |r| in heatmap [default: %default]"),
  make_option("--top-n-pathways",type = "integer",   default = 50L,  dest = "top_n_pathways",
              help = "Top pathways by |r| in heatmap [default: %default]"),
  make_option("--focus-ct",      type = "character", default = NULL, dest = "focus_ct",
              help = "Partial string match for focus CT; default = most depleted CT"),
  make_option("--n-perm",        type = "integer",   default = 10000L, dest = "n_perm",
              help = "Permutations for empirical p-value [default: %default]"),
  make_option("--n-procs",       type = "integer",   default = 60L,  dest = "n_procs",
              help = "mc.cores for mclapply [default: %default]"),
  make_option("--force",         action = "store_true", default = FALSE,
              help = "Proceed with intersection on manifest/crumblr mismatch"),
  make_option("--hsp-genes",     type = "character",
              default = "/mnt/accessory/seq_data/pd-freeze/HSP-genes.txt", dest = "hsp_genes",
              help = "One gene per line [default: %default]"),
  make_option("--pd-gwas-genes", type = "character",
              default = "/mnt/accessory/seq_data/pd-freeze/GP2-pd-genes.txt", dest = "pd_gwas_genes",
              help = "One gene per line [default: %default]"),
  make_option("--n-label",       type = "integer",   default = 15L,  dest = "n_label",
              help = "Genes/pathways labelled per direction on scatter [default: %default]"),
  make_option("--padj-threshold",type = "double",    default = 0.05, dest = "padj_threshold",
              help = "FDR threshold for significance marking in heatmap borders [default: %default]"),
  make_option("--keep-bad-annotations", action = "store_true", default = FALSE,
              dest = "keep_bad_annotations",
              help = "Disable annotation filter (keep ENSG IDs, LINC, pseudogenes, etc.) [default: %default]"),
  make_option("--rank-padj-max",   type = "double",  default = Inf,   dest = "rank_padj_max",
              help = "NA-mask gene×CT where per-gene-CT padj >= this (derived col: padj_* from logFC_*); use 0.05 for marker analyses [default: %default]"),
  make_option("--intercept-floor", type = "double",  default = -15,   dest = "intercept_floor",
              help = "NA-mask gene×CT where logFC_(Intercept) < this; set to -Inf to disable [default: %default]"),
  make_option("--min-stat",        type = "double",  default = -10,   dest = "min_stat",
              help = "NA-mask gene×CT where rank stat < this (unstable extreme estimates) [default: %default]"),
  make_option("--max-stat",        type = "double",  default =  10,   dest = "max_stat",
              help = "NA-mask gene×CT where rank stat > this (unstable extreme estimates) [default: %default]"),
  make_option("--gsea-padj-max",   type = "double",  default = Inf,   dest = "gsea_padj_max",
              help = "NA-mask pathway×CT where GSEA padj >= this; disabled by default since most pathways are non-significant per CT [default: %default]")
)

parser <- OptionParser(
  option_list  = option_list,
  description  = "Correlate per-gene/pathway DE statistics with crumblr cell-type depletion"
)
args <- parse_args(parser)

for (req in c("crumblr", "manifest", "output_dir")) {
  if (is.null(args[[req]]))
    stop(sprintf("--%s is required. Run with --help for usage.",
                 gsub("_", "-", req)))
}
for (f in c(args$crumblr, args$manifest)) {
  if (!file.exists(f))
    stop("File not found: ", f)
}
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Section 2 — Libraries
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(stringr)
  library(scales)
  library(svglite)
  library(parallel)
})

# ══════════════════════════════════════════════════════════════════════════════
# Section 3 — Color palettes and annotation filters
# ══════════════════════════════════════════════════════════════════════════════

# Blue-GRAY-red (NOT blue-white-red): positive r = depletion = red
FILL_R   <- c("#2166ac", "#92c5de", "#d4d4d4", "#f4a582", "#b2182b")
# Purple-white-green: annotation strip (crumblr logFC)
FILL_LFC <- c("#762a83", "#f7f7f7", "#1b7837")

BAD_PATTERNS <- c(
  "^ENSG[0-9]",    # Ensembl IDs used as gene names
  "^LINC[0-9]",    # Long intergenic non-coding RNA
  "-AS[0-9]+$",    # Antisense transcripts
  "^LOC[0-9]",     # Uncharacterized loci
  "^MIR[0-9]",     # microRNA
  "^SNOR[A-Z]",    # Small nucleolar RNA
  "^RNU[0-9]",     # RNA U snoRNA
  "^AL[0-9]{6}",   # Uncharacterized AL... clones
  "^AC[0-9]{6}",   # Uncharacterized AC... clones
  "^RP[0-9]+-",    # Pseudogene RP11- etc.
  "\\.[0-9]+$",    # Versioned transcript names (gene.1)
  "-IT[0-9]+$",    # Intronic transcripts
  "-DT[0-9]?$",    # Divergent transcripts
  "^RNA5",         # Ribosomal RNA
  "^Z[0-9]{5}"     # Unknown Z-number genes
)
is_good_annotation <- function(genes)
  !Reduce(`|`, lapply(BAD_PATTERNS, function(p) grepl(p, genes, perl = TRUE)))

cap_symmetric <- function(x, q = 0.99) {
  lim <- quantile(abs(x), q, na.rm = TRUE)
  pmax(pmin(x, lim), -lim)
}

# ══════════════════════════════════════════════════════════════════════════════
# Section 4 — Logging
# ══════════════════════════════════════════════════════════════════════════════
log_con <- file(file.path(args$output_dir, "run.log"), open = "wt")
on.exit(close(log_con), add = TRUE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  message(msg)
  writeLines(msg, log_con)
  flush(log_con)
}

log_msg("Starting de-depletion analysis")
log_msg("=== Arguments ===")
for (arg_name in sort(names(args))) {
  log_msg(sprintf("  %-25s %s", arg_name, args[[arg_name]]))
}
log_msg("=== End Arguments ===")

# ══════════════════════════════════════════════════════════════════════════════
# Section 5 — Load crumblr results
# ══════════════════════════════════════════════════════════════════════════════
crumblr_df <- read.csv(args$crumblr, check.names = TRUE) |>
  filter(coefficient == "case_controlpd")

if (nrow(crumblr_df) == 0)
  stop("No rows with coefficient == 'case_controlpd' in ", args$crumblr)

for (col in c("cell_type", "logFC", "adj.P.Val")) {
  if (!col %in% colnames(crumblr_df))
    stop("Required column '", col, "' not found in crumblr CSV")
}

depletion_lfc <- setNames(crumblr_df$logFC, crumblr_df$cell_type)
log_msg(sprintf("Loaded %d cell types from crumblr", length(depletion_lfc)))

# ══════════════════════════════════════════════════════════════════════════════
# Section 6 — Load and validate manifest
# ══════════════════════════════════════════════════════════════════════════════
manifest <- read.csv(args$manifest, check.names = TRUE, stringsAsFactors = FALSE)
colnames(manifest) <- gsub("[.-]", "_", colnames(manifest))   # handle ct-name → ct_name

required_manifest_cols <- c("ct_name", "de_csv_path", "gsea_csv_path", "de_rank_col")
missing_cols <- setdiff(required_manifest_cols, colnames(manifest))
if (length(missing_cols) > 0)
  stop("Manifest missing required columns: ", paste(missing_cols, collapse = ", "))

# gsea_rank_col defaults to "NES" if column absent or NA
if (!"gsea_rank_col" %in% colnames(manifest)) {
  manifest$gsea_rank_col <- "NES"
} else {
  manifest$gsea_rank_col[is.na(manifest$gsea_rank_col) | !nzchar(manifest$gsea_rank_col)] <- "NES"
}

manifest_cts <- manifest$ct_name
crumblr_cts  <- crumblr_df$cell_type

missing_from_crumblr  <- setdiff(manifest_cts, crumblr_cts)
missing_from_manifest <- setdiff(crumblr_cts,  manifest_cts)

if (length(missing_from_crumblr) > 0 || length(missing_from_manifest) > 0) {
  if (length(missing_from_crumblr) > 0)
    message("In manifest, not in crumblr: ",
            paste(missing_from_crumblr, collapse = ", "))
  if (length(missing_from_manifest) > 0)
    message("In crumblr, not in manifest: ",
            paste(missing_from_manifest, collapse = ", "))
  if (!args$force)
    stop("ct_name / crumblr mismatch. Use --force to proceed with intersection.")
  keep          <- intersect(manifest_cts, crumblr_cts)
  manifest      <- filter(manifest, ct_name %in% keep)
  depletion_lfc <- depletion_lfc[keep]
  log_msg("--force: restricting to ", length(keep), " cell types in intersection")
}

# Validate all DE and GSEA paths exist
de_missing   <- manifest$de_csv_path[!file.exists(manifest$de_csv_path)]
gsea_missing <- manifest$gsea_csv_path[!file.exists(manifest$gsea_csv_path)]
if (length(de_missing) > 0)
  stop("DE CSV files not found:\n  ", paste(de_missing, collapse = "\n  "))
if (length(gsea_missing) > 0)
  stop("GSEA CSV files not found:\n  ", paste(gsea_missing, collapse = "\n  "))

ct_order <- manifest$ct_name
log_msg(sprintf("Manifest: %d cell types", length(ct_order)))
for (ct in ct_order)
  log_msg(sprintf("  %-30s logFC = %.3f", ct, depletion_lfc[ct]))

# ══════════════════════════════════════════════════════════════════════════════
# Section 7 — Load special gene sets
# ══════════════════════════════════════════════════════════════════════════════
load_gene_list <- function(path, label) {
  if (!file.exists(path)) {
    warning(label, " file not found: ", path)
    return(character(0))
  }
  genes <- trimws(readLines(path))
  genes <- genes[nzchar(genes)]
  log_msg(sprintf("%s: %d genes from %s", label, length(genes), path))
  if (length(genes) == 0)
    warning(label, " file has 0 genes: ", path)
  genes
}

hsp_genes    <- load_gene_list(args$hsp_genes,     "HSP genes")
pdgwas_genes <- load_gene_list(args$pd_gwas_genes, "PD GWAS genes")

# ══════════════════════════════════════════════════════════════════════════════
# Section 8 — Build stat matrix (genes × CTs)
# ══════════════════════════════════════════════════════════════════════════════
log_msg("Building gene stat matrix …")

gene_col_candidates <- c("gene_short_name", "gene", "gene_id")

load_de_vector <- function(row) {
  d <- read.csv(row$de_csv_path, check.names = FALSE, stringsAsFactors = FALSE)

  # Identify gene column
  gene_col <- NULL
  for (cand in gene_col_candidates) {
    if (cand %in% colnames(d)) { gene_col <- cand; break }
  }
  if (is.null(gene_col))
    stop(sprintf("No gene column (%s) found in: %s",
                 paste(gene_col_candidates, collapse = "/"), row$de_csv_path))

  # Support "col1+col2" syntax to sum two columns (e.g. intercept + marker logFC)
  rank_col_spec <- row$de_rank_col
  rank_parts    <- trimws(strsplit(rank_col_spec, "+", fixed = TRUE)[[1]])
  for (part in rank_parts) {
    if (!part %in% colnames(d))
      stop(sprintf("de_rank_col part '%s' not found in: %s", part, row$de_csv_path))
    if (!is.numeric(d[[part]]))
      stop(sprintf("de_rank_col part '%s' is not numeric in: %s", part, row$de_csv_path))
  }
  rank_vals <- Reduce("+", lapply(rank_parts, function(p) d[[p]]))

  # NA-mask: intercept expression floor
  intercept_col <- "logFC_(Intercept)"
  if (is.finite(args$intercept_floor) && intercept_col %in% colnames(d))
    rank_vals[d[[intercept_col]] < args$intercept_floor] <- NA_real_

  # NA-mask: extreme stat values (model instability)
  rank_vals[rank_vals < args$min_stat | rank_vals > args$max_stat] <- NA_real_

  setNames(rank_vals, d[[gene_col]])
}

load_padj_vector <- function(row) {
  d <- read.csv(row$de_csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- NULL
  for (cand in gene_col_candidates) {
    if (cand %in% colnames(d)) { gene_col <- cand; break }
  }
  rank_parts <- trimws(strsplit(row$de_rank_col, "+", fixed = TRUE)[[1]])
  # Derive padj col from first rank part (e.g. logFC_is_ct_* -> padj_is_ct_*)
  padj_col <- sub("^logFC_", "padj_", rank_parts[1])
  if (!padj_col %in% colnames(d)) return(setNames(rep(NA_real_, nrow(d)), d[[gene_col]]))
  setNames(d[[padj_col]], d[[gene_col]])
}

de_vectors <- lapply(seq_len(nrow(manifest)), function(i) load_de_vector(manifest[i, ]))
names(de_vectors) <- ct_order

# Gene universe: present in >= min_de_sets files and passes annotation filter
gene_counts  <- table(unlist(lapply(de_vectors, names)))
gene_universe <- names(gene_counts)[gene_counts >= args$min_de_sets]
log_msg(sprintf("Genes in >= %d DE files: %d", args$min_de_sets, length(gene_universe)))

if (!args$keep_bad_annotations) {
  gene_universe <- gene_universe[is_good_annotation(gene_universe)]
  log_msg(sprintf("After annotation filter: %d genes", length(gene_universe)))
} else {
  log_msg(sprintf("Annotation filter disabled (--keep-bad-annotations): retaining all %d genes",
          length(gene_universe)))
}

stat_mat <- vapply(ct_order, function(ct) {
  v <- de_vectors[[ct]]
  v[match(gene_universe, names(v))]
}, numeric(length(gene_universe)))
rownames(stat_mat) <- gene_universe
colnames(stat_mat) <- ct_order

# Build padj matrix and apply as post-hoc NA overlay on stat_mat
padj_vectors <- lapply(seq_len(nrow(manifest)), function(i) load_padj_vector(manifest[i, ]))
names(padj_vectors) <- ct_order

padj_mat <- vapply(ct_order, function(ct) {
  v <- padj_vectors[[ct]]
  v[match(gene_universe, names(v))]
}, numeric(length(gene_universe)))
rownames(padj_mat) <- gene_universe
colnames(padj_mat) <- ct_order

if (is.finite(args$rank_padj_max)) {
  mask <- is.na(padj_mat) | padj_mat >= args$rank_padj_max
  stat_mat[mask] <- NA_real_
  log_msg(sprintf("After padj < %.2f mask: %.1f%% of gene×CT entries retained",
                  args$rank_padj_max,
                  100 * mean(!is.na(stat_mat))))
}

# ══════════════════════════════════════════════════════════════════════════════
# Section 9 — Build pathway matrix (pathways × CTs)
# ══════════════════════════════════════════════════════════════════════════════
log_msg("Building pathway NES matrix …")

pw_col_candidates <- c("pathway", "Pathway", "name")

load_gsea_vector <- function(row) {
  d <- read.csv(row$gsea_csv_path, check.names = FALSE, stringsAsFactors = FALSE)

  # Identify pathway column
  pw_col <- NULL
  for (cand in pw_col_candidates) {
    if (cand %in% colnames(d)) { pw_col <- cand; break }
  }
  if (is.null(pw_col))
    stop(sprintf("No pathway column (%s) found in: %s",
                 paste(pw_col_candidates, collapse = "/"), row$gsea_csv_path))

  rank_col <- row$gsea_rank_col
  if (!rank_col %in% colnames(d))
    stop(sprintf("gsea_rank_col '%s' not found in: %s", rank_col, row$gsea_csv_path))
  if (!is.numeric(d[[rank_col]]))
    stop(sprintf("gsea_rank_col '%s' is not numeric in: %s", rank_col, row$gsea_csv_path))

  # Create unique pathway ID: pathway__gene_set
  if ("gene_set" %in% colnames(d)) {
    pw_ids <- paste0(d[[pw_col]], "__", d$gene_set)
    pathways  <- d[[pw_col]]
    gene_sets <- d$gene_set
  } else {
    pw_ids    <- d[[pw_col]]
    pathways  <- d[[pw_col]]
    gene_sets <- rep(NA_character_, nrow(d))
  }

  nes_vals <- d[[rank_col]]

  # NA-mask: GSEA padj filter
  if (is.finite(args$gsea_padj_max) && "padj" %in% colnames(d))
    nes_vals[is.na(d$padj) | d$padj >= args$gsea_padj_max] <- NA_real_

  list(
    values    = setNames(nes_vals, pw_ids),
    pathway   = setNames(pathways,  pw_ids),
    gene_set  = setNames(gene_sets, pw_ids)
  )
}

gsea_lists <- lapply(seq_len(nrow(manifest)), function(i) load_gsea_vector(manifest[i, ]))
names(gsea_lists) <- ct_order

# Pathway universe: present in >= min_de_sets GSEA files
pw_counts  <- table(unlist(lapply(gsea_lists, function(x) names(x$values))))
pw_universe <- names(pw_counts)[pw_counts >= args$min_de_sets]
log_msg(sprintf("Pathways in >= %d GSEA files: %d", args$min_de_sets, length(pw_universe)))

nes_mat <- vapply(ct_order, function(ct) {
  v <- gsea_lists[[ct]]$values
  v[match(pw_universe, names(v))]
}, numeric(length(pw_universe)))
rownames(nes_mat) <- pw_universe
colnames(nes_mat) <- ct_order

# Pathway metadata
first_pw  <- gsea_lists[[ct_order[1]]]
pw_meta <- data.frame(
  pw_id    = pw_universe,
  pathway  = {
    p <- first_pw$pathway[pw_universe]
    ifelse(is.na(p), sub("__.*$", "", pw_universe), p)
  },
  gene_set = {
    gs <- first_pw$gene_set[pw_universe]
    ifelse(is.na(gs), sub("^.*__", "", pw_universe), gs)
  },
  stringsAsFactors = FALSE
)

# ══════════════════════════════════════════════════════════════════════════════
# Section 10 — Pearson correlation (genes)
# ══════════════════════════════════════════════════════════════════════════════
n_ct            <- length(ct_order)
#effective_min_n <- max(3L, floor(n_ct * 2/3))
#log_msg(sprintf("effective_min_n = %d (floor(%d * 2/3))", effective_min_n, n_ct))
effective_min_n <- max(3L, floor(n_ct * 0.5))
log_msg(sprintf("effective_min_n = %d (floor(%d * 0.5))", effective_min_n, n_ct))

log_msg(sprintf("Running gene permutations: %d genes × %d perms, %d cores …",
                nrow(stat_mat), args$n_perm, args$n_procs))

gene_corr_list <- mclapply(seq_len(nrow(stat_mat)), function(gi) {
  x  <- stat_mat[gi, ]
  ok <- !is.na(x)
  n  <- sum(ok)
  if (n < effective_min_n) return(c(r = NA_real_, n_ct = n, emp_p = NA_real_))
  y      <- -depletion_lfc[names(x)[ok]]   # negate: positive r = vulnerability
  r_obs  <- cor(x[ok], y, method = "pearson")
  null_r <- vapply(seq_len(args$n_perm), function(i)
    cor(x[ok], sample(y), method = "pearson"), numeric(1))
  c(r = r_obs, n_ct = n, emp_p = mean(abs(null_r) >= abs(r_obs)))
}, mc.cores = args$n_procs)

gene_corr <- do.call(rbind, gene_corr_list)
rownames(gene_corr) <- gene_universe

gene_df <- data.frame(
  gene          = gene_universe,
  r             = gene_corr[, "r"],
  n_ct          = gene_corr[, "n_ct"],
  emp_p         = gene_corr[, "emp_p"],
  mean_abs_stat = rowMeans(abs(stat_mat), na.rm = TRUE),
  stringsAsFactors = FALSE
) |>
  filter(!is.na(r)) |>
  mutate(
    emp_p_floored = pmax(emp_p, 1 / args$n_perm),
    padj          = p.adjust(emp_p_floored, method = "BH"),
    signed_nlp    = sign(r) * -log10(padj)
  ) |>
  arrange(desc(signed_nlp), desc(abs(r)))

log_msg(sprintf("Genes with r: %d", nrow(gene_df)))

# ══════════════════════════════════════════════════════════════════════════════
# Section 11 — Pearson correlation (pathways)
# ══════════════════════════════════════════════════════════════════════════════
log_msg(sprintf("Running pathway permutations: %d pathways × %d perms, %d cores …",
                nrow(nes_mat), args$n_perm, args$n_procs))

pw_corr_list <- mclapply(seq_len(nrow(nes_mat)), function(pi) {
  x  <- nes_mat[pi, ]
  ok <- !is.na(x)
  n  <- sum(ok)
  if (n < effective_min_n) return(c(r = NA_real_, n_ct = n, emp_p = NA_real_))
  y      <- -depletion_lfc[names(x)[ok]]   # negate: positive r = vulnerability
  r_obs  <- cor(x[ok], y, method = "pearson")
  null_r <- vapply(seq_len(args$n_perm), function(i)
    cor(x[ok], sample(y), method = "pearson"), numeric(1))
  c(r = r_obs, n_ct = n, emp_p = mean(abs(null_r) >= abs(r_obs)))
}, mc.cores = args$n_procs)
pw_corr <- do.call(rbind, pw_corr_list)
rownames(pw_corr) <- pw_universe

pw_df <- pw_meta |>
  mutate(
    r           = pw_corr[pw_id, "r"],
    n_ct        = pw_corr[pw_id, "n_ct"],
    emp_p       = pw_corr[pw_id, "emp_p"],
    max_abs_nes = apply(abs(nes_mat), 1, max, na.rm = TRUE)
  ) |>
  filter(!is.na(r)) |>
  mutate(emp_p_floored = pmax(emp_p, 1 / args$n_perm)) |>
  group_by(gene_set) |>
  mutate(padj = p.adjust(emp_p_floored, method = "BH")) |>
  ungroup() |>
  arrange(desc(abs(r)))

log_msg(sprintf("Pathways with r: %d", nrow(pw_df)))

# ══════════════════════════════════════════════════════════════════════════════
# Section 12 — Focus CT detection
# ══════════════════════════════════════════════════════════════════════════════
if (is.null(args$focus_ct)) {
  focus_ct <- ct_order[which.min(depletion_lfc[ct_order])]
  log_msg(sprintf("Focus CT (most depleted): %s  logFC = %.3f",
                  focus_ct, depletion_lfc[focus_ct]))
} else {
  matches <- ct_order[grepl(args$focus_ct, ct_order, fixed = TRUE)]
  if (length(matches) == 0)
    stop("--focus-ct '", args$focus_ct, "' matched no cell types. Available: ",
         paste(ct_order, collapse = ", "))
  if (length(matches) > 1)
    stop("--focus-ct '", args$focus_ct, "' matched multiple cell types: ",
         paste(matches, collapse = ", "))
  focus_ct <- matches
  log_msg(sprintf("Focus CT (--focus-ct '%s'): %s", args$focus_ct, focus_ct))
}

focus_suffix <- paste0("_focus_", gsub("[^A-Za-z0-9_]", "_", focus_ct))

# ══════════════════════════════════════════════════════════════════════════════
# Section 13 — Plot helpers
# ══════════════════════════════════════════════════════════════════════════════

# ── Annotation strip (top of heatmap) ─────────────────────────────────────────
make_col_annotation_strip <- function(ct_ord, dep_lfc) {
  lfc_lim <- max(abs(dep_lfc[ct_ord]))
  df <- data.frame(ct    = factor(ct_ord, levels = ct_ord),
                   logFC = dep_lfc[ct_ord])
  ggplot(df, aes(x = ct, y = 1, fill = logFC)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_gradientn(
      colors = FILL_LFC, limits = c(-lfc_lim, lfc_lim),
      name   = "crumblr\nlogFC",
      guide  = guide_colorbar(title.position = "top", barwidth = 0.8, barheight = 3)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void(base_size = 9) +
    theme(legend.position = "right")
}

assemble_heatmap <- function(p_ann, p_main, n_rows, label = NULL) {
  leg_lfc  <- get_legend(p_ann  + theme(legend.position = "right"))
  leg_stat <- get_legend(p_main + theme(legend.position = "right"))

  p_ann_bare  <- p_ann  + theme(legend.position = "none")
  p_main_bare <- p_main + theme(legend.position = "none")

  aligned <- align_plots(p_ann_bare, p_main_bare, align = "v", axis = "lr")
  body    <- plot_grid(aligned[[1]], aligned[[2]], ncol = 1, rel_heights = c(2, n_rows))
  legends <- plot_grid(leg_lfc, NULL, leg_stat, ncol = 1, rel_heights = c(3, 1, n_rows))
  p_out   <- plot_grid(body, legends, ncol = 2, rel_widths = c(8, 1.5))

  if (!is.null(label) && nzchar(label)) {
    title_row <- ggdraw() + draw_label(label, size = 10, fontface = "bold", hjust = 0.5)
    p_out <- plot_grid(title_row, p_out, ncol = 1, rel_heights = c(0.04, 1))
  }
  p_out
}

save_svg_csv <- function(p, base_path, width, height, csv_df = NULL) {
  svglite(paste0(base_path, ".svg"), width = width, height = height)
  print(p)
  dev.off()
  if (!is.null(csv_df))
    write.csv(csv_df, paste0(base_path, ".csv"), row.names = FALSE)
}

# ── Gene scatter ───────────────────────────────────────────────────────────────
make_gene_scatter <- function(gene_df, n_ct, n_label, output_dir,
                              suffix = "", label = NULL) {
  if (nrow(gene_df) == 0) { log_msg(sprintf("[gene scatter] 0 rows, skipping%s", suffix)); return() }

  sig_genes   <- filter(gene_df, emp_p < 0.05)
  label_pool  <- if (nrow(sig_genes) > 0) sig_genes else gene_df
  label_genes <- unique(c(
    label_pool |> slice_max(r, n = n_label, with_ties = FALSE) |> pull(gene),
    label_pool |> slice_min(r, n = n_label, with_ties = FALSE) |> pull(gene)
  ))
  subtitle <- sprintf("%d genes · %d cell types · r > 0 = more DE in depleted types",
                      nrow(gene_df), n_ct)
  if (!is.null(label) && nzchar(label)) subtitle <- paste0(subtitle, "\n", label)

  p <- ggplot(gene_df, aes(x = r, y = mean_abs_stat)) +
    geom_point(aes(colour = r), alpha = 0.55, size = 0.9) +
    geom_text_repel(
      data = filter(gene_df, gene %in% label_genes), aes(label = gene),
      size = 2.4, max.overlaps = 40, segment.size = 0.3
    ) +
    scale_colour_gradientn(colors = FILL_R, limits = c(-1, 1), name = "Pearson r",
                           oob = squish) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    labs(
      x     = "Pearson r  (DE stat ~ crumblr logFC)",
      y     = "Mean |DE stat| across cell types",
      title = "Gene-level DE vs. depletion correlation",
      subtitle = subtitle
    ) +
    theme_bw(base_size = 11)

  base    <- file.path(output_dir, paste0("gene_depletion_correlation_scatter", suffix))
  csv_out <- gene_df |> rename(gene_short_name = gene)  # run-gsea.R compatibility
  save_svg_csv(p, base, 8, 6, csv_out)
  log_msg(sprintf("[gene scatter + csv]%s  (%d genes)", suffix, nrow(gene_df)))
}

# ── Gene volcano ───────────────────────────────────────────────────────────────
make_gene_volcano <- function(gene_df, n_label, output_dir,
                              suffix = "", label = NULL) {
  if (nrow(gene_df) == 0) { log_msg(sprintf("[gene volcano] 0 rows, skipping%s", suffix)); return() }

  sig_genes   <- filter(gene_df, emp_p < 0.05)
  label_pool  <- if (nrow(sig_genes) > 0) sig_genes else gene_df
  label_genes <- unique(c(
    label_pool |> slice_max(r,             n = n_label, with_ties = FALSE) |> pull(gene),
    label_pool |> slice_min(r,             n = n_label, with_ties = FALSE) |> pull(gene),
    label_pool |> slice_max(-log10(padj),  n = n_label, with_ties = FALSE) |> pull(gene)
  ))

  n_ct_max <- max(gene_df$n_ct, na.rm = TRUE)
  subtitle  <- sprintf("%d genes · up to %d cell types", nrow(gene_df), n_ct_max)
  if (!is.null(label) && nzchar(label)) subtitle <- paste0(subtitle, "\n", label)

  plot_df <- gene_df |>
    mutate(neg_log_padj = pmin(-log10(padj), 4),
           is_label     = gene %in% label_genes)

  p <- ggplot(plot_df, aes(x = r, y = neg_log_padj)) +
    geom_point(aes(colour = r), alpha = 0.45, size = 0.7) +
    geom_text_repel(
      data = filter(plot_df, is_label), aes(label = gene),
      size = 2.2, max.overlaps = 50, segment.size = 0.3
    ) +
    scale_colour_gradientn(colors = FILL_R, limits = c(-1, 1), name = "Pearson r",
                           oob = squish) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(
      x        = "Pearson r  (DE stat ~ crumblr logFC)",
      y        = expression(-log[10](p[adj])),
      title    = "Gene-level DE vs. depletion correlation",
      subtitle = subtitle
    ) +
    theme_bw(base_size = 11)

  base    <- file.path(output_dir, paste0("gene_depletion_volcano", suffix))
  csv_out <- gene_df |> select(gene, r, n_ct, emp_p, padj, signed_nlp, mean_abs_stat)
  save_svg_csv(p, base, 8, 6, csv_out)
  log_msg(sprintf("[gene volcano + csv]%s  (%d genes)", suffix, nrow(gene_df)))
}

# ── Gene heatmap ───────────────────────────────────────────────────────────────
make_gene_heatmap <- function(stat_mat, gene_df, dep_lfc, ct_ord,
                              top_n, padj_threshold, output_dir,
                              suffix = "", label = NULL) {
  ct_ord <- ct_ord[order(dep_lfc[ct_ord])]   # most depleted (most negative logFC) first
  avail <- intersect(gene_df$gene, rownames(stat_mat))
  gene_df <- filter(gene_df, gene %in% avail)
  if (nrow(gene_df) == 0) { log_msg(sprintf("[gene heatmap] 0 rows, skipping%s", suffix)); return() }

  top_genes <- gene_df |> slice_max(abs(r), n = top_n, with_ties = FALSE) |>
    arrange(r) |> pull(gene)
  hm     <- stat_mat[top_genes, ct_ord, drop = FALSE]
  hm_cap <- apply(hm, 2, cap_symmetric)
  val_lim <- max(abs(hm_cap), na.rm = TRUE)
  fsize_row <- max(4, min(7, 350 / length(top_genes)))

  df <- as.data.frame(hm_cap) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "cell_type", values_to = "stat") |>
    mutate(gene      = factor(gene,      levels = top_genes),
           cell_type = factor(cell_type, levels = ct_ord))

  # Significance borders: genes where padj < padj_threshold in each CT
  # We need padj per gene per CT — recompute from stat_mat subset via permutation
  # For simplicity: mark cells where |stat| is non-NA and gene has global padj < threshold
  sig_df <- gene_df |>
    filter(gene %in% top_genes) |>
    select(gene, padj) |>
    crossing(cell_type = ct_ord) |>
    mutate(gene      = factor(gene,      levels = top_genes),
           cell_type = factor(cell_type, levels = ct_ord),
           sig       = padj < padj_threshold)

  p_main <- ggplot(df, aes(x = cell_type, y = gene)) +
    geom_tile(aes(fill = stat), colour = "white", linewidth = 0.1) +
    geom_tile(data = filter(sig_df, sig),
              fill = NA, colour = "black", linewidth = 0.35) +
    scale_fill_gradientn(
      colors = FILL_R, limits = c(-val_lim, val_lim), na.value = "grey90",
      name   = "DE stat",
      guide  = guide_colorbar(title.position = "top", barwidth = 0.8, barheight = 6)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y  = element_text(size = fsize_row),
      axis.title   = element_blank(),
      panel.grid   = element_blank(),
      legend.position = "right"
    )

  p_ann <- make_col_annotation_strip(ct_ord, dep_lfc)
  p_out <- assemble_heatmap(p_ann, p_main, length(top_genes), label = label)

  w <- max(6, 2.5 + length(ct_ord) * 0.65)
  h <- max(5, 1.5 + length(top_genes) * 0.13) +
       if (!is.null(label) && nzchar(label)) 0.3 else 0

  base    <- file.path(output_dir, paste0("gene_depletion_heatmap", suffix))
  hm_csv  <- as.data.frame(hm) |>
    rownames_to_column("gene") |>
    left_join(select(gene_df, gene, r), by = "gene") |>
    select(gene, r, all_of(ct_ord))
  save_svg_csv(p_out, base, w, h, hm_csv)
  log_msg(sprintf("[gene heatmap + csv]%s  (%d genes)", suffix, length(top_genes)))
}

# ── Pathway scatter ────────────────────────────────────────────────────────────
make_pathway_scatter <- function(pw_df, n_ct, n_label, output_dir,
                                 suffix = "", label = NULL) {
  if (nrow(pw_df) == 0) { log_msg(sprintf("[pathway scatter] 0 rows, skipping%s", suffix)); return() }

  sig_pw     <- filter(pw_df, emp_p < 0.05)
  label_pool <- if (nrow(sig_pw) > 0) sig_pw else pw_df
  label_pw <- union(
    label_pool |> group_by(gene_set) |> slice_max(r, n = n_label, with_ties = FALSE) |> pull(pw_id),
    label_pool |> group_by(gene_set) |> slice_min(r, n = n_label, with_ties = FALSE) |> pull(pw_id)
  )
  pw_plot <- pw_df |>
    mutate(label_short = ifelse(nchar(pathway) > 55, paste0(substr(pathway, 1, 52), "..."), pathway))

  subtitle <- sprintf("%d pathways · %d cell types", nrow(pw_df), n_ct)
  if (!is.null(label) && nzchar(label)) subtitle <- paste0(subtitle, "\n", label)

  n_dbs   <- length(unique(pw_df$gene_set))
  p_width  <- min(22, 7 * min(n_dbs, 3))
  p_height <- 5 * ceiling(n_dbs / 3)

  p <- ggplot(pw_plot, aes(x = r, y = max_abs_nes)) +
    geom_point(aes(colour = r), alpha = 0.6, size = 1.1) +
    geom_text_repel(
      data = filter(pw_plot, pw_id %in% label_pw), aes(label = label_short),
      size = 2.0, max.overlaps = 25, segment.size = 0.3, min.segment.length = 0
    ) +
    scale_colour_gradientn(colors = FILL_R, limits = c(-1, 1), name = "Pearson r",
                           oob = squish) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    facet_wrap(~ gene_set, scales = "free_y", ncol = 3) +
    labs(
      x        = "Pearson r  (NES ~ crumblr logFC)",
      y        = "Max |NES| across cell types",
      title    = "Pathway-level enrichment vs. depletion correlation",
      subtitle = subtitle
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text       = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey95"))

  base    <- file.path(output_dir, paste0("pathway_depletion_correlation_scatter", suffix))
  csv_out <- select(pw_plot, -label_short)
  save_svg_csv(p, base, p_width, p_height, csv_out)
  log_msg(sprintf("[pathway scatter + csv]%s  (%d pathways)", suffix, nrow(pw_df)))
}

# ── Pathway volcano ────────────────────────────────────────────────────────────
make_pathway_volcano <- function(pw_df, n_ct, output_dir, suffix = "", label = NULL) {
  if (nrow(pw_df) == 0) { log_msg(sprintf("[pathway volcano] 0 rows, skipping%s", suffix)); return() }

  eps <- 1 / args$n_perm
  pw_plot <- pw_df |>
    mutate(
      neg_log10_padj = pmin(-log10(padj + eps), 4),
      label_short    = ifelse(nchar(pathway) > 55, paste0(substr(pathway, 1, 52), "..."), pathway)
    )

  sig_pw     <- filter(pw_df, emp_p < 0.05)
  label_pool <- if (nrow(sig_pw) > 0) sig_pw else pw_df
  label_pw <- label_pool |>
    group_by(gene_set) |>
    slice_max(abs(r), n = 10, with_ties = FALSE) |>
    pull(pw_id)

  subtitle <- sprintf("%d pathways · %d cell types", nrow(pw_df), n_ct)
  if (!is.null(label) && nzchar(label)) subtitle <- paste0(subtitle, "\n", label)

  n_dbs   <- length(unique(pw_df$gene_set))
  p_width  <- min(22, 7 * min(n_dbs, 3))
  p_height <- 5 * ceiling(n_dbs / 3)

  p <- ggplot(pw_plot, aes(x = r, y = neg_log10_padj)) +
    geom_point(aes(colour = r), alpha = 0.6, size = 1.1) +
    geom_text_repel(
      data = filter(pw_plot, pw_id %in% label_pw), aes(label = label_short),
      size = 2.0, max.overlaps = 25, segment.size = 0.3, min.segment.length = 0
    ) +
    scale_colour_gradientn(colors = FILL_R, limits = c(-1, 1), name = "Pearson r",
                           oob = squish) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    facet_wrap(~ gene_set, scales = "free_x", ncol = 3) +
    labs(
      x        = "Pearson r  (NES ~ crumblr logFC)",
      y        = expression(-log[10](p[adj])),
      title    = "Pathway-level enrichment vs. depletion: volcano",
      subtitle = subtitle
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text       = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey95"))

  base    <- file.path(output_dir, paste0("pathway_depletion_volcano", suffix))
  csv_out <- select(pw_plot, -label_short)
  save_svg_csv(p, base, p_width, p_height, csv_out)
  log_msg(sprintf("[pathway volcano + csv]%s  (%d pathways)", suffix, nrow(pw_df)))
}

# ── Pathway heatmap ────────────────────────────────────────────────────────────
make_pathway_heatmap <- function(nes_mat, pw_df, dep_lfc, ct_ord,
                                 top_n, output_dir, suffix = "", label = NULL) {
  ct_ord <- ct_ord[order(dep_lfc[ct_ord])]   # most depleted (most negative logFC) first
  avail <- intersect(pw_df$pw_id, rownames(nes_mat))
  pw_df <- filter(pw_df, pw_id %in% avail)
  if (nrow(pw_df) == 0) { log_msg(sprintf("[pathway heatmap] 0 rows, skipping%s", suffix)); return() }

  sel <- pw_df |>
    slice_max(abs(r), n = top_n, with_ties = FALSE) |>
    arrange(gene_set, r) |>
    mutate(
      label    = ifelse(nchar(pathway) > 50, paste0(substr(pathway, 1, 47), "..."), pathway),
      db_short = str_remove(gene_set, "_[0-9]{4}.*$")
    ) |>
    mutate(label = make.unique(label, sep = " #"))

  hm      <- nes_mat[sel$pw_id, ct_ord, drop = FALSE]
  nes_lim <- max(abs(hm), na.rm = TRUE)
  fsize_row <- max(4, min(7, 500 / nrow(sel)))

  df <- as.data.frame(hm) |>
    rownames_to_column("pw_id") |>
    pivot_longer(-pw_id, names_to = "cell_type", values_to = "nes") |>
    left_join(select(sel, pw_id, label, db_short), by = "pw_id") |>
    mutate(
      label     = factor(label,     levels = sel$label),
      cell_type = factor(cell_type, levels = ct_ord),
      db_short  = factor(db_short,  levels = unique(sel$db_short))
    )

  p_main <- ggplot(df, aes(x = cell_type, y = label)) +
    geom_tile(aes(fill = nes), colour = "white", linewidth = 0.1) +
    scale_fill_gradientn(
      colors = FILL_R, limits = c(-nes_lim, nes_lim), na.value = "grey90",
      name   = "NES",
      guide  = guide_colorbar(title.position = "top", barwidth = 0.8, barheight = 6)
    ) +
    facet_grid(db_short ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y       = element_text(size = fsize_row),
      axis.title        = element_blank(),
      panel.grid        = element_blank(),
      legend.position   = "right",
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7, face = "bold"),
      strip.background  = element_rect(fill = "grey95"),
      panel.spacing     = unit(0.3, "lines")
    )

  p_ann <- make_col_annotation_strip(ct_ord, dep_lfc)
  p_out <- assemble_heatmap(p_ann, p_main, nrow(sel), label = label)

  w <- max(6, 2.5 + length(ct_ord) * 0.7)
  h <- max(5, 2.5 + nrow(sel) * 0.2) +
       if (!is.null(label) && nzchar(label)) 0.3 else 0

  base    <- file.path(output_dir, paste0("pathway_depletion_heatmap", suffix))
  hm_csv  <- as.data.frame(hm) |>
    rownames_to_column("pw_id") |>
    left_join(select(sel, pw_id, pathway, gene_set, r), by = "pw_id") |>
    select(pathway, gene_set, r, all_of(ct_ord))
  save_svg_csv(p_out, base, w, h, hm_csv)
  log_msg(sprintf("[pathway heatmap + csv]%s  (%d pathways)", suffix, nrow(sel)))
}

# ── Emit helpers ───────────────────────────────────────────────────────────────
emit_gene_outputs <- function(g_df, s_mat, suffix, label = NULL) {
  make_gene_scatter(g_df, n_ct, args$n_label, args$output_dir, suffix, label)
  make_gene_volcano(g_df, args$n_label,       args$output_dir, suffix, label)
  make_gene_heatmap(s_mat, g_df, depletion_lfc, ct_order,
                    args$top_n_genes, args$padj_threshold,
                    args$output_dir, suffix, label)
}

emit_pathway_outputs <- function(p_df, n_mat, suffix, label = NULL) {
  make_pathway_scatter(p_df, n_ct, args$n_label, args$output_dir, suffix, label)
  make_pathway_volcano(p_df, n_ct,               args$output_dir, suffix, label)
  make_pathway_heatmap(n_mat, p_df, depletion_lfc, ct_order,
                       args$top_n_pathways, args$output_dir, suffix, label)
}

# ══════════════════════════════════════════════════════════════════════════════
# Section 13 — Output generation
# ══════════════════════════════════════════════════════════════════════════════

# Helper to subset gene_df to a gene set
subset_gene_df <- function(g_df, gene_set) {
  filter(g_df, gene %in% gene_set)
}

# ── Union set ──────────────────────────────────────────────────────────────────
log_msg("=== Generating union set outputs ===")

hsp_df  <- subset_gene_df(gene_df, hsp_genes)
gwas_df <- subset_gene_df(gene_df, pdgwas_genes)

invisible(emit_gene_outputs(gene_df,  stat_mat, suffix = "",        label = "All genes"))
invisible(emit_gene_outputs(hsp_df,   stat_mat, suffix = "_hsp",    label = "HSP genes"))
invisible(emit_gene_outputs(gwas_df,  stat_mat, suffix = "_pdgwas", label = "PD GWAS genes"))
invisible(emit_pathway_outputs(pw_df, nes_mat,  suffix = ""))

# ── Focus set ──────────────────────────────────────────────────────────────────
log_msg(sprintf("=== Generating focus set outputs (focus CT: %s) ===", focus_ct))

focus_label <- sprintf("Focus: %s  (logFC = %.3f)", focus_ct, depletion_lfc[focus_ct])

# Focus gene set: genes significant in focus CT (using global padj on gene_df)
# Since we compute global padj (not per-CT), we use the DE stat for focus CT
focus_genes <- gene_df |>
  filter(padj < args$padj_threshold) |>
  pull(gene)
log_msg(sprintf("Focus gene set: %d genes with padj < %.2f", length(focus_genes), args$padj_threshold))

focus_gene_df  <- subset_gene_df(gene_df,  focus_genes)
focus_hsp_df   <- subset_gene_df(hsp_df,   focus_genes)
focus_gwas_df  <- subset_gene_df(gwas_df,  focus_genes)

invisible(emit_gene_outputs(focus_gene_df,  stat_mat, suffix = focus_suffix,
                            label = focus_label))
invisible(emit_gene_outputs(focus_hsp_df,   stat_mat, suffix = paste0(focus_suffix, "_hsp"),
                            label = paste0(focus_label, "\nHSP genes")))
invisible(emit_gene_outputs(focus_gwas_df,  stat_mat, suffix = paste0(focus_suffix, "_pdgwas"),
                            label = paste0(focus_label, "\nPD GWAS genes")))

# Focus pathway set
focus_pw_df <- pw_df |> filter(padj < args$padj_threshold)
log_msg(sprintf("Focus pathway set: %d pathways with padj < %.2f",
                nrow(focus_pw_df), args$padj_threshold))
invisible(emit_pathway_outputs(focus_pw_df, nes_mat, suffix = focus_suffix, label = focus_label))

# ══════════════════════════════════════════════════════════════════════════════
# Section 14 — run-gsea.R call
# ══════════════════════════════════════════════════════════════════════════════
args_full  <- commandArgs(trailingOnly = FALSE)
file_arg   <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg)) {
  script_path  <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
  sandbox_path <- normalizePath(file.path(dirname(script_path), ".."))
} else {
  sandbox_path <- normalizePath(file.path(getwd()))
}

run_gsea_script <- file.path(sandbox_path, "scripts", "run-gsea.R")
if (!file.exists(run_gsea_script)) {
  log_msg("WARNING: run-gsea.R not found at: ", run_gsea_script, " — skipping GSEA step")
} else {
  scatter_csv <- file.path(args$output_dir, "gene_depletion_correlation_scatter.csv")
  if (!file.exists(scatter_csv)) {
    log_msg("WARNING: scatter CSV not found for GSEA step: ", scatter_csv)
  } else {
    log_msg("Running run-gsea.R on gene r-values …")
    # Use Rscript directly since run-gsea.R may not be executable
    status <- system2("Rscript",
                      args = c("--no-save", "--no-restore",
                               shQuote(run_gsea_script),
                               "--path", shQuote(scatter_csv),
                               "--rank-col", "r"))
    if (status != 0)
      log_msg(sprintf("WARNING: run-gsea.R exited with status %d", status))
    else
      log_msg("run-gsea.R completed successfully")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Done
# ══════════════════════════════════════════════════════════════════════════════
log_msg(sprintf("All outputs written to: %s", normalizePath(args$output_dir)))
log_msg("Done.")
