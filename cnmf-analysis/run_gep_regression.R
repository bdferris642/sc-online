#!/usr/bin/env Rscript
# run_gep_regression.R
#
# For each row in a manifest TSV (Seurat .qs path + feature-loadings path):
#   1. Joins GEP scores onto cell metadata
#   2. Pseudo-bulks to participant level (mean aggregation)
#   3. Fits NNLS regression for raw GEP type; beta regression for normalized types
#   4. Collinearity detection (VIF + GEP~covariate R²)
#   5. Spearman correlation matrices and heatmaps (one per normalization type)
#   6. Sample heatmaps (absolute + z-scored) for id_col + multi-level categoricals
#   7. Full plot suite: volcano, distributions, coef heatmap, R² bar, PCA biplot
#   8. Writes gep_regression_results.csv and all plots per manifest row
#
# Self-activates varpart-venv — no manual env activation needed.
#
# Usage:
#   Rscript gep-regression/run_gep_regression.R \
#     --input gep-regression/test_manifest.tsv \
#     --output-dir gep-regression/test_output/

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
                           "cnmf-analysis-env", "bin", "Rscript")
  env_r_home  <- file.path(sandbox, "micromamba_root", "envs",
                           "cnmf-analysis-env", "lib", "R")

  if (file.exists(env_r_home) &&
      normalizePath(R.home()) == normalizePath(env_r_home)) return()

  if (!file.exists(env_rscript)) {
    message("ERROR: cnmf-analysis-env not found at: ", env_rscript)
    message("See cnmf-analysis/README.md for setup instructions.")
    quit(status = 1, save = "no")
  }

  env_r_lib <- file.path(sandbox, "micromamba_root", "envs",
                         "cnmf-analysis-env", "lib", "R", "library")
  Sys.unsetenv("R_HOME")
  Sys.setenv(R_LIBS = env_r_lib, R_LIBS_SITE = "", R_LIBS_USER = "")

  status <- system2(env_rscript,
                    c("--no-save", "--no-restore", "--no-site-file", "--no-init-file",
                      shQuote(script),
                      sapply(commandArgs(trailingOnly = TRUE), shQuote)))
  quit(status = status, save = "no")
}))

# ══════════════════════════════════════════════════════════════════════════════
# Section 1 — Argument parsing
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--input",      type = "character", default = NULL,
              help = "Path to manifest TSV [required]"),
  make_option("--output-dir", type = "character", default = NULL,
              help = "Base output directory [required]"),
  make_option("--ncores",     type = "integer",   default = NULL,
              help = "Parallel cores (default: detectCores()-1, capped at nrow(manifest))")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input))           stop("ERROR: --input is required")
if (is.null(opt[["output-dir"]])) stop("ERROR: --output-dir is required")

opt$output_dir <- opt[["output-dir"]]

if (!file.exists(opt$input)) stop(sprintf("ERROR: manifest not found: %s", opt$input))

manifest <- tryCatch(
  read.table(opt$input, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             quote = "", comment.char = ""),
  error = function(e) stop(sprintf("ERROR reading manifest '%s': %s",
                                   opt$input, conditionMessage(e)))
)

required_manifest_cols <- c("qs_path", "feature_loadings_path", "qc_geps", "id_col", "covs")
missing_mcols <- setdiff(required_manifest_cols, colnames(manifest))
if (length(missing_mcols) > 0) {
  stop(sprintf("ERROR: manifest missing columns: %s", paste(missing_mcols, collapse = ", ")))
}

n_rows <- nrow(manifest)
cat(sprintf("Manifest: %d row(s)\n", n_rows))

n_available <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
n_cores <- if (is.null(opt$ncores)) n_available else opt$ncores
n_cores <- min(n_cores, n_rows)
cat(sprintf("Using %d parallel core(s)\n", n_cores))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Section 2 — Load libraries
# ══════════════════════════════════════════════════════════════════════════════
required_pkgs <- c("SeuratObject", "dplyr", "betareg", "nnls",
                   "ggplot2", "svglite", "tidyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("ERROR: Required package '%s' not installed. Run setup_varpart_env.sh.", pkg))
  }
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(dplyr)
  library(betareg)
  library(nnls)
  library(ggplot2)
  library(svglite)
  library(tidyr)
})

# ══════════════════════════════════════════════════════════════════════════════
# Section 3 — Top-level helper functions
# ══════════════════════════════════════════════════════════════════════════════

# ── Utility: rownames to a column (avoids tibble dependency) ─────────────────
rn_to_col <- function(df, col = "rowname") {
  df2 <- df
  df2[[col]] <- rownames(df2)
  df2[, c(col, setdiff(colnames(df2), col)), drop = FALSE]
}

# ── Utility: save ggplot as SVG + PNG ────────────────────────────────────────
save_plot <- function(p, base_path, width_px, height_px, dpi = 150) {
  width_in  <- width_px  / dpi
  height_in <- height_px / dpi
  tryCatch({
    svglite::svglite(paste0(base_path, ".svg"), width = width_in, height = height_in)
    print(p)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    warning(sprintf("SVG save failed for %s: %s", base_path, conditionMessage(e)))
  })
  tryCatch({
    png(paste0(base_path, ".png"), width = width_px, height = height_px,
        res = dpi, units = "px")
    print(p)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    warning(sprintf("PNG save failed for %s: %s", base_path, conditionMessage(e)))
  })
}

# ── Utility: encode categorical covariates as numeric ────────────────────────
encode_covs <- function(df, covs) {
  df_enc <- df
  for (cv in covs) {
    if (!cv %in% colnames(df_enc)) next
    if (!is.numeric(df_enc[[cv]])) {
      df_enc[[cv]] <- as.numeric(as.factor(df_enc[[cv]]))
    }
  }
  df_enc
}

# ── Utility: majority vote for categorical aggregation ───────────────────────
majority_vote <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab        <- sort(table(x), decreasing = TRUE)
  candidates <- names(tab)[tab == tab[1]]
  sort(candidates)[1]
}

# ── NNLS regression for raw GEP type ─────────────────────────────────────────
run_nnls_regressions <- function(df, covs, gep_cols_new, row_idx) {
  n <- nrow(df)
  if (n < 3) {
    warning(sprintf("[row %d] gep_type='raw': only %d participants; skipping NNLS",
                    row_idx, n))
    return(NULL)
  }

  df_enc <- encode_covs(df, covs)

  A <- tryCatch(
    cbind(1, as.matrix(df_enc[, covs, drop = FALSE])),
    error = function(e) {
      warning(sprintf("[row %d] NNLS: failed to build design matrix: %s",
                      row_idx, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(A)) return(NULL)

  results_list <- vector("list", length(gep_cols_new))

  for (gi in seq_along(gep_cols_new)) {
    gep_name <- gep_cols_new[gi]
    if (!gep_name %in% colnames(df_enc)) next

    y            <- df_enc[[gep_name]]
    complete_idx <- which(!is.na(y) & rowSums(is.na(A)) == 0)
    if (length(complete_idx) < 3) {
      warning(sprintf("[row %d] raw %s: <3 complete cases; skipping", row_idx, gep_name))
      next
    }

    A_sub <- A[complete_idx, , drop = FALSE]
    y_sub <- y[complete_idx]

    fit <- tryCatch(
      nnls::nnls(A = A_sub, b = y_sub),
      error = function(e) {
        warning(sprintf("[row %d] raw %s NNLS failed: %s", row_idx, gep_name,
                        conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) next

    coefs     <- fit$x
    total_ss  <- sum((y_sub - mean(y_sub))^2)
    r_squared <- if (total_ss > 0) 1 - sum(fit$residuals^2) / total_ss else NA_real_

    results_list[[gi]] <- data.frame(
      gep_name       = gep_name,
      gep_type       = "raw",
      term           = c("(Intercept)", covs),
      estimate       = coefs,
      std_error      = NA_real_,
      z_value        = NA_real_,
      p_value        = NA_real_,
      pseudo_r2      = r_squared,
      n_participants = length(complete_idx),
      row.names      = NULL,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}

# ── Beta regression for normalized GEP types ─────────────────────────────────
run_beta_regressions <- function(df, gep_type, covs, gep_cols_new, row_idx) {
  n <- nrow(df)
  if (n < 3) {
    warning(sprintf("[row %d] gep_type='%s': only %d participants; skipping beta regression",
                    row_idx, gep_type, n))
    return(NULL)
  }

  df_enc      <- encode_covs(df, covs)
  formula_rhs <- paste(covs, collapse = " + ")
  results_list <- vector("list", length(gep_cols_new))

  for (gi in seq_along(gep_cols_new)) {
    gep_name <- gep_cols_new[gi]
    if (!gep_name %in% colnames(df_enc)) next

    y_raw <- df_enc[[gep_name]]
    if (all(is.na(y_raw))) {
      warning(sprintf("[row %d] %s %s: all NA; skipping", row_idx, gep_type, gep_name))
      next
    }

    # Smithson & Verkuilen boundary squeeze: y_adj = (y*(n-1)+0.5)/n
    y_adj <- (y_raw * (n - 1) + 0.5) / n

    if (any(y_adj <= 0 | y_adj >= 1, na.rm = TRUE)) {
      warning(sprintf("[row %d] %s %s: y_adj out of (0,1); clamping",
                      row_idx, gep_type, gep_name))
      y_adj <- pmax(pmin(y_adj, 1 - .Machine$double.eps * 100),
                    .Machine$double.eps * 100)
    }

    df_fit       <- df_enc[, covs, drop = FALSE]
    df_fit$y_adj <- y_adj
    df_fit       <- df_fit[complete.cases(df_fit), ]
    if (nrow(df_fit) < 3) {
      warning(sprintf("[row %d] %s %s: <3 complete cases; skipping", row_idx, gep_type, gep_name))
      next
    }

    frm <- as.formula(paste("y_adj ~", formula_rhs))

    fit <- tryCatch(
      betareg(frm, data = df_fit),
      error = function(e) {
        warning(sprintf("[row %d] %s %s betareg failed: %s", row_idx, gep_type, gep_name,
                        conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) next

    smry     <- summary(fit)
    coef_tbl <- as.data.frame(smry$coefficients$mean)

    results_list[[gi]] <- data.frame(
      gep_name       = gep_name,
      gep_type       = gep_type,
      term           = rownames(coef_tbl),
      estimate       = coef_tbl[, "Estimate"],
      std_error      = coef_tbl[, "Std. Error"],
      z_value        = coef_tbl[, "z value"],
      p_value        = coef_tbl[, "Pr(>|z|)"],
      pseudo_r2      = smry$pseudo.r.squared,
      n_participants = nrow(df_fit),
      row.names      = NULL,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}

# ── Collinearity detection ────────────────────────────────────────────────────
check_collinearity <- function(df, covs, gep_cols_new, context_label) {
  if (length(covs) < 2) return(invisible(NULL))

  df_enc <- encode_covs(df, covs)

  # 1. Covariate VIF
  for (j in seq_along(covs)) {
    cov_j      <- covs[j]
    other_covs <- setdiff(covs, cov_j)
    cols_need  <- intersect(c(cov_j, other_covs), colnames(df_enc))
    if (length(cols_need) < 2) next

    df_sub <- df_enc[complete.cases(df_enc[, cols_need, drop = FALSE]), cols_need, drop = FALSE]
    if (nrow(df_sub) <= length(other_covs)) next

    frm  <- as.formula(paste(cov_j, "~", paste(other_covs, collapse = " + ")))
    lm_j <- tryCatch(lm(frm, data = df_sub), error = function(e) NULL)
    if (is.null(lm_j)) next

    r2    <- summary(lm_j)$r.squared
    if (is.na(r2)) next
    vif_j <- if (r2 < 1) 1 / (1 - r2) else Inf
    if (vif_j > 10) {
      cat(sprintf("  [collinearity] %s: %s VIF=%.1f (threshold 10)\n",
                  context_label, cov_j, vif_j))
    }
  }

  # 2. GEP ~ covariates R²
  avail_geps <- intersect(gep_cols_new, colnames(df_enc))
  for (gep_name in avail_geps) {
    cols_need <- intersect(c(gep_name, covs), colnames(df_enc))
    if (length(cols_need) < 2) next

    df_sub <- df_enc[complete.cases(df_enc[, cols_need, drop = FALSE]), cols_need, drop = FALSE]
    if (nrow(df_sub) <= length(covs)) next

    frm  <- as.formula(paste(gep_name, "~", paste(covs, collapse = " + ")))
    lm_g <- tryCatch(lm(frm, data = df_sub), error = function(e) NULL)
    if (is.null(lm_g)) next

    r2 <- summary(lm_g)$r.squared
    if (!is.na(r2) && r2 > 0.9) {
      cat(sprintf("  [collinearity] %s: %s R\u00b2=%.3f against covariates \u2014 may be covariate-driven\n",
                  context_label, gep_name, r2))
    }
  }

  invisible(NULL)
}

# ── Identify eligible covariates for Spearman correlation ─────────────────────
# Returns list: $numeric_covs, $binary_covs, $eligible_names
eligible_spearman_covs <- function(df, covs) {
  numeric_covs <- character(0)
  binary_covs  <- character(0)

  for (cv in covs) {
    if (!cv %in% colnames(df)) next
    vals     <- df[[cv]]
    vals_nna <- vals[!is.na(vals)]
    n_unique <- length(unique(vals_nna))

    if (is.numeric(vals)) {
      numeric_covs <- c(numeric_covs, cv)
    } else if (n_unique == 2) {
      binary_covs <- c(binary_covs, cv)
    }
    # > 2 levels: excluded from Spearman
  }

  list(numeric_covs  = numeric_covs,
       binary_covs   = binary_covs,
       eligible_names = c(numeric_covs, binary_covs))
}

# ── Spearman correlation CSVs + heatmap ───────────────────────────────────────
run_spearman_analysis <- function(df_x, gep_cols_new, covs, type_label, out_subdir, row_idx) {
  elig <- eligible_spearman_covs(df_x, covs)

  gep_vars <- intersect(gep_cols_new, colnames(df_x))
  var_names <- c(gep_vars, elig$eligible_names)
  if (length(var_names) < 2) {
    cat(sprintf("  Spearman %s: <2 eligible variables; skipping\n", type_label))
    return(invisible(NULL))
  }

  # Encode binary covariates as 0/1
  df_plot <- df_x
  for (cv in elig$binary_covs) {
    lv <- sort(unique(df_plot[[cv]][!is.na(df_plot[[cv]])]))
    df_plot[[cv]] <- as.integer(df_plot[[cv]] == lv[2])
  }

  var_mat <- as.matrix(df_plot[, var_names, drop = FALSE])
  storage.mode(var_mat) <- "numeric"

  n_vars <- length(var_names)
  n_geps <- length(gep_vars)
  n_cov  <- length(elig$eligible_names)

  # Spearman rho matrix
  rho_mat <- cor(var_mat, method = "spearman", use = "pairwise.complete.obs")

  # P-value matrix (upper triangle, mirrored)
  pval_mat <- matrix(NA_real_, nrow = n_vars, ncol = n_vars,
                     dimnames = list(var_names, var_names))
  diag(pval_mat) <- 1
  for (i in seq_len(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      x  <- var_mat[, i]; y <- var_mat[, j]
      ok <- !is.na(x) & !is.na(y)
      if (sum(ok) < 5) next
      res <- tryCatch(cor.test(x[ok], y[ok], method = "spearman"),
                      error = function(e) NULL)
      if (!is.null(res)) {
        pval_mat[i, j] <- res$p.value
        pval_mat[j, i] <- res$p.value
      }
    }
  }

  # Write CSVs
  rho_df  <- rn_to_col(as.data.frame(rho_mat),  "variable")
  pval_df <- rn_to_col(as.data.frame(pval_mat), "variable")
  write.csv(rho_df,  file.path(out_subdir, sprintf("spearman_corr_%s.csv",  type_label)),
            row.names = FALSE)
  write.csv(pval_df, file.path(out_subdir, sprintf("spearman_pval_%s.csv",  type_label)),
            row.names = FALSE)
  cat(sprintf("  Saved: spearman_corr_%s.csv + spearman_pval_%s.csv\n", type_label, type_label))

  # Heatmap
  tryCatch({
    rho_long <- tidyr::pivot_longer(
      rn_to_col(as.data.frame(rho_mat), "row_var"),
      cols = -"row_var", names_to = "col_var", values_to = "rho"
    )
    pval_long <- tidyr::pivot_longer(
      rn_to_col(as.data.frame(pval_mat), "row_var"),
      cols = -"row_var", names_to = "col_var", values_to = "pval"
    )
    merged <- merge(rho_long, pval_long, by = c("row_var", "col_var"))

    # Factor ordering: GEPs at left/top, covariates at right/bottom
    cov_vars  <- elig$eligible_names
    var_order <- c(gep_vars, cov_vars)

    # col_var: GEPs at positions 1..n_geps (left), covariates at n_geps+1..n_vars (right)
    merged$col_var <- factor(merged$col_var, levels = var_order)
    # row_var: reversed so GEPs appear at top, covariates at bottom
    # rev(var_order) = c(rev(cov_vars), rev(gep_vars))
    # → bottom levels = last_cov..first_cov at y=1..n_cov
    # → top levels   = last_gep..first_gep at y=n_cov+1..n_vars
    merged$row_var <- factor(merged$row_var, levels = rev(var_order))

    # Significance stars in GEP×covariate cells
    row_is_gep <- as.character(merged$row_var) %in% gep_vars
    col_is_gep <- as.character(merged$col_var) %in% gep_vars
    merged$is_gep_cov <- (row_is_gep & !col_is_gep) | (!row_is_gep & col_is_gep)

    merged$stars <- ifelse(!is.na(merged$pval) & merged$pval < 0.001, "***",
                    ifelse(!is.na(merged$pval) & merged$pval < 0.01,  "**",
                    ifelse(!is.na(merged$pval) & merged$pval < 0.05,  "*",  "")))
    merged$label <- ifelse(merged$is_gep_cov, merged$stars, "")

    # Suppress text in GEP×GEP block when n_geps > 15
    if (n_geps > 15) {
      merged$label[row_is_gep & col_is_gep] <- ""
    }

    # Separation line positions (numeric on discrete axis: level i → position i)
    sep_x <- n_geps + 0.5   # between GEPs and covariates on x-axis
    sep_y <- n_cov  + 0.5   # between covariates (bottom) and GEPs (top) on y-axis

    p <- ggplot(merged, aes(x = col_var, y = row_var, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.2) +
      geom_text(aes(label = label), size = 2.5, color = "black") +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                           midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho",
                           na.value = "grey80") +
      labs(title = sprintf("Spearman correlation \u2014 %s", type_label),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            panel.grid   = element_blank())

    if (n_cov > 0 && n_geps > 0) {
      p <- p +
        annotate("segment", x = sep_x, xend = sep_x,
                 y = 0.5, yend = n_vars + 0.5, color = "black", linewidth = 0.8) +
        annotate("segment", x = 0.5, xend = n_vars + 0.5,
                 y = sep_y, yend = sep_y, color = "black", linewidth = 0.8)
    }

    w_px <- max(600L, as.integer(n_vars * 40))
    h_px <- max(600L, as.integer(n_vars * 40))
    save_plot(p, file.path(out_subdir, sprintf("spearman_heatmap_%s", type_label)),
              w_px, h_px)
    cat(sprintf("  Saved: spearman_heatmap_%s.svg + .png\n", type_label))
  }, error = function(e) {
    warning(sprintf("[row %d] Spearman heatmap %s failed: %s",
                    row_idx, type_label, conditionMessage(e)))
  })

  invisible(NULL)
}

# ── Sample heatmaps (absolute + z-scored) ─────────────────────────────────────
run_sample_heatmaps <- function(df_x, target_var, gep_cols_new, type_label, out_subdir,
                                row_idx, max_label_groups = 60) {
  avail_geps <- intersect(gep_cols_new, colnames(df_x))
  if (length(avail_geps) == 0 || !target_var %in% colnames(df_x)) return(invisible(NULL))

  df_grp <- df_x[!is.na(df_x[[target_var]]), c(target_var, avail_geps), drop = FALSE]
  if (nrow(df_grp) == 0) return(invisible(NULL))

  # Group mean per level (for id_col each participant is already 1 row)
  mat_df <- tryCatch({
    df_grp %>%
      group_by(.data[[target_var]]) %>%
      summarise(across(all_of(avail_geps), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
      as.data.frame()
  }, error = function(e) {
    warning(sprintf("[row %d] sample_heatmap groupby %s/%s failed: %s",
                    row_idx, target_var, type_label, conditionMessage(e)))
    NULL
  })
  if (is.null(mat_df)) return(invisible(NULL))

  row_labels <- as.character(mat_df[[target_var]])
  mat_vals   <- as.matrix(mat_df[, avail_geps, drop = FALSE])
  rownames(mat_vals) <- row_labels

  n_rows_mat  <- nrow(mat_vals)
  n_geps_plot <- ncol(mat_vals)

  if (n_rows_mat < 2) {
    cat(sprintf("  Sample heatmap %s/%s: only %d group(s); skipping\n",
                target_var, type_label, n_rows_mat))
    return(invisible(NULL))
  }

  # Hierarchical clustering of rows
  hc_row <- tryCatch(
    hclust(dist(mat_vals, method = "euclidean"), method = "average"),
    error = function(e) NULL
  )
  if (!is.null(hc_row)) {
    row_order <- order.dendrogram(as.dendrogram(hc_row))
    mat_vals  <- mat_vals[row_order, , drop = FALSE]
  }

  # Column z-score
  z_mat <- mat_vals
  for (j in seq_len(ncol(z_mat))) {
    col_mean <- mean(z_mat[, j], na.rm = TRUE)
    col_sd   <- sd(z_mat[, j],   na.rm = TRUE)
    z_mat[, j] <- if (!is.na(col_sd) && col_sd > 0) {
      (z_mat[, j] - col_mean) / col_sd
    } else {
      0
    }
  }

  show_labels <- n_rows_mat <= max_label_groups
  w_px <- max(500L, as.integer(n_geps_plot * 45))
  h_px <- max(400L, as.integer(n_rows_mat * 20 + 100))

  varname_safe <- gsub("[^A-Za-z0-9_.-]", "_", target_var)

  make_mat_heatmap <- function(m, fill_type) {
    df_long <- tidyr::pivot_longer(
      rn_to_col(as.data.frame(m), "group"),
      cols = -"group", names_to = "GEP", values_to = "value"
    )
    df_long$group <- factor(df_long$group, levels = rev(rownames(m)))
    df_long$GEP   <- factor(df_long$GEP,   levels = colnames(m))

    p <- ggplot(df_long, aes(x = GEP, y = group, fill = value)) +
      geom_tile(color = "white", linewidth = 0.2) +
      labs(title = sprintf("%s GEP loading by %s \u2014 %s", fill_type, target_var, type_label),
           x = "GEP", y = target_var) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid  = element_blank())

    if (!show_labels) {
      p <- p + theme(axis.text.y = element_blank())
    }

    if (fill_type == "Absolute") {
      p <- p + scale_fill_viridis_c(name = "GEP\nscore", na.value = "grey90")
    } else {
      p <- p + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                                    midpoint = 0, limits = c(-2.5, 2.5),
                                    oob = scales::squish,
                                    name = "z-score", na.value = "grey90")
    }
    p
  }

  # Absolute heatmap
  tryCatch({
    p_abs    <- make_mat_heatmap(mat_vals, "Absolute")
    base_abs <- file.path(out_subdir,
                          sprintf("sample_heatmap_abs_%s_%s", varname_safe, type_label))
    save_plot(p_abs, base_abs, w_px, h_px)
    write.csv(rn_to_col(as.data.frame(mat_vals), "group"),
              paste0(base_abs, ".csv"), row.names = FALSE)
    cat(sprintf("  Saved: sample_heatmap_abs_%s_%s.svg+png+csv\n", varname_safe, type_label))
  }, error = function(e) {
    warning(sprintf("[row %d] abs heatmap %s/%s failed: %s",
                    row_idx, target_var, type_label, conditionMessage(e)))
  })

  # Z-scored heatmap
  tryCatch({
    p_z   <- make_mat_heatmap(z_mat, "Z-score")
    base_z <- file.path(out_subdir,
                        sprintf("sample_heatmap_zscore_%s_%s", varname_safe, type_label))
    save_plot(p_z, base_z, w_px, h_px)
    write.csv(rn_to_col(as.data.frame(z_mat), "group"),
              paste0(base_z, ".csv"), row.names = FALSE)
    cat(sprintf("  Saved: sample_heatmap_zscore_%s_%s.svg+png+csv\n", varname_safe, type_label))
  }, error = function(e) {
    warning(sprintf("[row %d] zscore heatmap %s/%s failed: %s",
                    row_idx, target_var, type_label, conditionMessage(e)))
  })

  invisible(NULL)
}

# ── Volcano plot (beta regression types only) ─────────────────────────────────
make_volcano_plot <- function(results, type_label, out_subdir, row_idx) {
  df_v <- results[results$gep_type == type_label &
                    !is.na(results$p_value) &
                    results$term != "(Intercept)", ]
  if (nrow(df_v) == 0) return(invisible(NULL))

  df_v$neg_log10_p <- -log10(pmax(df_v$p_value, .Machine$double.eps))

  n_terms <- length(unique(df_v$term))
  p <- ggplot(df_v, aes(x = estimate, y = neg_log10_p, color = gep_name)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    facet_wrap(~ term, scales = "free_x") +
    labs(title = sprintf("Volcano plot \u2014 %s", type_label),
         x = "Estimate", y = "-log10(p-value)", color = "GEP") +
    theme_minimal(base_size = 9)

  w_px <- max(600L, as.integer(n_terms * 250))
  save_plot(p, file.path(out_subdir, sprintf("volcano_%s", type_label)), w_px, 500L)
  cat(sprintf("  Saved: volcano_%s.svg + .png\n", type_label))
  invisible(NULL)
}

# ── GEP distribution boxplots (one per binary covariate per type) ─────────────
make_distribution_plots <- function(df_x, gep_cols_new, covs, type_label, out_subdir, row_idx) {
  avail_geps <- intersect(gep_cols_new, colnames(df_x))
  if (length(avail_geps) == 0) return(invisible(NULL))

  for (cv in covs) {
    if (!cv %in% colnames(df_x)) next
    vals_nna <- df_x[[cv]][!is.na(df_x[[cv]])]
    if (length(unique(vals_nna)) != 2) next

    df_long <- df_x[, c(cv, avail_geps), drop = FALSE] %>%
      tidyr::pivot_longer(all_of(avail_geps), names_to = "GEP", values_to = "score")
    df_long$GEP    <- factor(df_long$GEP, levels = avail_geps)
    df_long[[cv]]  <- as.character(df_long[[cv]])

    p <- ggplot(df_long, aes(x = GEP, y = score, fill = .data[[cv]])) +
      geom_boxplot(outlier.size = 0.8, linewidth = 0.4) +
      labs(title = sprintf("GEP scores by %s \u2014 %s", cv, type_label),
           x = "GEP", y = "Score", fill = cv) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    w_px     <- max(600L, as.integer(length(avail_geps) * 40 + 100))
    cov_safe <- gsub("[^A-Za-z0-9_.-]", "_", cv)
    save_plot(p, file.path(out_subdir, sprintf("gep_distribution_%s_%s", cov_safe, type_label)),
              w_px, 450L)
    cat(sprintf("  Saved: gep_distribution_%s_%s.svg + .png\n", cov_safe, type_label))
  }
  invisible(NULL)
}

# ── Coefficient heatmap (beta regression types only) ──────────────────────────
make_coef_heatmap <- function(results, type_label, out_subdir, row_idx) {
  df_c <- results[results$gep_type == type_label &
                    !is.na(results$estimate) &
                    results$term != "(Intercept)", ]
  if (nrow(df_c) == 0) return(invisible(NULL))

  terms_order <- unique(df_c$term)
  geps_order  <- unique(df_c$gep_name)
  df_c$term     <- factor(df_c$term,     levels = terms_order)
  df_c$gep_name <- factor(df_c$gep_name, levels = geps_order)

  max_abs <- max(abs(df_c$estimate), na.rm = TRUE)
  if (max_abs == 0) max_abs <- 1

  p <- ggplot(df_c, aes(x = term, y = gep_name, fill = estimate)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, limits = c(-max_abs, max_abs),
                         name = "Estimate") +
    labs(title = sprintf("Coefficient heatmap \u2014 %s", type_label),
         x = "Term", y = "GEP") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())

  n_terms <- length(terms_order)
  n_geps  <- length(geps_order)
  w_px <- max(400L, as.integer(n_terms * 60 + 150))
  h_px <- max(400L, as.integer(n_geps  * 25 + 100))

  base_path <- file.path(out_subdir, sprintf("coef_heatmap_%s", type_label))
  save_plot(p, base_path, w_px, h_px)

  out_wide <- tidyr::pivot_wider(df_c[, c("gep_name", "term", "estimate")],
                                  names_from = "term", values_from = "estimate")
  write.csv(out_wide, paste0(base_path, ".csv"), row.names = FALSE)
  cat(sprintf("  Saved: coef_heatmap_%s.svg+png+csv\n", type_label))
  invisible(NULL)
}

# ── Pseudo-R² / R² bar chart (all types combined) ────────────────────────────
make_r2_barplot <- function(results, out_subdir, row_idx) {
  # One row per GEP per type: use the (Intercept) row since pseudo_r2 is replicated
  df_r2 <- results[results$term == "(Intercept)", c("gep_name", "gep_type", "pseudo_r2")]
  if (nrow(df_r2) == 0) {
    df_r2 <- results %>%
      group_by(gep_name, gep_type) %>%
      summarise(pseudo_r2 = first(pseudo_r2), .groups = "drop") %>%
      as.data.frame()
  }
  df_r2 <- df_r2[!is.na(df_r2$pseudo_r2), ]
  if (nrow(df_r2) == 0) return(invisible(NULL))

  df_r2$gep_name <- factor(df_r2$gep_name)
  df_r2$gep_type <- factor(df_r2$gep_type,
                            levels = c("raw", "normalized", "normalized_clean"))

  p <- ggplot(df_r2, aes(x = gep_name, y = pseudo_r2, fill = gep_type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c("raw"              = "#4393c3",
                                 "normalized"       = "#d6604d",
                                 "normalized_clean" = "#74c476"),
                      name = "GEP type") +
    labs(title = "R\u00b2 / Pseudo-R\u00b2 by GEP and normalization",
         x = "GEP", y = "R\u00b2 / Pseudo-R\u00b2") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  n_geps <- length(unique(df_r2$gep_name))
  w_px   <- max(600L, as.integer(n_geps * 50 + 150))
  base_path <- file.path(out_subdir, "r2_barplot")
  save_plot(p, base_path, w_px, 450L)
  write.csv(df_r2, paste0(base_path, ".csv"), row.names = FALSE)
  cat(sprintf("  Saved: r2_barplot.svg+png+csv\n"))
  invisible(NULL)
}

# ── Participant PCA biplot ─────────────────────────────────────────────────────
make_pca_biplot <- function(df_pca, gep_cols_new, covs, type_label, out_subdir, row_idx) {
  avail_geps <- intersect(gep_cols_new, colnames(df_pca))
  if (length(avail_geps) < 2) return(invisible(NULL))

  gep_mat <- as.matrix(df_pca[, avail_geps, drop = FALSE])
  ok_rows <- complete.cases(gep_mat)
  if (sum(ok_rows) < 3) return(invisible(NULL))
  gep_mat <- gep_mat[ok_rows, , drop = FALSE]

  pca <- tryCatch(prcomp(gep_mat, center = TRUE, scale. = TRUE),
                  error = function(e) NULL)
  if (is.null(pca)) return(invisible(NULL))

  n_pcs  <- min(2L, ncol(pca$x))
  scores <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
  if (n_pcs < 2) scores$PC2 <- 0
  scores$participant <- rownames(gep_mat)

  # Identify a binary covariate for coloring (first one found)
  color_col <- NULL
  for (cv in covs) {
    if (!cv %in% colnames(df_pca)) next
    cv_sub   <- df_pca[[cv]][ok_rows]
    n_unique <- length(unique(cv_sub[!is.na(cv_sub)]))
    if (n_unique == 2) {
      color_col      <- cv
      scores[[cv]]   <- as.character(cv_sub)
      break
    }
  }

  # Variance explained
  pct_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pc1_lab <- sprintf("PC1 (%.1f%%)", pct_var[1])
  pc2_lab <- if (length(pct_var) >= 2) sprintf("PC2 (%.1f%%)", pct_var[2]) else "PC2"

  # Top 5 GEP loadings by combined PC1+PC2 magnitude
  load_mat <- pca$rotation[, seq_len(n_pcs), drop = FALSE]
  if (n_pcs < 2) load_mat <- cbind(load_mat, PC2 = 0)
  load_mag  <- rowSums(load_mat[, c("PC1", "PC2"), drop = FALSE]^2)
  top_geps  <- names(sort(load_mag, decreasing = TRUE))[seq_len(min(5L, length(load_mag)))]

  score_scale <- max(abs(scores[, c("PC1", "PC2")]), na.rm = TRUE)
  load_scale  <- max(sqrt(rowSums(load_mat[top_geps, , drop = FALSE]^2)), na.rm = TRUE)
  arrow_scale <- score_scale * 0.7 / max(load_scale, 1e-6)

  load_df <- data.frame(
    GEP = top_geps,
    PC1 = load_mat[top_geps, "PC1"] * arrow_scale,
    PC2 = load_mat[top_geps, "PC2"] * arrow_scale
  )

  # Base scatter
  if (!is.null(color_col)) {
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = .data[[color_col]])) +
      geom_point(size = 2, alpha = 0.8)
  } else {
    p <- ggplot(scores, aes(x = PC1, y = PC2)) +
      geom_point(size = 2, alpha = 0.8, color = "#4393c3")
  }

  p <- p +
    geom_segment(data = load_df, inherit.aes = FALSE,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "grey30", linewidth = 0.6) +
    geom_text(data = load_df, inherit.aes = FALSE,
              aes(x = PC1 * 1.12, y = PC2 * 1.12, label = GEP),
              size = 2.5, color = "grey20") +
    labs(title   = sprintf("PCA biplot \u2014 %s", type_label),
         x       = pc1_lab,
         y       = pc2_lab,
         color   = color_col) +
    theme_minimal(base_size = 9)

  base_path <- file.path(out_subdir, sprintf("pca_biplot_%s", type_label))
  save_plot(p, base_path, 700L, 600L)
  write.csv(scores, paste0(base_path, ".csv"), row.names = FALSE)
  cat(sprintf("  Saved: pca_biplot_%s.svg+png+csv\n", type_label))
  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# Section 4 — Per-row processing function
# ══════════════════════════════════════════════════════════════════════════════
process_row <- function(row_idx) {
  row         <- manifest[row_idx, ]
  qs_path     <- row$qs_path
  fl_path     <- row$feature_loadings_path
  qc_geps_raw <- row$qc_geps
  qc_geps <- if (is.na(qc_geps_raw) || !nzchar(trimws(qc_geps_raw))) {
    character(0)
  } else {
    trimws(strsplit(qc_geps_raw, ",")[[1]])
  }
  id_col <- trimws(row$id_col)
  covs   <- trimws(strsplit(row$covs, ",")[[1]])

  qs_basename <- tools::file_path_sans_ext(basename(qs_path))
  out_subdir  <- file.path(opt$output_dir, sprintf("%d_%s", row_idx, qs_basename))
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n[Row %d] %s\n", row_idx, qs_basename))

  # ── 3a — Load Seurat object ────────────────────────────────────────────────
  if (!file.exists(qs_path))
    stop(sprintf("ERROR [row %d]: qs_path not found: %s", row_idx, qs_path))
  qs_ext <- tolower(tools::file_ext(qs_path))
  if (!qs_ext %in% c("qs", "qs2"))
    stop(sprintf("ERROR [row %d]: qs_path must be .qs or .qs2, got: %s", row_idx, qs_path))

  seurat_obj <- tryCatch({
    if (qs_ext == "qs2") qs2::qs_read(qs_path) else qs::qread(qs_path)
  }, error = function(e)
    stop(sprintf("ERROR [row %d] reading Seurat: %s", row_idx, conditionMessage(e))))

  if (!inherits(seurat_obj, "Seurat"))
    stop(sprintf("ERROR [row %d]: loaded object is not Seurat (got %s)",
                 row_idx, class(seurat_obj)[1]))
  cat(sprintf("  Loaded Seurat: %d cells\n", ncol(seurat_obj)))

  # ── 3b — Load feature loadings ─────────────────────────────────────────────
  if (!file.exists(fl_path))
    stop(sprintf("ERROR [row %d]: feature_loadings_path not found: %s", row_idx, fl_path))

  fl <- tryCatch(
    read.table(fl_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
               quote = "", comment.char = "", check.names = FALSE),
    error = function(e)
      stop(sprintf("ERROR [row %d] reading feature loadings: %s", row_idx, conditionMessage(e)))
  )

  n_geps <- ncol(fl) - 1L
  if (n_geps < 1L) stop(sprintf("ERROR [row %d]: feature loadings has < 2 columns", row_idx))
  gep_cols_new <- sprintf("GEP_%02d", seq_len(n_geps))
  colnames(fl) <- c("library__barcode", gep_cols_new)

  # Resolve qc_geps: integer index or GEP_XX name
  qc_geps_resolved <- character(0)
  for (qc in qc_geps) {
    idx <- suppressWarnings(as.integer(qc))
    if (!is.na(idx)) {
      if (idx < 1L || idx > n_geps)
        stop(sprintf("ERROR [row %d]: qc_gep index %d out of range (1..%d)",
                     row_idx, idx, n_geps))
      qc_geps_resolved <- c(qc_geps_resolved, sprintf("GEP_%02d", idx))
    } else if (qc %in% gep_cols_new) {
      qc_geps_resolved <- c(qc_geps_resolved, qc)
    } else {
      stop(sprintf("ERROR [row %d]: qc_gep '%s' is not a valid integer index or GEP name (GEP_01..GEP_%02d)",
                   row_idx, qc, n_geps))
    }
  }
  qc_geps     <- unique(qc_geps_resolved)
  has_qc_geps <- length(qc_geps) > 0

  cat(sprintf("  Feature loadings: %d cells x %d GEPs  (QC GEPs: %s)\n",
              nrow(fl), n_geps,
              if (has_qc_geps) paste(qc_geps, collapse = ", ")
              else "none \u2014 normalized_clean will be skipped"))

  # ── 3c — Join feature loadings onto meta.data ──────────────────────────────
  meta_df <- seurat_obj@meta.data
  meta_df$library__barcode <- rownames(meta_df)

  joined <- inner_join(meta_df, fl, by = "library__barcode")
  if (nrow(joined) == 0)
    stop(sprintf("ERROR [row %d]: inner join produced 0 rows.", row_idx))
  cat(sprintf("  After inner join: %d cells retained (of %d in Seurat, %d in loadings)\n",
              nrow(joined), nrow(meta_df), nrow(fl)))

  # ── 4 — Validate columns ────────────────────────────────────────────────────
  joined_cols <- colnames(joined)
  if (!id_col %in% joined_cols)
    stop(sprintf("ERROR [row %d]: id_col '%s' not found. Available: %s",
                 row_idx, id_col, paste(joined_cols, collapse = ", ")))

  missing_covs <- setdiff(covs, joined_cols)
  if (length(missing_covs) > 0)
    stop(sprintf("ERROR [row %d]: covariate(s) not found: %s\nAvailable: %s",
                 row_idx, paste(missing_covs, collapse = ", "),
                 paste(joined_cols, collapse = ", ")))

  key_cols <- c(id_col, covs)
  na_mask  <- rowSums(is.na(joined[, key_cols, drop = FALSE])) > 0
  n_na     <- sum(na_mask)
  if (n_na > 0) {
    warning(sprintf("[row %d] Dropping %d / %d cells (%.1f%%) with NA in id_col or covariates",
                    row_idx, n_na, nrow(joined), 100 * n_na / nrow(joined)))
    joined <- joined[!na_mask, ]
  }
  if (nrow(joined) == 0)
    stop(sprintf("ERROR [row %d]: 0 cells remain after dropping NA rows", row_idx))

  n_participants <- length(unique(joined[[id_col]]))
  if (n_participants < 3)
    stop(sprintf("ERROR [row %d]: only %d participants remain; need >= 3",
                 row_idx, n_participants))

  # ── 5 — Cell-level GEP normalizations ──────────────────────────────────────
  gep_mat  <- as.matrix(joined[, gep_cols_new, drop = FALSE])
  row_sums <- rowSums(gep_mat)

  if (any(row_sums == 0)) {
    warning(sprintf("[row %d] %d cells have GEP row sum = 0; replacing with uniform 1/n_geps",
                    row_idx, sum(row_sums == 0)))
    gep_mat[row_sums == 0, ] <- 1 / n_geps
    row_sums[row_sums == 0] <- 1
  }
  norm_gep_mat <- sweep(gep_mat, 1, row_sums, "/")

  if (has_qc_geps) {
    clean_mat <- gep_mat
    clean_mat[, qc_geps] <- 0
    clean_row_sums <- rowSums(clean_mat)

    joined$clean_sum_zero <- clean_row_sums == 0
    zero_parts <- joined %>%
      group_by(.data[[id_col]]) %>%
      summarise(all_zero = all(clean_sum_zero), .groups = "drop") %>%
      filter(all_zero) %>%
      pull(.data[[id_col]])

    if (length(zero_parts) > 0)
      warning(sprintf("[row %d] Excluding %d participant(s) where all cells have clean GEP sum = 0: %s",
                      row_idx, length(zero_parts), paste(zero_parts, collapse = ", ")))

    n_clean_geps     <- max(1L, n_geps - length(qc_geps))
    zero_clean_idx   <- which(clean_row_sums == 0)
    if (length(zero_clean_idx) > 0) {
      non_qc <- setdiff(gep_cols_new, qc_geps)
      clean_mat[zero_clean_idx, non_qc] <- 1 / n_clean_geps
      clean_row_sums[zero_clean_idx]    <- 1
    }
    norm_clean_mat             <- sweep(clean_mat, 1, clean_row_sums, "/")
    colnames(norm_clean_mat)   <- paste0(gep_cols_new, "_norm_clean")
    joined$clean_sum_zero      <- NULL
    joined                     <- cbind(joined, norm_clean_mat)
  } else {
    zero_parts <- character(0)
  }

  colnames(norm_gep_mat) <- paste0(gep_cols_new, "_norm")
  joined <- cbind(joined, norm_gep_mat)

  # ── 6 — Aggregate to participant level ──────────────────────────────────────
  is_numeric_cov <- sapply(covs, function(cv) is.numeric(joined[[cv]]))
  num_covs <- covs[is_numeric_cov]
  cat_covs <- covs[!is_numeric_cov]

  all_gep_cols <- intersect(
    c(gep_cols_new, paste0(gep_cols_new, "_norm"), paste0(gep_cols_new, "_norm_clean")),
    colnames(joined)
  )

  df_agg_num <- joined %>%
    group_by(.data[[id_col]]) %>%
    summarise(across(all_of(c(num_covs, all_gep_cols)),
                     \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  if (length(cat_covs) > 0) {
    df_agg_cat <- joined %>%
      group_by(.data[[id_col]]) %>%
      summarise(across(all_of(cat_covs), majority_vote), .groups = "drop")
    df_participant <- left_join(df_agg_num, df_agg_cat, by = id_col)
  } else {
    df_participant <- df_agg_num
  }

  cat(sprintf("  Participant-level data: %d rows\n", nrow(df_participant)))

  # Build three participant-level data frames with clean GEP column names
  raw_gep_cols   <- gep_cols_new
  norm_gep_cols  <- paste0(gep_cols_new, "_norm")
  clean_gep_cols <- paste0(gep_cols_new, "_norm_clean")

  df_raw  <- df_participant[, c(id_col, covs, raw_gep_cols), drop = FALSE]
  df_norm <- df_participant[, c(id_col, covs, norm_gep_cols), drop = FALSE]
  colnames(df_norm)[colnames(df_norm) %in% norm_gep_cols] <- gep_cols_new

  if (has_qc_geps) {
    df_clean_rows <- if (length(zero_parts) > 0) {
      !df_participant[[id_col]] %in% zero_parts
    } else {
      rep(TRUE, nrow(df_participant))
    }
    df_clean <- df_participant[df_clean_rows, c(id_col, covs, clean_gep_cols), drop = FALSE]
    colnames(df_clean)[colnames(df_clean) %in% clean_gep_cols] <- gep_cols_new
  }

  # ── 7 — Collinearity detection ─────────────────────────────────────────────
  cat("  Checking collinearity ...\n")
  check_collinearity(df_raw,  covs, gep_cols_new, "raw")
  check_collinearity(df_norm, covs, gep_cols_new, "normalized")
  if (has_qc_geps) check_collinearity(df_clean, covs, gep_cols_new, "normalized_clean")

  # ── 7 — Regression ─────────────────────────────────────────────────────────
  # Raw: NNLS (non-negative constrained; appropriate for cNMF loadings)
  results_raw <- run_nnls_regressions(df_raw, covs, gep_cols_new, row_idx)

  # Normalized: beta regression
  results_norm  <- run_beta_regressions(df_norm,  "normalized",       covs, gep_cols_new, row_idx)
  results_clean <- if (has_qc_geps) {
    run_beta_regressions(df_clean, "normalized_clean", covs, gep_cols_new, row_idx)
  } else {
    cat("  Skipping normalized_clean (no QC GEPs specified)\n")
    NULL
  }

  all_results <- do.call(rbind, list(results_raw, results_norm, results_clean))

  if (is.null(all_results) || nrow(all_results) == 0)
    stop(sprintf("ERROR [row %d]: no regression results produced", row_idx))

  # ── Write gep_regression_results.csv ────────────────────────────────────────
  out_csv <- file.path(out_subdir, "gep_regression_results.csv")
  write.csv(all_results, out_csv, row.names = FALSE)
  cat(sprintf("  Saved: %s  (%d rows)\n", out_csv, nrow(all_results)))

  # ── Section 8 — Spearman correlations ──────────────────────────────────────
  cat("  Computing Spearman correlations ...\n")
  run_spearman_analysis(df_raw,  gep_cols_new, covs, "raw",              out_subdir, row_idx)
  run_spearman_analysis(df_norm, gep_cols_new, covs, "normalized",       out_subdir, row_idx)
  if (has_qc_geps)
    run_spearman_analysis(df_clean, gep_cols_new, covs, "normalized_clean", out_subdir, row_idx)

  # ── Section 9 — Sample heatmaps ────────────────────────────────────────────
  # Targets: id_col + categorical covariates with > 2 unique levels
  cat("  Generating sample heatmaps ...\n")
  sample_targets <- id_col
  for (cv in covs) {
    if (!cv %in% colnames(df_raw)) next
    vals_nna <- df_raw[[cv]][!is.na(df_raw[[cv]])]
    if (!is.numeric(vals_nna) && length(unique(vals_nna)) > 2)
      sample_targets <- c(sample_targets, cv)
  }

  for (target in sample_targets) {
    run_sample_heatmaps(df_raw,  target, gep_cols_new, "raw",              out_subdir, row_idx)
    run_sample_heatmaps(df_norm, target, gep_cols_new, "normalized",       out_subdir, row_idx)
    if (has_qc_geps)
      run_sample_heatmaps(df_clean, target, gep_cols_new, "normalized_clean", out_subdir, row_idx)
  }

  # ── Section 10 — Additional plots ──────────────────────────────────────────
  cat("  Generating additional plots ...\n")

  # 10a: Volcano plots (beta reg types only — NNLS has no p-values)
  for (type_label in c("normalized", "normalized_clean")) {
    tryCatch(make_volcano_plot(all_results, type_label, out_subdir, row_idx),
             error = function(e) warning(sprintf("[row %d] volcano %s: %s",
                                                 row_idx, type_label, conditionMessage(e))))
  }

  # 10b: GEP distribution boxplots (binary covariates, all types)
  for (type_label_df in list(list(df_raw, "raw"), list(df_norm, "normalized"))) {
    tryCatch(make_distribution_plots(type_label_df[[1]], gep_cols_new, covs,
                                     type_label_df[[2]], out_subdir, row_idx),
             error = function(e) warning(sprintf("[row %d] distributions %s: %s",
                                                 row_idx, type_label_df[[2]], conditionMessage(e))))
  }
  if (has_qc_geps) {
    tryCatch(make_distribution_plots(df_clean, gep_cols_new, covs,
                                     "normalized_clean", out_subdir, row_idx),
             error = function(e) warning(sprintf("[row %d] distributions normalized_clean: %s",
                                                 row_idx, conditionMessage(e))))
  }

  # 10c: Coefficient heatmaps (beta reg types only)
  for (type_label in c("normalized", "normalized_clean")) {
    tryCatch(make_coef_heatmap(all_results, type_label, out_subdir, row_idx),
             error = function(e) warning(sprintf("[row %d] coef_heatmap %s: %s",
                                                 row_idx, type_label, conditionMessage(e))))
  }

  # 10d: R² bar chart (all types combined)
  tryCatch(make_r2_barplot(all_results, out_subdir, row_idx),
           error = function(e) warning(sprintf("[row %d] r2_barplot: %s",
                                               row_idx, conditionMessage(e))))

  # 10e: PCA biplot — use normalized type; fall back to raw if normalized has < 3 rows
  pca_df   <- if (nrow(df_norm) >= 3) df_norm else df_raw
  pca_type <- if (nrow(df_norm) >= 3) "normalized" else "raw"
  tryCatch(make_pca_biplot(pca_df, gep_cols_new, covs, pca_type, out_subdir, row_idx),
           error = function(e) warning(sprintf("[row %d] pca_biplot: %s",
                                               row_idx, conditionMessage(e))))

  invisible(NULL)
}

# ══════════════════════════════════════════════════════════════════════════════
# Section 5 — Run all rows in parallel
# ══════════════════════════════════════════════════════════════════════════════
cat(sprintf("\n\u2500\u2500 Processing %d manifest row(s) with %d core(s) \u2500\u2500\n",
            n_rows, n_cores))

errors <- parallel::mclapply(
  seq_len(n_rows),
  function(i) tryCatch(
    { process_row(i); NULL },
    error = function(e) conditionMessage(e)
  ),
  mc.cores = n_cores
)

had_error <- FALSE
for (i in seq_along(errors)) {
  if (!is.null(errors[[i]])) {
    cat(sprintf("ERROR in row %d: %s\n", i, errors[[i]]))
    had_error <- TRUE
  }
}

if (had_error) quit(status = 1, save = "no")

cat(sprintf("\n\u2500\u2500 Done.  Outputs in: %s\n", normalizePath(opt$output_dir)))
