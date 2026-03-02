#!/usr/bin/env Rscript
# crumblr_analysis.R
#
# Runs crumblr proportion analysis on a (samples × cell_types) count matrix.
# Called by run_workflow.py; all arguments passed via command line.
#
# Args:
#   1  counts_csv      – samples × cell_types count matrix (row names = sample IDs)
#   2  metadata_csv    – samples × covariates (row names = sample IDs, same order)
#   3  formula_str     – R formula string, e.g. "~ case_control + age + (1|donor_id)"
#   4  contrast_col    – column name of the main contrast variable, e.g. "case_control"
#   5  output_dir      – directory for result CSVs and R plots

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage: crumblr_analysis.R counts_csv metadata_csv",
    "formula contrast_col output_dir"
  ))
}

counts_csv   <- args[1]
metadata_csv <- args[2]
formula_str  <- args[3]
contrast_col <- args[4]
output_dir   <- args[5]

# ── Library path set by caller via R_LIBS_CUSTOM env var ─────────────────────
custom_lib <- Sys.getenv("R_LIBS_CUSTOM", unset = NA_character_)
if (!is.na(custom_lib) && nchar(custom_lib) > 0) {
  .libPaths(c(custom_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(crumblr)
  library(variancePartition)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("── Loading data ────────────────────────────────────────────────────────\n")
counts   <- as.matrix(read.csv(counts_csv,   row.names = 1, check.names = FALSE))
metadata <- read.csv(metadata_csv, row.names = 1, check.names = FALSE)

# Guarantee matching row order
if (!identical(rownames(counts), rownames(metadata))) {
  stop("Row names of counts and metadata do not match after loading.")
}

cat(sprintf("  Samples: %d  |  Cell types: %d\n", nrow(counts), ncol(counts)))
cat(sprintf("  Formula: %s\n", formula_str))
cat(sprintf("  Contrast variable: %s\n", contrast_col))

# ── crumblr object ───────────────────────────────────────────────────────────
cat("── Running crumblr ─────────────────────────────────────────────────────\n")
cobj <- crumblr(counts)

# Export CLR-transformed values (features × samples in EList convention)
clr_df <- as.data.frame(t(cobj$E))    # → samples × cell_types
write.csv(clr_df, file.path(output_dir, "clr_values.csv"))

# Export precision weights
wt_df <- as.data.frame(t(cobj$weights))
write.csv(wt_df, file.path(output_dir, "precision_weights.csv"))

# ── dream model ──────────────────────────────────────────────────────────────
cat("── Fitting dream model ─────────────────────────────────────────────────\n")
form <- tryCatch(as.formula(formula_str), error = function(e) {
  stop("Invalid formula string '", formula_str, "': ", conditionMessage(e))
})

# dream needs the formula to contain a random effect for lme4; if none present,
# fall back to lmFit (fixed-effects only)
has_random <- grepl("\\|", formula_str)

# ── Collinearity pre-check (auto-drop redundant terms) ───────────────────────
# Strip random-effect terms to get the fixed-effects-only formula for checking.
fixed_form <- as.formula(sub("\\+\\s*\\([^)]+\\)", "", formula_str))
mm_check   <- tryCatch(model.matrix(fixed_form, data = metadata),
                       error = function(e) NULL)
if (!is.null(mm_check)) {
  qr_obj  <- qr(mm_check)
  rank_mm <- qr_obj$rank
  if (rank_mm < ncol(mm_check)) {
    term_labels     <- attr(terms(fixed_form), "term.labels")
    term_assignment <- attr(mm_check, "assign")

    # ── Identify structural cause before blindly trusting QR pivot ────────────
    # QR column pivoting finds *a* redundant column, but its choice is arbitrary
    # among equivalent solutions.  A common case it misidentifies: two binary
    # indicator columns that are perfect complements (col_i + col_j == 1 for all
    # samples, e.g. region_SN + region_VTA when every sample is in one region).
    # QR may pivot out a brain_bank dummy instead of the region column.
    # Check for this explicitly so we drop the right term.

    intercept_vec <- rep(1, nrow(mm_check))
    non_int_idx   <- which(colnames(mm_check) != "(Intercept)")

    # Binary columns: only 0s and 1s
    is_binary <- sapply(non_int_idx, function(j) all(mm_check[, j] %in% c(0, 1)))
    binary_idx <- non_int_idx[is_binary]

    drop_col_idx   <- NULL
    complement_msg <- NULL

    # Pairwise complement check: col_i + col_j == intercept
    if (length(binary_idx) >= 2) {
      for (ii in seq_len(length(binary_idx) - 1)) {
        for (jj in seq(ii + 1, length(binary_idx))) {
          ci <- binary_idx[ii]; cj <- binary_idx[jj]
          if (all(mm_check[, ci] + mm_check[, cj] == 1)) {
            # Drop the term associated with the LATER column (jj); the earlier
            # term is kept as it was listed first in --fixed-effects.
            drop_col_idx   <- cj
            complement_msg <- sprintf(
              "  Cause: '%s' + '%s' == 1 for all samples — perfect complement\n  indicators (every sample belongs to exactly one group).\n  Dropping '%s'; keep '%s' to represent the contrast.\n",
              colnames(mm_check)[ci], colnames(mm_check)[cj],
              colnames(mm_check)[cj], colnames(mm_check)[ci]
            )
            break
          }
        }
        if (!is.null(drop_col_idx)) break
      }
    }

    # Fall back to QR pivot if no complement pair found
    if (is.null(drop_col_idx)) {
      redundant_col_idx <- qr_obj$pivot[(rank_mm + 1L):ncol(mm_check)]
      drop_col_idx      <- redundant_col_idx[1]
      complement_msg    <- sprintf(
        "  QR-identified redundant column(s): %s\n",
        paste(colnames(mm_check)[redundant_col_idx], collapse = ", ")
      )
    }

    drop_term_idx <- term_assignment[drop_col_idx]
    drop_term     <- term_labels[drop_term_idx]

    # ── Print the exact linear dependency (null-space vector) ────────────────
    # For a rank-r matrix with n columns, the QR gives us R (upper triangular).
    # The last pivot column is the redundant one; solving R[1:r,1:r]*c = -R[1:r,r+1]
    # yields coefficients for: c1*col_1 + ... + cr*col_r + 1*col_redundant = 0
    tryCatch({
      R_mat <- qr.R(qr_obj)
      r     <- rank_mm
      dep_c <- backsolve(R_mat[seq_len(r), seq_len(r)],
                         -R_mat[seq_len(r), r + 1L])
      full_dep <- numeric(ncol(mm_check))
      full_dep[qr_obj$pivot[seq_len(r)]]  <- dep_c
      full_dep[qr_obj$pivot[r + 1L]]      <- 1.0
      nz      <- which(abs(full_dep) > 1e-8)
      dep_str <- paste(sprintf("%+.3f * %s", full_dep[nz], colnames(mm_check)[nz]),
                       collapse = "  ")
      cat("  Exact linear dependency:  ", dep_str, "= 0\n")
    }, error = function(e) NULL)

    cat("WARNING: Design matrix is rank-deficient — automatically dropping one term.\n")
    cat("  Rank:", rank_mm, "/", ncol(mm_check), "columns\n")
    cat(complement_msg)
    cat("  Dropping formula term:", drop_term, "\n")
    cat("  Tip: remove '", drop_term, "' from --fixed-effects to suppress this warning.\n\n",
        sep = "")

    # Rebuild formula: keep remaining fixed terms + all random terms.
    formula_rhs  <- as.character(form)[3]
    random_parts <- regmatches(formula_rhs,
                               gregexpr("\\([^)]+\\)", formula_rhs))[[1]]
    fixed_kept   <- setdiff(term_labels, drop_term)
    new_rhs      <- paste(c(fixed_kept, random_parts), collapse = " + ")
    form         <- as.formula(paste("~", new_rhs))
    cat("  Updated formula:", deparse(form), "\n\n")

    # Verify the cleaned matrix is now full rank.
    mm_verify <- tryCatch(
      model.matrix(as.formula(paste("~", paste(fixed_kept, collapse = " + "))),
                   data = metadata),
      error = function(e) NULL
    )
    if (!is.null(mm_verify) && qr(mm_verify)$rank < ncol(mm_verify)) {
      stop("Design matrix still rank-deficient after removing identified term. ",
           "There may be additional collinear covariates — inspect manually.")
    }
  }
}

if (has_random) {
  fit <- dream(cobj, form, metadata, BPPARAM = SerialParam())
} else {
  # Build model matrix manually for fixed-effects-only case
  mm  <- model.matrix(form, data = metadata)
  fit <- lmFit(cobj, mm)
}
fit <- eBayes(fit)

# ── Extract results for each level of the contrast variable ──────────────────
cat("── Extracting results ──────────────────────────────────────────────────\n")
all_coef   <- colnames(coef(fit))
test_coefs <- all_coef[grepl(contrast_col, all_coef, fixed = TRUE)]

if (length(test_coefs) == 0) {
  stop(sprintf(
    "No coefficients matching '%s' found. Available: %s",
    contrast_col, paste(all_coef, collapse = ", ")
  ))
}
cat("  Test coefficients:", paste(test_coefs, collapse = ", "), "\n")

results_list <- lapply(test_coefs, function(coeff) {
  res            <- topTable(fit, coef = coeff, number = Inf, sort.by = "none")
  res$cell_type  <- rownames(res)
  res$coefficient <- coeff
  res
})
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL
write.csv(results_df, file.path(output_dir, "crumblr_results.csv"), row.names = FALSE)

# Also write a tidy summary sorted by p-value
top_hits <- results_df[order(results_df$P.Value), ]
write.csv(top_hits, file.path(output_dir, "crumblr_results_sorted.csv"), row.names = FALSE)

# ── R-side forest / coefficient plot (per contrast coefficient) ───────────────
cat("── Generating R plots ──────────────────────────────────────────────────\n")

for (coeff in test_coefs) {
  sub <- results_df[results_df$coefficient == coeff, ]
  sub <- sub[order(sub$logFC), ]
  sub$cell_type <- factor(sub$cell_type, levels = sub$cell_type)
  sub$sig       <- ifelse(sub$adj.P.Val < 0.05, "FDR < 0.05", "n.s.")

  # SE is not a standard topTable column; derive it from logFC / t-statistic
  sub$SE <- with(sub, ifelse(abs(t) > 1e-10, abs(logFC / t), NA_real_))

  p <- ggplot(sub, aes(x = logFC, y = cell_type, colour = sig,
                       xmin = logFC - 1.96 * SE, xmax = logFC + 1.96 * SE)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbar(orientation = "y", width = 0.3, linewidth = 0.6) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("FDR < 0.05" = "#c0392b", "n.s." = "#7f8c8d")) +
    labs(
      title   = paste("crumblr effect sizes –", coeff),
      x       = "log2 fold-change (CLR scale)",
      y       = NULL,
      colour  = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  safe_coeff <- gsub("[^A-Za-z0-9_]", "_", coeff)
  ggsave(file.path(output_dir, paste0("forest_", safe_coeff, ".pdf")),
         p, width = 8, height = max(5, nrow(sub) * 0.30 + 2))
}

cat(sprintf(
  "── Done.  Results written to %s\n", normalizePath(output_dir)
))
