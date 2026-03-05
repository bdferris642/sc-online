# add_rank_purity.R
#
# Computes partial-residual rank-purity scores for each gene.
#
# Workflow
# --------
# 1. Pseudobulk logCPM per (donor × group), then pool groups per donor.
# 2. Build a donor × covariate design matrix from Seurat metadata, matching
#    base covariate names to CSV columns by left-prefix:
#      exact match  → numeric covariate  (e.g. "age"     → donor_meta$age)
#      suffix match → dummy variable     (e.g. "sexMale" → donor_meta$sex == "Male")
# 3. Partial residuals: subtract fitted covariate contributions (all except
#    contrast_col) from pseudobulked expression using CSV coefficients.
#    Residuals are median-centred per gene.
# 4. Rank-purity per gene:
#      categorical contrast → Mann-Whitney U, normalised to [-1, 1]
#      continuous contrast  → Spearman correlation
#    Sign check: if purity opposes sign of the CSV contrast coefficient, set to 0.
#    Report abs(purity).
#
# Returns a data.frame with columns: gene, rank_purity.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
})


# ── Main function ─────────────────────────────────────────────────────────────

compute_rank_purity = function(
    coef_df,       # data.frame: must have columns 'gene', <covariates>, <contrast_col>
                   #   covariate columns are log2FC per unit (numeric) or per dummy level
    covariates,    # character vector: base covariate names as they appear in seurat_obj metadata
                   #   CSV columns are matched by left-prefix (e.g. "sex" matches "sexMale")
    contrast_col,  # exact column name for the contrast of interest (in coef_df and metadata)
    seurat_obj,    # Seurat object: logCPM in layers["data"], covariates in meta.data
    groupby_col = NULL,  # metadata column defining cell groups for pseudobulking (e.g. cell_type)
                         #   NULL (default) pseudobulks by donor only
    donor_col,           # metadata column identifying donors/participants
    case_value = NULL  # for categorical contrast: which metadata level is "case" (mapped to 1)
) {
    # ── 1. Validate inputs ────────────────────────────────────────────────────
    if (!"gene" %in% colnames(coef_df)) {
        stop("coef_df must have a 'gene' column.")
    }

    # if contrast_col in covariates, remove it 
    covariates = setdiff(covariates, contrast_col)
    for (cov in covariates) {
        matched = colnames(coef_df)[startsWith(colnames(coef_df), cov)]
        matched = matched[matched != contrast_col]
        if (length(matched) == 0) {
            stop("No coef_df column found with prefix '", cov, "'.\n",
                 "  coef_df columns: ", paste(colnames(coef_df), collapse = ", "))
        }
        if (cov %in% matched && length(matched) > 1L) {
            stop(
                "Ambiguous: '", cov, "' matches both an exact and suffixed coef_df column: ",
                paste(matched, collapse = ", "), ".\n",
                "  Remove or rename the exact-match column."
            )
        }
    }
    

    meta = seurat_obj@meta.data
    if (is.null(groupby_col)) groupby_col = donor_col

    needed_meta = c(donor_col, groupby_col, contrast_col, covariates)
    missing_meta = setdiff(needed_meta, colnames(meta))
    if (length(missing_meta) > 0) {
        stop("Seurat metadata missing columns: ", paste(missing_meta, collapse = ", "))
    }

    # ── 2. Extract logCPM ─────────────────────────────────────────────────────
    logcpm = tryCatch(
        LayerData(seurat_obj, layer = "data"),
        error = function(e) GetAssayData(seurat_obj, slot = "data")
    )

    shared_cells = intersect(colnames(logcpm), rownames(meta))
    if (length(shared_cells) == 0) stop("No cells shared between logcpm and metadata.")
    logcpm = logcpm[, shared_cells, drop = FALSE]
    meta   = meta[shared_cells, , drop = FALSE]

    # ── 3. Pseudobulk ─────────────────────────────────────────────────────────
    message("Pseudobulking ...")

    meta$._group_key = paste(meta[[donor_col]], meta[[groupby_col]], sep = "|||")
    group_keys = unique(meta$._group_key)
    n_cells    = ncol(logcpm)

    grp_idx = match(meta$._group_key, group_keys)
    A = sparseMatrix(
        i = grp_idx, j = seq_len(n_cells), x = 1.0,
        dims = c(length(group_keys), n_cells),
        dimnames = list(group_keys, colnames(logcpm))
    )
    A = Diagonal(x = 1.0 / Matrix::rowSums(A)) %*% A

    pb_by_group    = logcpm %*% t(A)
    colnames(pb_by_group) = group_keys

    donor_of_group = sub("\\|\\|\\|.*$", "", group_keys)
    donors = unique(meta[[donor_col]])

    pb_by_donor = do.call(cbind, lapply(donors, function(d) {
        idx = which(donor_of_group == d)
        if (length(idx) == 1L) return(pb_by_group[, idx, drop = FALSE])
        Matrix::rowMeans(pb_by_group[, idx, drop = FALSE])
    }))
    colnames(pb_by_donor) = donors
    rm(pb_by_group); gc()

    # ── 4. Donor-level metadata ───────────────────────────────────────────────
    donor_meta = meta %>%
        group_by(across(all_of(donor_col))) %>%
        summarise(
            across(all_of(c(covariates, contrast_col)), function(x) {
                if (is.numeric(x)) mean(x, na.rm = TRUE) else x[[1L]]
            }),
            .groups = "drop"
        ) %>%
        as.data.frame()
    rownames(donor_meta) = donor_meta[[donor_col]]

    common_donors = intersect(colnames(pb_by_donor), donor_meta[[donor_col]])
    if (length(common_donors) == 0) stop("No donors in common between expression and metadata.")
    pb_by_donor = pb_by_donor[, common_donors, drop = FALSE]
    donor_meta  = donor_meta[common_donors, , drop = FALSE]

    # ── 5. Subset to genes present in both coef_df and expression data ────────
    common_genes = intersect(coef_df$gene, rownames(pb_by_donor))
    if (length(common_genes) == 0) {
        stop("No genes in common between coef_df and Seurat feature names.")
    }
    message(sprintf("  %d / %d genes found in Seurat data", length(common_genes), nrow(coef_df)))

    pb_sub   = as.matrix(pb_by_donor[common_genes, , drop = FALSE])
    coef_sub = coef_df[match(common_genes, coef_df$gene), , drop = FALSE]

    # ── 6. Partial residuals ──────────────────────────────────────────────────
    message("Computing partial residuals ...")

    design     = .expand_covariate_design(covariates, coef_sub, donor_meta, common_donors, contrast_col)
    fitted_cov = design$B %*% t(design$X)
    residuals  = pb_sub - fitted_cov

    gene_medians = apply(residuals, 1L, median, na.rm = TRUE)
    residuals    = sweep(residuals, 1L, gene_medians, FUN = "-")

    # ── 7. Rank-purity score ──────────────────────────────────────────────────
    message("Computing purity scores ...")

    contrast_vals  = donor_meta[[contrast_col]]
    is_categorical = !is.numeric(contrast_vals)

    if (is_categorical) {
        if (is.null(case_value)) {
            stop(
                "contrast_col '", contrast_col, "' is categorical but case_value is NULL.\n",
                "  Unique values: ", paste(sort(unique(contrast_vals)), collapse = ", ")
            )
        }
        if (!case_value %in% contrast_vals) {
            stop(
                "case_value '", case_value, "' not found in '", contrast_col, "'.\n",
                "  Unique values: ", paste(sort(unique(contrast_vals)), collapse = ", ")
            )
        }
        contrast_binary = as.integer(contrast_vals == case_value)
        rank_purity = apply(residuals, 1L, function(expr_vec) {
            .purity_score_categorical(expr_vec, contrast_binary)
        })
    } else {
        rank_purity = apply(residuals, 1L, function(expr_vec) {
            cor(expr_vec, contrast_vals, method = "spearman", use = "complete.obs")
        })
    }

    logfc_contrast = coef_sub[[contrast_col]]
    rank_purity[rank_purity * logfc_contrast < 0] = 0
    rank_purity = abs(rank_purity)

    # ── 8. Return ─────────────────────────────────────────────────────────────
    data.frame(gene = common_genes, rank_purity = rank_purity, row.names = NULL)
}


# ── Helpers ───────────────────────────────────────────────────────────────────

# Build expanded B (genes × k) and X (donors × k) matrices by matching base
# covariate names to coef_df columns via left-prefix:
#   suffix == ""  →  numeric: X column = donor_meta[[cov]]
#   suffix != ""  →  dummy:   X column = as.integer(donor_meta[[cov]] == suffix)
.expand_covariate_design = function(covariates, coef_sub, donor_meta, common_donors, contrast_col) {
    B_cols = list()
    X_cols = list()

    for (cov in covariates) {
        matched = colnames(coef_sub)[startsWith(colnames(coef_sub), cov)]
        matched = matched[matched != contrast_col]

        for (csv_col in matched) {
            suffix = substring(csv_col, nchar(cov) + 1L)

            if (suffix == "") {
                x_vec = donor_meta[[cov]]
                if (!is.numeric(x_vec)) {
                    stop(
                        "'", cov, "' has an exact-match coef_df column but is non-numeric ",
                        "in metadata (values: ", paste(head(unique(x_vec), 5), collapse = ", "), ").\n",
                        "  For categorical covariates use dummy columns in coef_df (e.g. '", cov, "<level>')."
                    )
                }
                message(sprintf("  [numeric]  %-20s  <- meta$%s", csv_col, cov))
                X_cols[[csv_col]] = as.numeric(x_vec)
            } else {
                meta_vals = donor_meta[[cov]]
                if (!suffix %in% meta_vals) {
                    warning(
                        "Level '", suffix, "' (from column '", csv_col, "') not found in ",
                        "metadata column '", cov, "' — indicator will be all-zero.\n",
                        "  Unique metadata values: ", paste(sort(unique(meta_vals)), collapse = ", ")
                    )
                }
                message(sprintf("  [dummy]    %-20s  <- (meta$%s == '%s')", csv_col, cov, suffix))
                X_cols[[csv_col]] = as.integer(meta_vals == suffix)
            }

            B_cols[[csv_col]] = coef_sub[[csv_col]]
        }
    }

    if (length(B_cols) == 0) stop("No covariate columns could be matched in coef_df.")

    B = do.call(cbind, B_cols)
    X = do.call(cbind, X_cols)
    rownames(X) = common_donors

    list(B = B, X = X)
}


# Mann-Whitney U purity score, normalised to [-1, 1].
.purity_score_categorical = function(x, group) {
    stopifnot(length(x) == length(group), all(group %in% c(0L, 1L)))

    rnks          = rank(x)
    case_ranks    = rnks[group == 1L]
    control_ranks = rnks[group == 0L]
    n_case        = length(case_ranks)
    n_control     = length(control_ranks)

    if (n_case == 0L || n_control == 0L) return(NA_real_)

    U = sum(outer(case_ranks, control_ranks, FUN = function(a, b) as.numeric(a > b)))
    2.0 * (U / (n_case * n_control)) - 1.0
}
