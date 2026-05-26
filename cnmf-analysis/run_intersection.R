print("**************** LOADING LIBRARIES ****************")
# Detect script path when running via Rscript
args = commandArgs(trailingOnly = FALSE)
script_path = sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 1) {
  script_dir = dirname(normalizePath(script_path))
  message("Script located in directory: ", script_dir)
} else {
  stop("Cannot determine script path. Are you running via Rscript?")
}

# run-intersection.R
# Given an intersection input TSV (gene_name, gene_id, GEP_01, ...) produced by
# get-cnmf-spectra-intersection, runs Fisher's exact test + BH correction for each
# gene set collection, per GEP specified by --rank-col.
# Outputs per-collection CSVs and SVG dot plots.

# Paths to shared resources
GSEA_R       = file.path(script_dir, "gsea.R")
GENESETS_DIR = "/home/ferris/sc-online/genesets/latest"

suppressWarnings(suppressMessages({
    library(dplyr)
    library(ggplot2)
    library(glue)
    library(getopt)
    library(stringr)
    library(scales)
    source(GSEA_R)
}))

# Gene set collections — same as run-gsea.R
gene_sets = list(
    kegg_2026_human     = file.path(GENESETS_DIR, "KEGG_2026.gmt"),
    reactome_2024       = file.path(GENESETS_DIR, "Reactome_Pathways_2024.gmt"),
    gwas_2025           = file.path(GENESETS_DIR, "GWAS_Catalog_2025.gmt"),
    jensen_disease_2025 = file.path(GENESETS_DIR, "Jensen_DISEASES_Curated_2025.gmt"),
    trrust_2019         = file.path(GENESETS_DIR, "TRRUST_Transcription_Factors_2019.gmt"),
    jaspar_tf_2025      = file.path(GENESETS_DIR, "JASPAR_PWM_Human_2025.gmt"),
    go_bp_2025          = file.path(GENESETS_DIR, "GO_Biological_Process_2025.gmt"),
    go_mf_2025          = file.path(GENESETS_DIR, "GO_Molecular_Function_2025.gmt")
)

# ── Argument parsing ──────────────────────────────────────────────────────────
spec <- matrix(c(
    'path',     'p', 1, "character",
    'rank-col', 'r', 1, "character",
    'top-n',    'n', 1, "integer",
    'zscore',   'z', 1, "double",
    'min-size', 'i', 1, "integer",
    'max-size', 'x', 1, "integer"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)

if (is.null(opt[['path']]))     stop("--path is required")
if (is.null(opt[['rank-col']])) stop("--rank-col is required")

PATH      = opt[['path']]
RANK_COL  = opt[['rank-col']]
TOP_N     = if (is.null(opt[['top-n']]))    150L  else opt[['top-n']]
ZSCORE    = if (is.null(opt[['zscore']]))   3.0   else opt[['zscore']]
MIN_SIZE  = if (is.null(opt[['min-size']])) 15L   else opt[['min-size']]
MAX_SIZE  = if (is.null(opt[['max-size']])) 250L  else opt[['max-size']]

message("Path:     ", PATH)
message("Rank col: ", RANK_COL)
message("Top-N:    ", TOP_N)
message("Z-score:  ", ZSCORE)
message("Min size: ", MIN_SIZE)
message("Max size: ", MAX_SIZE)

# ── Load input ────────────────────────────────────────────────────────────────
input_df = read.delim(PATH, stringsAsFactors = FALSE, check.names = FALSE)

if (!(RANK_COL %in% colnames(input_df))) {
    stop(glue("Column '{RANK_COL}' not found in {PATH}"))
}

# Background: all rows with non-NA gene_name
background_df = input_df[!is.na(input_df$gene_name) & input_df$gene_name != "", ]
n_background  = nrow(background_df)
message("Background genes: ", n_background)

# ── Gene selection: z-score + top-n cap ───────────────────────────────────────
values  = background_df[[RANK_COL]]
z_vals  = (values - mean(values, na.rm = TRUE)) / sd(values, na.rm = TRUE)

passing_idx = which(z_vals > ZSCORE)
message("Genes passing z-score threshold (>", ZSCORE, "): ", length(passing_idx))

if (length(passing_idx) > TOP_N) {
    # Keep only the top-N highest-loading genes among z-score passers
    ranked_passing = passing_idx[order(values[passing_idx], decreasing = TRUE)]
    passing_idx    = ranked_passing[seq_len(TOP_N)]
    message("Capped to top-N: ", TOP_N)
}

top_genes = background_df$gene_name[passing_idx]
top_genes = unique(top_genes[!is.na(top_genes) & top_genes != ""])
message("Final top gene count: ", length(top_genes))

if (length(top_genes) == 0) {
    message("No genes selected for ", RANK_COL, ". Skipping.")
    quit(status = 0)
}

# ── Derive output directories from input path ─────────────────────────────────
# Input path:  <cache_dir>/<basename_noext>-intersection-input.tsv
# GEP outputs: <cache_dir>/intersection/<RANK_COL>/
cache_dir      = dirname(PATH)
input_basename = tools::file_path_sans_ext(basename(PATH))
# strip trailing "-intersection-input"
base_stem      = sub("-intersection-input$", "", input_basename)

gep_outdir       = file.path(cache_dir, "intersection", RANK_COL)
fig_svg_dir      = file.path(gep_outdir, "figures", "svg")
fig_svg_filt_dir = file.path(fig_svg_dir, "redundancy_filtered")

dir.create(gep_outdir,       recursive = TRUE, showWarnings = FALSE)
dir.create(fig_svg_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(fig_svg_filt_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: wrap long pathway names ──────────────────────────────────────────
wrap_pathway = function(s, width = 45) {
    stringr::str_wrap(s, width = width)
}

# ── Helper: dot plot ──────────────────────────────────────────────────────────
plot_intersection_dotplot = function(res_df, title, svg_path, top_n = 15) {
    sig_df = res_df[!is.na(res_df$padj) & res_df$padj < 0.05, , drop = FALSE]
    if (nrow(sig_df) == 0) {
        message("  No significant pathways (padj<0.05) for plot: ", title)
        return(invisible(NULL))
    }

    # Sort by odds_ratio descending, keep top_n
    sig_df = sig_df[order(sig_df$odds_ratio, decreasing = TRUE), ]
    if (nrow(sig_df) > top_n) sig_df = sig_df[seq_len(top_n), ]

    sig_df$log2_OR        = log2(pmax(sig_df$odds_ratio, 1e-6))
    sig_df$neg_log10_padj = pmin(-log10(pmax(sig_df$padj, 1e-10)), 10)
    sig_df$pathway_wrap   = sapply(sig_df$pathway, wrap_pathway)
    sig_df$pathway_wrap   = factor(
        sig_df$pathway_wrap,
        levels = rev(unique(sig_df$pathway_wrap))
    )

    p = ggplot(sig_df, aes(x = log2_OR, y = pathway_wrap,
                           size = n_overlap, color = neg_log10_padj)) +
        geom_point(alpha = 0.85) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        scale_color_viridis_c(
            name = "-log10(padj)", limits = c(0, 10),
            oob = scales::squish
        ) +
        scale_size_continuous(name = "n overlap", range = c(2, 8)) +
        labs(title = title, x = "log2(Odds Ratio)", y = NULL) +
        theme_bw(base_size = 10) +
        theme(
            plot.title      = element_text(size = 9, face = "bold"),
            axis.text.y     = element_text(size = 7),
            axis.text.x     = element_text(size = 8),
            legend.position = "right"
        )

    n_paths = nrow(sig_df)
    fig_h   = max(3, min(n_paths * 0.45 + 2, 20))
    fig_w   = 8

    ggsave(svg_path, plot = p, width = fig_w, height = fig_h,
           device = "svg", limitsize = FALSE)
    message("  Saved: ", svg_path)
}

# ── Run Fisher's exact test per collection ────────────────────────────────────
all_bg_genes = background_df$gene_name

for (collection_name in names(gene_sets)) {
    geneset_path = gene_sets[[collection_name]]

    if (!file.exists(geneset_path)) {
        message("Gene set file not found, skipping: ", geneset_path)
        next
    }

    message("\n── Collection: ", collection_name, " ──")

    gs = .sconline.GSEA.readGMT(
        file        = geneset_path,
        bkg_genes   = all_bg_genes,
        min.gs.size = MIN_SIZE,
        max.gs.size = MAX_SIZE
    )
    message("  Pathways after size filter: ", length(gs))

    if (length(gs) == 0) {
        message("  No pathways remaining after size filter, skipping.")
        next
    }

    # Fisher's exact test for each pathway
    results = lapply(names(gs), function(pw_name) {
        pathway_genes = gs[[pw_name]]
        pathway_genes = pathway_genes[!is.na(pathway_genes) & pathway_genes != ""]

        a = length(intersect(top_genes,   pathway_genes))  # top ∩ pathway
        b = length(top_genes)    - a                        # top \ pathway
        c = length(pathway_genes) - a                       # pathway \ top
        d = n_background - a - b - c                        # rest

        if (a == 0) {
            return(data.frame(
                pathway       = pw_name,
                n_overlap     = 0L,
                n_top_genes   = length(top_genes),
                n_pathway     = length(pathway_genes),
                n_background  = n_background,
                odds_ratio    = 0,
                ci_lower      = 0,
                ci_upper      = Inf,
                pvalue        = 1,
                padj          = NA_real_,
                overlap_genes = "",
                stringsAsFactors = FALSE
            ))
        }

        ft = tryCatch(
            fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater"),
            error = function(e) NULL
        )
        if (is.null(ft)) {
            return(data.frame(
                pathway       = pw_name,
                n_overlap     = a,
                n_top_genes   = length(top_genes),
                n_pathway     = length(pathway_genes),
                n_background  = n_background,
                odds_ratio    = NA_real_,
                ci_lower      = NA_real_,
                ci_upper      = NA_real_,
                pvalue        = NA_real_,
                padj          = NA_real_,
                overlap_genes = paste(intersect(top_genes, pathway_genes), collapse = ","),
                stringsAsFactors = FALSE
            ))
        }

        data.frame(
            pathway       = pw_name,
            n_overlap     = a,
            n_top_genes   = length(top_genes),
            n_pathway     = length(pathway_genes),
            n_background  = n_background,
            odds_ratio    = ft$estimate,
            ci_lower      = ft$conf.int[1],
            ci_upper      = ft$conf.int[2],
            pvalue        = ft$p.value,
            padj          = NA_real_,
            overlap_genes = paste(intersect(top_genes, pathway_genes), collapse = ","),
            stringsAsFactors = FALSE
        )
    })

    res_df = do.call(rbind, results)

    # BH correction across all pathways in this collection
    res_df$padj = p.adjust(res_df$pvalue, method = "BH")
    res_df      = res_df[order(res_df$padj, na.last = TRUE), ]

    n_sig = sum(!is.na(res_df$padj) & res_df$padj < 0.05)
    message("  Significant pathways (padj<0.05): ", n_sig)

    # Save CSV
    out_csv = file.path(
        gep_outdir,
        glue("{base_stem}_{RANK_COL}_{collection_name}_intersection.csv")
    )
    write.csv(res_df, out_csv, row.names = FALSE)
    message("  Saved CSV: ", out_csv)

    # Per-GEP dot plot (all significant)
    out_svg = file.path(
        fig_svg_dir,
        glue("{base_stem}_{RANK_COL}_{collection_name}_intersection.svg")
    )
    plot_intersection_dotplot(
        res_df,
        title    = glue("{RANK_COL} | {collection_name} | Fisher intersection"),
        svg_path = out_svg
    )

    # ── Redundancy filtering ──────────────────────────────────────────────────
    # Adapt res_df to match filter_redundant_gene_sets() expectations:
    #   leadingEdge <- overlap_genes (comma-string)
    #   NES         <- log2(odds_ratio)  (all positive; only enrichment tested)
    #   padj        <- padj
    #   pathway     <- pathway
    if (n_sig > 0) {
        filt_input             = res_df
        filt_input$leadingEdge = filt_input$overlap_genes
        filt_input$NES         = log2(pmax(filt_input$odds_ratio, 1e-6))
        filt_input$NES[!is.finite(filt_input$NES)] = 0

        filt_df = tryCatch(
            filter_redundant_gene_sets(filt_input),
            error = function(e) {
                message("  Redundancy filter error: ", e$message, ". Using unfiltered.")
                filt_input
            }
        )

        # Drop helper columns added by filter_redundant_gene_sets or adapter
        drop_cols = c("leadingEdge", "NES", "._comp_degree")
        filt_df   = filt_df[, setdiff(colnames(filt_df), drop_cols), drop = FALSE]

        # Save filtered CSV
        out_csv_filt = file.path(
            gep_outdir,
            glue("{base_stem}_{RANK_COL}_{collection_name}_intersection_redundancy_filtered.csv")
        )
        write.csv(filt_df, out_csv_filt, row.names = FALSE)
        message("  Saved filtered CSV: ", out_csv_filt)

        # Plot filtered
        out_svg_filt = file.path(
            fig_svg_filt_dir,
            glue("{base_stem}_{RANK_COL}_{collection_name}_intersection_redundancy_filtered.svg")
        )
        plot_intersection_dotplot(
            filt_df,
            title    = glue("{RANK_COL} | {collection_name} | Fisher intersection (redundancy filtered)"),
            svg_path = out_svg_filt
        )
    }
}

message("\nDone: ", RANK_COL)
