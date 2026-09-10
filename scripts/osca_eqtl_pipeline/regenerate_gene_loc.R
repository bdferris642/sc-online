#!/usr/bin/env Rscript
# regenerate_gene_loc.R — Rebuild gene_loc_v2.txt from current Ensembl annotation.
#
# v2: Keyed by Ensembl gene ID (ensg_id). No Entrez ID required — zero gene loss
# from the Entrez join that affected v1. Gene symbol is best-effort (falls back
# to ensg_id if hgnc_symbol is absent).
#
# Output columns (tab-delimited, with header):
#   ensg_id     — Ensembl gene ID (primary key)
#   chr         — chromosome (autosomes 1-22 only)
#   TSS         — transcription start site (bp)
#   gene_symbol — HGNC symbol, or ensg_id as fallback if unmapped
#   strand      — "+" or "-"
#
# Usage:
#   Rscript regenerate_gene_loc.R \
#     --expr-csv  /path/to/any_expression_matrix_ds.csv \
#     --out       /path/to/gene_loc_v2.txt \
#     [--ensembl-version 109]
#
# Requirements: biomaRt (bioconductor-biocmart in the env)

suppressMessages(library(getopt))
suppressMessages(library(dplyr))

spec <- matrix(c(
    'expr-csv',        'e', 1, "character",
    'out',             'o', 1, "character",
    'ensembl-version', 'v', 2, "integer"   # optional; omit for current release
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

EXPR_CSV        <- opt[["expr-csv"]]
OUT_PATH        <- opt[["out"]]
ENSEMBL_VERSION <- opt[["ensembl-version"]]  # NULL → current

if (is.null(EXPR_CSV) || is.null(OUT_PATH)) {
    cat(getopt(spec, usage = TRUE))
    quit(status = 1)
}

# ── 1. Collect all Ensembl IDs from the expression CSV header ─────────────────
cat("Reading Ensembl IDs from expression CSV ...\n")
header <- read.csv(EXPR_CSV, nrows = 0)
ensg_ids <- colnames(header)
ensg_ids <- ensg_ids[grepl("^ENSG", ensg_ids)]
cat(sprintf("  %d Ensembl IDs found\n", length(ensg_ids)))

if (length(ensg_ids) == 0) {
    stop("No ENSG* column names found in the expression CSV. ",
         "Ensure the CSV has ENSG IDs as column names (participant_id in first col).")
}

# ── 2. Query annotation source ────────────────────────────────────────────────
# Prefer EnsDb.Hsapiens.v86 (local, no network). Fall back to biomaRt if not installed.

use_ensdb <- requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) &&
             requireNamespace("ensembldb",           quietly = TRUE)

if (use_ensdb) {
    suppressMessages(library(EnsDb.Hsapiens.v86))
    suppressMessages(library(ensembldb))
    cat("Querying EnsDb.Hsapiens.v86 (local — no network required) ...\n")
    txs <- transcripts(EnsDb.Hsapiens.v86,
        filter      = GeneIdFilter(ensg_ids),
        return.type = "data.frame",
        columns     = c("gene_id", "gene_name", "tx_id",
                        "tx_seq_start", "tx_seq_end", "seq_name", "seq_strand"))
    cat(sprintf("  Retrieved %d transcripts for %d unique genes\n",
                nrow(txs), length(unique(txs$gene_id))))

    txs <- txs %>%
        mutate(TSS = ifelse(seq_strand == 1, tx_seq_start, tx_seq_end))

    gene_loc <- txs %>%
        filter(seq_name %in% as.character(1:22)) %>%
        group_by(gene_id) %>%
        summarise(
            chr            = first(seq_name),
            strand_num     = first(seq_strand),
            gene_name_raw  = first(gene_name),
            TSS            = ifelse(first(seq_strand) == 1, min(TSS), max(TSS)),
            .groups = "drop"
        ) %>%
        mutate(
            ensg_id     = gene_id,
            strand      = ifelse(strand_num == 1, "+", "-"),
            gene_symbol = ifelse(!is.na(gene_name_raw) & gene_name_raw != "",
                                 gene_name_raw, gene_id)
        ) %>%
        select(ensg_id, chr, TSS, gene_symbol, strand) %>%
        distinct(ensg_id, .keep_all = TRUE) %>%
        arrange(as.integer(chr), TSS)

} else {
    suppressMessages(library(biomaRt))
    cat("EnsDb.Hsapiens.v86 not found — falling back to biomaRt (requires network) ...\n")
    if (!is.null(ENSEMBL_VERSION)) {
        host <- sprintf("https://e%d.ensembl.org", ENSEMBL_VERSION)
        mart <- tryCatch(
            useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host),
            error = function(e) useEnsembl("genes", dataset = "hsapiens_gene_ensembl",
                                           version = ENSEMBL_VERSION)
        )
        cat(sprintf("  Using Ensembl version %d\n", ENSEMBL_VERSION))
    } else {
        mart <- tryCatch(
            useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                    host = "https://www.ensembl.org"),
            error = function(e) useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
        )
        cat("  Using current Ensembl release\n")
    }

    res <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol",
                       "chromosome_name", "transcription_start_site", "strand"),
        filters    = "ensembl_gene_id",
        values     = ensg_ids,
        mart       = mart
    )
    cat(sprintf("  Retrieved %d rows for %d unique Ensembl IDs\n",
                nrow(res), length(unique(res$ensembl_gene_id))))

    gene_loc <- res %>%
        filter(chromosome_name %in% as.character(1:22)) %>%
        mutate(
            ensg_id     = ensembl_gene_id,
            chr         = chromosome_name,
            TSS         = transcription_start_site,
            gene_symbol = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "",
                                 hgnc_symbol, ensembl_gene_id),
            strand      = ifelse(strand == 1, "+", "-")
        ) %>%
        select(ensg_id, chr, TSS, gene_symbol, strand) %>%
        group_by(ensg_id) %>%
        arrange(ifelse(strand == "+", TSS, -TSS)) %>%
        slice(1) %>%
        ungroup() %>%
        arrange(as.integer(chr), TSS)
}

n_retained  <- nrow(gene_loc)
n_with_sym  <- sum(gene_loc$gene_symbol != gene_loc$ensg_id)
n_fallback  <- n_retained - n_with_sym
cat(sprintf("Final gene_loc_v2: %d genes across chromosomes 1-22\n", n_retained))
cat(sprintf("  %d with HGNC symbol, %d using ENSG fallback for gene_symbol\n",
            n_with_sym, n_fallback))

# Coverage vs input gene list
n_input   <- length(ensg_ids)
n_covered <- sum(ensg_ids %in% gene_loc$ensg_id)
cat(sprintf("Coverage: %d / %d input genes retained (%.1f%%) — no Entrez filter applied\n",
            n_covered, n_input, 100 * n_covered / n_input))

# ── 4. Write output ───────────────────────────────────────────────────────────
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
write.table(gene_loc, OUT_PATH, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Written to: %s\n", OUT_PATH))
