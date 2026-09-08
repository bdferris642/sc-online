#!/usr/bin/env Rscript
# regenerate_gene_loc.R — Rebuild gene_loc.txt from current Ensembl annotation.
#
# Queries biomaRt for all Ensembl IDs present in the expression data and
# writes a gene_loc.txt compatible with OSCA's --gene-anot format:
#   probe  chr  TSS  NAME  strand
# where probe = Entrez gene ID, TSS = transcription start site, chr = 1-22 only.
#
# Usage:
#   Rscript regenerate_gene_loc.R \
#     --expr-csv  /path/to/any_expression_matrix_ds.csv \
#     --ensg-map  /mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv \
#     --out       /mnt/accessory/analysis/eqtl/gene_loc.txt   # or a new writable path
#
# Requirements: biomaRt (install via bioconductor-biocmart in the env if absent)

suppressMessages(library(getopt))
suppressMessages(library(dplyr))

spec <- matrix(c(
    'expr-csv',  'e', 1, "character",
    'ensg-map',  'm', 1, "character",
    'out',       'o', 1, "character",
    'ensembl-version', 'v', 2, "integer"   # optional; omit for current release
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

EXPR_CSV        <- opt[["expr-csv"]]
ENSG_MAP        <- opt[["ensg-map"]]
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

# ── 2. Query biomaRt ──────────────────────────────────────────────────────────
suppressMessages(library(biomaRt))

cat("Connecting to Ensembl biomaRt ...\n")
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

cat("Querying gene annotations ...\n")
res <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id",
                   "chromosome_name", "transcription_start_site", "strand"),
    filters    = "ensembl_gene_id",
    values     = ensg_ids,
    mart       = mart
)
cat(sprintf("  Retrieved %d rows for %d unique Ensembl IDs\n",
            nrow(res), length(unique(res$ensembl_gene_id))))

# ── 3. Clean and format ───────────────────────────────────────────────────────
gene_loc <- res %>%
    # keep only autosomes
    filter(chromosome_name %in% as.character(1:22)) %>%
    # require Entrez ID and symbol
    filter(!is.na(entrezgene_id) & hgnc_symbol != "") %>%
    mutate(
        probe  = entrezgene_id,
        chr    = as.integer(chromosome_name),
        TSS    = transcription_start_site,
        NAME   = hgnc_symbol,
        strand = ifelse(strand == 1, "+", "-")
    ) %>%
    select(probe, chr, TSS, NAME, strand) %>%
    # one row per symbol: if multiple TSS entries, take the most 5' TSS per strand
    group_by(NAME) %>%
    arrange(ifelse(strand == "+", TSS, -TSS)) %>%
    slice(1) %>%
    ungroup() %>%
    distinct(probe, .keep_all = TRUE) %>%   # one row per Entrez ID
    arrange(chr, TSS)

cat(sprintf("Final gene_loc: %d genes across chromosomes 1-22\n", nrow(gene_loc)))

# ── 4. Write output ───────────────────────────────────────────────────────────
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
write.table(gene_loc, OUT_PATH, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Written to: %s\n", OUT_PATH))

# ── 5. Coverage report ────────────────────────────────────────────────────────
if (!is.null(ENSG_MAP)) {
    ensg_sym <- read.csv(ENSG_MAP) %>%
        filter(hgnc_symbol != "", !is.na(hgnc_symbol)) %>%
        distinct(ensembl_gene_id, .keep_all = TRUE)
    mapped_syms  <- ensg_sym %>% filter(ensembl_gene_id %in% ensg_ids) %>% pull(hgnc_symbol)
    in_loc       <- sum(mapped_syms %in% gene_loc$NAME)
    cat(sprintf("\nCoverage: %d / %d expressed+mapped genes now in gene_loc (%.1f%%)\n",
                in_loc, length(mapped_syms), 100 * in_loc / length(mapped_syms)))
}
