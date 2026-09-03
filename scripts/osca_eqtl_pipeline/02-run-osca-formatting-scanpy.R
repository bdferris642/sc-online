# 02-run-osca-formatting-scanpy.R — Format pseudobulk CSVs into OSCA input files.
#
# Arguments:
#   --expression-dir  / -e  (required) Directory with {cell_class}_expression_matrix_ds.csv files
#   --output-dir      / -o  (required) OSCA input/output directory (also contains bfiles)
#   --vcf-slogan      / -v  (required) Prefix of .fam and _pca.eigenvec files
#   --gene-anot       / -g  (optional) Gene annotation table; default /mnt/accessory/analysis/eqtl/gene_loc.txt
#   --metadata        / -m  (optional) CSV with participant_id and all covariate columns
#   --participants    / -p  (required) Path to participants text file
#   --sva-formula     / -s  (optional) Override SVA formula; auto-built from --h5ad-cat-covars + --h5ad-quant-covars
#   --h5ad-cat-covars / -c  (required unless --sva-formula given) Space-delimited categorical covariate names
#   --h5ad-quant-covars / -q (required unless --sva-formula given) Space-delimited quantitative covariate names
#
# Processing per cell class:
#   1. Loads expression CSV; subsets to participants in --participants
#   2. Guards: fails fast if no participants remain; fails fast if no WGS intersection
#   3. Removes low-variance genes (SD < 1e-6)
#   4. Runs SVA (n.sv = 5 surrogates) using the formula built from covariate lists
#   5. Runs PCA (30 PCs) on WGS-intersected expression data
#   6. Removes SVs with Spearman |r| > 0.9 to any other covariate; drops PCs 21-30
#   7. Annotates genes via gene_annotation table; drops unannotated genes

suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library("sva")))
suppressMessages(suppressWarnings(library(tidyverse)))


spec <- matrix(c(
    'expression-dir',   'e', 1, "character",
    'output-dir',       'o', 1, "character",
    'vcf-slogan',       'v', 1, "character",
    'gene-anot',        'g', 1, "character",
    'metadata',         'm', 1, "character",
    'participants',     'p', 1, "character",
    'sva-formula',      's', 1, "character",
    'h5ad-cat-covars',  'c', 1, "character",
    'h5ad-quant-covars','q', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
EXPRESSION_DIR = opt[['expression-dir']]
OUTPUT_DIR = opt[['output-dir']]
VCF_SLOGAN = opt[['vcf-slogan']]
PARTICIPANTS_PATH = opt[['participants']]

# Parse covariate lists (space-delimited strings → character vectors)
CAT_COVARS   = if (is.null(opt[['h5ad-cat-covars']]))   character(0) else strsplit(trimws(opt[['h5ad-cat-covars']]),   "\\s+")[[1]]
QUANT_COVARS = if (is.null(opt[['h5ad-quant-covars']])) character(0) else strsplit(trimws(opt[['h5ad-quant-covars']]), "\\s+")[[1]]
cat("Categorical covariates:", if (length(CAT_COVARS))   paste(CAT_COVARS,   collapse=", ") else "(none)", "\n")
cat("Quantitative covariates:", if (length(QUANT_COVARS)) paste(QUANT_COVARS, collapse=", ") else "(none)", "\n")

# FIX: SVA formula built from covariate lists; --sva-formula overrides if provided
if (!is.null(opt[['sva-formula']])) {
    SVA_FORMULA_STR = opt[['sva-formula']]
} else if (length(c(CAT_COVARS, QUANT_COVARS)) == 0) {
    stop("Must provide --h5ad-cat-covars, --h5ad-quant-covars, or --sva-formula.")
} else {
    SVA_FORMULA_STR = paste("~", paste(c(CAT_COVARS, QUANT_COVARS), collapse = " + "))
}
cat("SVA formula:", SVA_FORMULA_STR, "\n")

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR)
}

if (is.null(opt[['gene-anot']])) {
    GENE_ANNOT_PATH = "/mnt/accessory/analysis/eqtl/gene_loc.txt"
} else {
    GENE_ANNOT_PATH = opt[['gene-anot']]
}
# If --metadata is given, load once; otherwise auto-discover per-cell-class CSV from EXPRESSION_DIR.
METADATA_PATH = if (!is.null(opt[['metadata']])) opt[['metadata']] else NULL

gene_annotation = read.table(GENE_ANNOT_PATH, header = TRUE) %>% distinct(NAME, .keep_all = TRUE)
participants = readLines(PARTICIPANTS_PATH)
print(participants)

wgs_subset = read.table(file.path(OUTPUT_DIR, paste0(VCF_SLOGAN, ".fam")), header = FALSE)
genotype_pcs = read.table(file.path(OUTPUT_DIR, paste0(VCF_SLOGAN, "_pca.eigenvec"))) %>% rename(
        participant_id = V2, G_PC1 = V3, G_PC2 = V4, G_PC3 = V5, G_PC4 = V6, G_PC5 = V7)
genotype_pcs$V1 = NULL
head(genotype_pcs)
fam_file = read.table(file.path(OUTPUT_DIR, paste0(VCF_SLOGAN, ".fam")), header = FALSE)


all_files = list.files(EXPRESSION_DIR)
expression_files = all_files[grep("expression_matrix", all_files)]
composition_files = all_files[grep("composition_matrix", all_files)]
common_prefixes = gsub("_expression_matrix_ds.csv", "",basename(expression_files))
names(expression_files) = common_prefixes
names(composition_files) = common_prefixes
cat("\nCommon prefixes:", common_prefixes, "\n")

# Process each file
for (cc in common_prefixes) {
    print(cc)

    # Load metadata: per-cell-class CSV written by step 01, or global --metadata if provided.
    if (is.null(METADATA_PATH)) {
        meta_file = file.path(EXPRESSION_DIR, paste0(cc, "_obs_metadata.csv"))
        if (!file.exists(meta_file)) stop(paste0(
            "No metadata found for ", cc, ": expected ", meta_file, ". ",
            "Either provide --metadata or run step 01 to generate per-cell-class metadata."))
        metadata = read.csv(meta_file)
    } else {
        metadata = read.csv(METADATA_PATH)
    }

    efile = file.path(EXPRESSION_DIR, expression_files[[cc]])
    cfile = file.path(EXPRESSION_DIR, composition_files[[cc]])
      
    cat(paste0("\nLoading expression data from ", efile))
    edata = t(as.matrix(read.csv(efile)) %>% as.data.frame() %>% column_to_rownames("X"))

    # subset to participants of interest, enumerated in file in PARTICIPANTS_PATH
    participant_cols = colnames(edata)[colnames(edata) %in% participants]
    edata = edata[, participant_cols]
    if (ncol(edata) == 0) stop(paste0(
        "No participants from ", PARTICIPANTS_PATH, " match expression data columns for ", cc, "."))

    rn = rownames(edata)
    edata = apply(edata, 2, as.numeric)
    rownames(edata) = rn

    # remove constant rows, which should have low sd
    low_sd_rows = apply(edata, 1, sd) < 1e-6
    cat("\nRemoving ", sum(low_sd_rows), " low sd rows\n")
    edata = edata[!low_sd_rows, ] 
    cat("\n", dim(edata), "\n")
    print(head(edata))


    cat("\nFiltering Phenotypic Data:\n")
    # Filter metadata to match participant_ids in expression data
    pheno = metadata %>%
        filter(participant_id %in% colnames(edata)) %>%
        arrange(match(participant_id, colnames(edata)))
    rownames(pheno) = pheno$participant_id

    cat("\nRunning SVA:\n")
    # FIX: validate that all variables in the formula are present in pheno before running SVA
    sva_formula = as.formula(SVA_FORMULA_STR)
    sva_vars = all.vars(sva_formula)
    missing_vars = setdiff(sva_vars, colnames(pheno))
    if (length(missing_vars) > 0) {
        stop(paste0(
            "SVA formula '", SVA_FORMULA_STR, "' references variables not present in metadata: ",
            paste(missing_vars, collapse = ", "),
            ". Available columns: ", paste(colnames(pheno), collapse = ", ")
        ))
    }
    mod = model.matrix(sva_formula, data = pheno)
    mod0 = model.matrix(~ 1, data = pheno)
    svobj = sva(edata, mod, mod0, n.sv = 5)  # n.sv = 5: estimate 5 surrogate variables

    sv_factors = as.data.frame(svobj$sv)
    colnames(sv_factors) = paste0("SV", seq_len(ncol(sv_factors)))
    sv_col_names = colnames(sv_factors)
    pheno_with_svs = cbind(pheno, sv_factors)

    # Match participant_ids with genotype data
    matching_columns = wgs_subset$V2[wgs_subset$V2 %in% colnames(edata)]
    if (length(matching_columns) == 0) stop(paste0(
        "No samples in expression data for ", cc, " have matching genotype data in ", VCF_SLOGAN, ".fam"))
    edata = edata[, matching_columns, drop = FALSE]

    cat("\nRunning PCA:\n")
    pca_result = prcomp(t(edata), scale. = TRUE)
    top_30_pcs = pca_result$x[, 1:30]
    rownames(top_30_pcs) = colnames(edata)

    # Create phenotype and covariate data
    phenotype = as.data.frame(t(edata)) %>%
        rownames_to_column("participant_id")

    cat("\nMerging metadata with SVs:\n")
    covs1 = metadata %>%
        filter(participant_id %in% row.names(top_30_pcs)) %>%
        select(participant_id, all_of(c(CAT_COVARS, QUANT_COVARS))) %>%
        left_join(pheno_with_svs %>% select(participant_id, all_of(sv_col_names)),
                  by = "participant_id")

    pcs = as.data.frame(top_30_pcs) %>% rownames_to_column("participant_id")
    clusters =  (
        read.csv(cfile) %>%
        filter(X %in% colnames(edata))%>%
        rename(participant_id=X))

    print("Merging metadata and svs with pcs:")
    scaled_quant_covars = paste0(QUANT_COVARS, "_scaled")
    masterdf = merge(
        merge(
            merge(
                merge(covs1, pcs, by = "participant_id"),
                genotype_pcs, by = "participant_id"),
            phenotype, by = "participant_id") %>%
            mutate(across(all_of(QUANT_COVARS),
                          ~ (. - min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE)),
                          .names = "{.col}_scaled")),
        clusters, by = "participant_id")

    # Drop covariates with a single unique value across participants in this cell class.
    # Constant covariates cause rank deficiency in model.matrix() and must be excluded.
    const_cat = CAT_COVARS[sapply(CAT_COVARS, function(v) length(unique(na.omit(masterdf[[v]]))) < 2)]
    if (length(const_cat) > 0) {
        warning(paste0("[", cc, "] Dropping constant categorical covariate(s) (one unique value): ",
                       paste(const_cat, collapse=", ")))
        CAT_COVARS = setdiff(CAT_COVARS, const_cat)
    }
    const_quant = QUANT_COVARS[sapply(QUANT_COVARS, function(v) length(unique(na.omit(masterdf[[v]]))) < 2)]
    if (length(const_quant) > 0) {
        warning(paste0("[", cc, "] Dropping constant quantitative covariate(s) (one unique value): ",
                       paste(const_quant, collapse=", ")))
        QUANT_COVARS = setdiff(QUANT_COVARS, const_quant)
        scaled_quant_covars = paste0(QUANT_COVARS, "_scaled")
    }

    valid_columns = intersect(names(phenotype)[-1], names(masterdf))

    phenotype = masterdf %>%
        mutate(FID = 0, IID = participant_id) %>%
        select(FID, IID, all_of(valid_columns)) %>%
        select(-FID) %>%
        column_to_rownames("IID") %>%
        t() %>%
        as.data.frame()

    merged_data = merge(
            phenotype %>% rownames_to_column("NAME"),
            gene_annotation,
            by = "NAME",
            all.x = TRUE) %>%
        select(probe, chr, TSS, NAME, strand) %>%
        filter(!is.na(probe) & !is.na(chr) & !is.na(TSS))

    # Update phenotype and match order
    phenotype = phenotype[rownames(phenotype) %in% merged_data$NAME, ]
    phenotype = phenotype[match(merged_data$NAME, rownames(phenotype)), ]

    merged_data = merged_data %>%
    filter(NAME %in% rownames(phenotype)) %>%
    arrange(match(NAME, rownames(phenotype)))

    # Write outputs
    write.table(as.data.frame(t(phenotype)) %>% rownames_to_column("participant_id") %>%
                    mutate(FID = 0, IID = participant_id) %>% select(FID, IID, rownames(phenotype)),
                paste0(OUTPUT_DIR, "/Phenotype_", cc, "_osca.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)

    write.table(merged_data %>% select(chr, NAME, TSS, probe, strand),
                    paste0(OUTPUT_DIR, "/Upprobe_", cc, ".opi"),
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    write.table(masterdf %>%
                    mutate(FID = 0, IID = participant_id) %>%
                    mutate(across(all_of(CAT_COVARS), as.factor)) %>%
                    select(FID, IID, all_of(CAT_COVARS)),
                paste0(OUTPUT_DIR, "/cov1_", cc, ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)

    write.table(masterdf %>%
                    mutate(FID = 0, IID = participant_id) %>%
                    select(FID, IID, names(pcs)[-1],
                        all_of(c("G_PC1", "G_PC2", "G_PC3", "G_PC4", "G_PC5")),
                        all_of(sv_col_names),
                        all_of(scaled_quant_covars),
                        names(clusters)[-1]),
                paste0(OUTPUT_DIR, "/cov2_", cc, ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
}


for (cc in common_prefixes) {
    cat("Processing", cc, "\n")
    cov2 = read.table(paste0(OUTPUT_DIR, "/cov2_", cc, ".txt"), header = T, sep = "\t") 
    
    cor_matrix = cor(cov2[,! colnames(cov2) %in% c("FID", "IID")], method = "spearman") %>% as.data.frame()

    # remove SVs that are highly correlated with other covariates
    # and PCs 21-30
    sv_cols_to_remove = c()
    for (i in 1:nrow(cor_matrix)) {
        for (j in 1:ncol(cor_matrix)) {
            if (i == j) {
                next
            }
            if (cor_matrix[i, j] > 0.9) {  # Spearman r > 0.9: remove SVs collinear with other covariates
                rowname = rownames(cor_matrix)[[i]]
                colname = colnames(cor_matrix)[[j]]
                two_cols = c(rowname, colname)
                print(paste(rowname, colname, cor_matrix[i, j]))
                sv_col = two_cols[grepl("SV", two_cols)]
                sv_cols_to_remove = c(sv_cols_to_remove, sv_col)
            }
        }
    }
    sv_cols_to_remove = unique(sv_cols_to_remove)

    cols_to_remove = c(
        "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30",
        sv_cols_to_remove
    )
    print(cols_to_remove)

    cov2_reduced = cov2[, !colnames(cov2) %in% cols_to_remove]
    write.table(cov2_reduced, paste0(OUTPUT_DIR, "/cov2_", cc, "_reduced.txt"), sep = "\t", quote = F, row.names = F)
}