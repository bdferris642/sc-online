suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library("sva")))
suppressMessages(suppressWarnings(library(tidyverse)))


spec <- matrix(c(
    'expression-dir', 'e', 1, "character",
    'output-dir', 'o', 1, "character",
    'vcf-slogan', 'v', 1, "character",
    'gene-anot', 'g', 1, "character",
    'metadata', 'm', 1, "character",
    'participants', 'p', 1, "character"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
EXPRESSION_DIR = opt[['expression-dir']]
OUTPUT_DIR = opt[['output-dir']]
VCF_SLOGAN = opt[['vcf-slogan']]
PARTICIPANTS_PATH = opt[['participants']]

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR)
}

if (is.null(opt[['gene-anot']])) {
    GENE_ANNOT_PATH = "/mnt/accessory/analysis/eqtl/gene_loc.txt"
} else {
    GENE_ANNOT_PATH = opt[['gene-anot']]
}
if (is.null(opt[['metadata']])) {
    METADATA_PATH = "/mnt/accessory/seq_data/pd_all/240514/dapi_nurr_merged_seurat_clean_subsets/dapi_nurr_metadata_new.csv"
} else {
    METADATA_PATH = opt[['metadata']]
}

metadata = read.csv(METADATA_PATH)
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
common_prefixes = gsub("dapi_nurr_merged_seurat_clean__" , "", common_prefixes)
names(expression_files) = common_prefixes
names(composition_files) = common_prefixes
cat("\nCommon prefixes:", common_prefixes, "\n")

# Process each file
for (cc in common_prefixes) {
    print(cc)
    # Subset only the NPH expression and composition files
    
    efile = file.path(EXPRESSION_DIR, expression_files[[cc]])
    cfile = file.path(EXPRESSION_DIR, composition_files[[cc]])
      
    cat(paste0("\nLoading expression data from ", efile))
    edata = t(as.matrix(read.csv(efile)) %>% as.data.frame() %>% column_to_rownames("X"))

    # subset to participants of interest, enumerated in file in PARTICIPANTS_PATH
    participant_cols = colnames(edata)[colnames(edata) %in% participants]
    edata = edata[, participant_cols]
    
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
    mod = model.matrix(~ age + sex + pmi, data = pheno)
    mod0 = model.matrix(~ 1, data = pheno)
    svobj = sva(edata, mod, mod0, n.sv = 5)

    sv_factors = as.data.frame(svobj$sv)
    colnames(sv_factors) = paste0("SV", seq_len(ncol(sv_factors)))
    pheno_with_svs = cbind(pheno, sv_factors)

    # Match participant_ids with genotype data
    matching_columns = wgs_subset$V2[wgs_subset$V2 %in% colnames(edata)]
    edata = edata[, matching_columns, drop = FALSE]

    cat("\nRunning PCA:\n")
    pca_result = prcomp(t(edata), scale. = TRUE)
    top_30_pcs = pca_result$x[, 1:30]
    rownames(top_30_pcs) = colnames(edata)

    # Create phenotype and covariate data
    phenotype = as.data.frame(t(edata)) %>%
        rownames_to_column("participant_id")

    cat("\nMerging metadata with SVs:\n")
    covs1 = merge(
                metadata %>% filter(participant_id %in% row.names(top_30_pcs)), 
                pheno_with_svs, 
                by = "participant_id") %>% 
        select(participant_id, sex.x, age.x, pmi.x, SV1, SV2, SV3, SV4, SV5) %>%
        rename(sex = sex.x, age = age.x, pmi=pmi.x)

    pcs = as.data.frame(top_30_pcs) %>% rownames_to_column("participant_id")
    clusters =  (
        read.csv(cfile) %>%
        filter(X %in% colnames(edata))%>%
        rename(participant_id=X))

    print("Merging metadata and svs with pcs:")
    masterdf = merge(
        merge(
            merge(
                merge(covs1, pcs, by = "participant_id"), 
                genotype_pcs, by = "participant_id"), 
            phenotype, by = "participant_id") %>% 
            mutate(
                age_scaled = (age - min(age)) / (max(age) - min(age)),
                pmi_scaled = (pmi - min(pmi)) / (max(pmi) - min(pmi))),
        clusters, by="participant_id")

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
                    select(FID, IID, sex) %>%
                    mutate(sex = as.factor(sex)),
                paste0(OUTPUT_DIR, "/cov1_", cc, ".txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)

    write.table(masterdf %>%
                    mutate(FID = 0, IID = participant_id) %>%
                    select(FID, IID, names(pcs)[-1],
                        "G_PC1", "G_PC2", "G_PC3", "G_PC4", "G_PC5",
                        "SV1", "SV2", "SV3", "SV4", "SV5", 
                        age_scaled, pmi_scaled, names(clusters)[-1]),
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
            if (cor_matrix[i, j] > 0.9) {
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
    write.table(cov2_reduced, paste0(OUTPUT_DIR, "cov2_", cc, "_reduced.txt"), sep = "\t", quote = F, row.names = F)
}

