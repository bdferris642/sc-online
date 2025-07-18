#!/bin/bash

# FILES=(
#   "/mnt/accessory/seq_data/gtex/score-eqtl/sox6_agtr1_eya4.mean_upr_sn.smr.txt"
#   "/mnt/accessory/seq_data/gtex/score-eqtl/upr_astro.mean_upr_sn.smr.txt"
#   "/mnt/accessory/seq_data/gtex/score-eqtl/upr_mg.mean_upr_sn.smr.txt"
#   "/mnt/accessory/seq_data/gtex/score-eqtl/upr_oligo.mean_upr_sn.smr.txt"
#   "/mnt/accessory/seq_data/gtex/score-eqtl/upr_opc.mean_upr_sn.smr.txt"
# )
# FILES=(
#     "/mnt/accessory/seq_data/gtex/score-eqtl/sox6_agtr1_eya4.mean_upr_sn.eqtl_for_besd.txt"
# )

# "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_astro.filtered.tsv"
#   "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_da.filtered.tsv"
#   "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_mg.filtered.tsv"
#   "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_nonda.filtered.tsv"

FILES=(
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_astro.filtered.tsv"
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_da.filtered.tsv"
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_da__SOX6_AGTR1.filtered.tsv"
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_mg.filtered.tsv"
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_nonda.filtered.tsv"
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_oligo.filtered.tsv"
    "/mnt/accessory/seq_data/gtex/eqtl_final_outs/gtex_osca_0.01/eqtl_opc.filtered.tsv"
)

BFILE="/mnt/accessory/seq_data/gtex/gtex_merged_v8_v9_normed_dbsnp_annot_filtered_dedupe_rs_only"
GWAS="/home/ferris/filtered_PD_GWAS.ma"

# Loop through each eQTL file and run SMR pipeline
for FILE in "${FILES[@]}"; do
    echo "Processing $(basename "$FILE")"

    # Strip .eqtl_for_besd.txt extension to get the prefix
    prefix="${FILE%.tsv}"

    # Step 1: Create BESD file if it doesn't exist
    if [[ ! -f "${prefix}.besd" ]]; then
        echo "Creating BESD file for $(basename "$FILE")"
        smr --qfile "$FILE" --make-besd --out "$prefix"
    else
        echo "BESD file already exists for $(basename "$FILE"), skipping creation."
    fi

    # Step 2: Run SMR with GWAS summary
    smr --bfile "$BFILE" \
        --gwas-summary "$GWAS" \
        --beqtl-summary "$prefix" \
        --out "${prefix}_PDgwaspublic_smr"

    echo "Finished: $(basename "$FILE")"
done

echo "All files processed successfully!"



# EQTL_FILE=$1
# BED_FILE_PREFIX=$2
# GWAS_FILE=$3

# # EQTL_FILE: Path to the eQTL summary statistics file. Best guess at format is a text file with headers
# # SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b se p


# # format GWAS into SMR compatible format
# # SNP A1 A2 Probe B SE P

# # pd_gwas = read.csv("~/gene_sets_ucell/kim-pd-gwas-2024-s3.csv")

# # pd_gwas_re = pd_gwas[pd_gwas$Method == "Random-effect",c("rsID", "CHR", "BP", "A1", "A2", "Nearest.Gene.Feature", "BETA.RE.", "P.RE.")]
# # pd_gwas_mega = pd_gwas[pd_gwas$Method == "MR-MEGA", c("rsID", "CHR", "BP", "A1", "A2", "Nearest.Gene.Feature", "BETA.MR.MEGA.", "P.MR-MEGA.")]

# # TODO: determine the name of the BFILE

# # Step 1: Create BESD file
# ./smr --qfile "$EQTL_FILE" --make-besd --out "$EQTL_FILE"

# # Step 2: Run SMR analysis
# ./smr --bfile "$BED_FILE_PREFIX" --gwas-summary $GWAS_FILE --beqtl-summary "$EQTL_FILE" --out "${EQTL_FILE}_smr"

# echo "SMR processed successfully!"

