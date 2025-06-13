#!/bin/bash

# start with h5ad files, eg in /mnt/accessory/seq_data/pd_all/240514/dapi_nurr_merged_seurat_clean_subsets
# make PB csvs, eg in /mnt/accessory/seq_data/pd_all/240514/dapi_nurr_merged_seurat_clean_subsets/{pseudobulk}
# put scanpy outputs into a different directory 
# containing OSCA inputs, eg in /mnt/accessory/seq_data/gtex/eqtl_inputs


# STEP 1: make pseudobulk csvs of gene expression and cell composition
# STEP 2: format these into OSCA inputs, along with associated metadata
# STEP 3: run OSCA on all cell classes simultaneously
# STEP 4: copy OSCA tsv outputs to Google Cloud Storage
# STEP 5: process and plot OSCA outputs
# STEP 6: subset to common SNP-probes
# STEP 7: run mashr
# STEP 8: copy all OSCA inputs and outputs to Google Cloud Storage

# Default values
CC_FILE="/mnt/analysis/eqtl/gtex/ccs"
CT_ID="cell_class"
GENE_LOG_EXPR_THRESHOLD=1.5
MASH_EPS=1e-6
MASH_NUM_RANDOM=1000000
MASH_PADJ_THRESH=0.01
MIN_NUM_CELLS=10
PB_OUTPUT_SUBDIR="pseudobulk"
SAMPLE_ID="participant_id"
STRS_TO_SKIP="endo,opc"
START_AT_STEP=1
STOP_AFTER_STEP=1000

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cc-file ) CC_FILE="$2"; shift ;;
        --ct-id) CT_ID="$2"; shift ;;
        --expr-input-dir) GENE_EXPR_INPUT_DIR="$2"; shift ;;
        --expr-output-subdir) PB_OUTPUT_SUBDIR="$2"; shift ;;
        --gene-log-expr-threshold) GENE_LOG_EXPR_THRESHOLD="$2"; shift ;;
        --mash-eps) MASH_EPS="$2"; shift ;;
        --mash-num-random) MASH_NUM_RANDOM="$2"; shift ;;
        --mash-padj-thresh) MASH_PADJ_THRESH="$2"; shift ;;
        --min-num-cells) MIN_NUM_CELLS="$2"; shift ;;
        --osca-input-dir) OSCA_INPUT_DIR="$2"; shift ;;
        --participants) PARTICIPANT_FNAME="$2"; shift ;;  
        --pipeline-slogan) PIPELINE_SLOGAN="$2"; shift ;;
        --sample-id) SAMPLE_ID="$2"; shift ;;
        --start-at-step) START_AT_STEP="$2"; shift ;;
        --stop-after-step) STOP_AFTER_STEP="$2"; shift ;;
        --strs-to-skip) STRS_TO_SKIP="$2"; shift ;;
        --vcf-slogan) VCF_SLOGAN="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Required argument checks
: "${GENE_EXPR_INPUT_DIR:?--expr-input-dir is required}"
: "${OSCA_INPUT_DIR:?--osca-input-dir is required}"
: "${PIPELINE_SLOGAN:?--pipeline-slogan is required}"
: "${VCF_SLOGAN:?--vcf-slogan is required}"
: "${PARTICIPANT_FNAME:?--participants is required}"


# include date in log filenames
DATE_TIME=$(date +'%Y-%m-%d_%H-%M-%S')
GOOGLE_BUCKET=gs://macosko_data/ferris/eqtl_outs/$PIPELINE_SLOGAN
LOG_FNAME="/mnt/accessory/analysis/eqtl/logs/osca_eqtl_pipeline_$PIPELINE_SLOGAN_${DATE_TIME}.log"
SCRIPT_DIR=/home/ferris/sc-online/scripts/osca_eqtl_pipeline

# redirect stdout and stderr to log file and terminal
touch "$LOG_FNAME"
exec > >(tee -i "$LOG_FNAME") 2>&1

conda activate mashr

# store git branch and commit hash to log file so pipeline can be replicated later
branch=$(git rev-parse --abbrev-ref HEAD)
commit_hash=$(git rev-parse HEAD)

echo "Branch: $branch"
echo "Commit: $commit_hash"

OSCA_OUTPUT_DIR=$OSCA_INPUT_DIR/eqtl_final_outs/$PIPELINE_SLOGAN
PB_OUTPUT_DIR="$GENE_EXPR_INPUT_DIR/$PB_OUTPUT_SUBDIR"

if [ ! -d "$PB_OUTPUT_DIR" ]; then
    mkdir -p "$PB_OUTPUT_DIR"
fi

if [ ! -d "$OSCA_INPUT_DIR" ] || \
    [ ! -f "$OSCA_INPUT_DIR/$VCF_SLOGAN.bed" ] || \
    [ ! -f "$OSCA_INPUT_DIR/$VCF_SLOGAN.bim" ] || \
    [ ! -f "$OSCA_INPUT_DIR/$VCF_SLOGAN.fam" ]; then
    mkdir -p $OSCA_INPUT_DIR
    echo "Must add .bed/.bim./.fam of VCF files to $OSCA_INPUT_DIR before running Pipeline."
    exit 1
fi

if [ ! -d "$OSCA_OUTPUT_DIR" ]; then
    mkdir -p "$OSCA_OUTPUT_DIR"
fi

if [ ! -f "$PARTICIPANT_FNAME" ]; then
    echo "Must add line-delimited participant file to $PARTICIPANT_FNAME before running Pipeline."
    exit 1
fi

if [ ! -f "$CC_FILE" ]; then
    echo "CC file not found: $CC_FILE"
    exit 1
fi

# write all function arguments to log file
echo "Running OSCA eQTL pipeline with the following arguments:"
echo "CC_FILE: $CC_FILE"
echo "CT_ID: $CT_ID"
echo "GENE_EXPR_INPUT_DIR: $GENE_EXPR_INPUT_DIR"
echo "PB_OUTPUT_SUBDIR: $PB_OUTPUT_SUBDIR"
echo "GENE_LOG_EXPR_THRESHOLD: $GENE_LOG_EXPR_THRESHOLD"
echo "MASH_EPS: $MASH_EPS"
echo "MASH_NUM_RANDOM: $MASH_NUM_RANDOM"
echo "MASH_PADJ_THRESH: $MASH_PADJ_THRESH"
echo "MIN_NUM_CELLS: $MIN_NUM_CELLS"
echo "OSCA_INPUT_DIR: $OSCA_INPUT_DIR"
echo "PARTICIPANT_FNAME: $PARTICIPANT_FNAME"
echo "PIPELINE_SLOGAN: $PIPELINE_SLOGAN"
echo "SAMPLE_ID: $SAMPLE_ID"
echo "START_AT_STEP: $START_AT_STEP"
echo "STOP_AFTER_STEP: $STOP_AFTER_STEP"
echo "STRS_TO_SKIP: $STRS_TO_SKIP"
echo "VCF_SLOGAN: $VCF_SLOGAN"

# check if all required arguments are provided 

# Pseudobulk h5ad files in GENE_EXPR_INPUT_DIR
# save _expression_matrix_ds.csv and _composition_matrix_ds.csv to GENE_EXPR_INPUT_DIR/PB_OUTPUT_SUBDIR

if [ $START_AT_STEP -le 1 ]; then
    echo "************************************* STEP 1 *************************************"
    echo "************************************* PSEUDOBULK *********************************"
    python $SCRIPT_DIR/01-make-eqtl-pseudobulk.py \
        --input-dir "$GENE_EXPR_INPUT_DIR" \
        --output-dir "$PB_OUTPUT_DIR" \
        --min-num-cells "$MIN_NUM_CELLS" \
        --gene-log-expr-threshold "$GENE_LOG_EXPR_THRESHOLD" \
        --ct-id "$CT_ID" \
        --sample-id "$SAMPLE_ID" \
        --strs-to-skip "$STRS_TO_SKIP" && {
            echo "STEP 1 SUCCESSFULLY created pseudobulk expression and composition matrices."
        } || {
            echo "STEP 1 FAILED to create pseudobulk expression and composition matrices."
            exit 1
        }
else
    echo "************************************* SKIPPING STEP 1 ****************************"
fi

# only run if start-at-step is <= 2 and stop-after-step is >= 2
if [ $START_AT_STEP -le 2 ] && [ $STOP_AFTER_STEP -ge 2 ]; then

    echo "************************************* STEP 2 *************************************"
    echo "************************************* FORMAT OSCA INPUTS *************************"
    # format the OSCA inputs
    Rscript $SCRIPT_DIR/02-run-osca-formatting-scanpy.R \
        --expression-dir="$PB_OUTPUT_DIR" \
        --output-dir="$OSCA_INPUT_DIR" \
        --vcf-slogan="$VCF_SLOGAN" \
        --participants="$PARTICIPANT_FNAME" && {
            echo "STEP 2 SUCCESSFULLY formatted OSCA inputs."
        } || {
            echo "STEP 2 FAILED to format OSCA inputs."
            exit 1
        }
else
    echo "************************************* SKIPPING STEP 2 ****************************"
fi

if [ $START_AT_STEP -le 3 ] && [ $STOP_AFTER_STEP -ge 3 ]; then
    echo "************************************* STEP 3 *************************************"
    echo "************************************* RUN OSCA ***********************************"
    echo "Running OSCA eQTL analysis with the following parameters:"
    echo "efile: $OSCA_INPUT_DIR/Phenotype_@_osca.txt"
    echo "befile_prefix: $OSCA_INPUT_DIR/befile_@"
    echo "bfile: $OSCA_INPUT_DIR/$VCF_SLOGAN"
    echo "update_opi: $OSCA_INPUT_DIR/Upprobe_@.opi"
    echo "covar_file: $OSCA_INPUT_DIR/cov1_@.txt"
    echo "qcovar_file: $OSCA_INPUT_DIR/cov2_@_reduced.txt"
    echo "cores: 6"
    echo "final_output: $OSCA_OUTPUT_DIR/eqtl_@.tsv"

    # run OSCA on all CCs simultaneously
    cat $CC_FILE | \
        xargs -I @ echo $SCRIPT_DIR/03-build-eqtl.sh \
            "$OSCA_INPUT_DIR/Phenotype_@_osca.txt" \
            "$OSCA_INPUT_DIR/befile_@" \
            "$OSCA_INPUT_DIR/$VCF_SLOGAN" \
            "$OSCA_INPUT_DIR/Upprobe_@.opi" \
            "$OSCA_INPUT_DIR/cov1_@.txt" \
            "$OSCA_INPUT_DIR/cov2_@_reduced.txt" \
            6 \
            "$OSCA_OUTPUT_DIR/eqtl_@.tsv" | \
        parallel -j 0 --tmpdir /mnt/accessory/tmp && {
            echo "STEP 3 SUCCESSFULLY ran OSCA eQTL pipeline."
        } || {
            echo "STEP 3 FAILED to run OSCA eQTL pipeline."
            exit 1
        }

    mv "$OSCA_INPUT_DIR/befile_*" "$OSCA_OUTPUT_DIR"
    mv "$OSCA_INPUT_DIR/eqtl*log" "$OSCA_OUTPUT_DIR"
    mv "$OSCA_INPUT_DIR/osca.log" "$OSCA_OUTPUT_DIR"
    mv "$OSCA_INPUT_DIR/*tsv" "$OSCA_OUTPUT_DIR"
    rm "$OSCA_INPUT_DIR/tempeqtl_*"
else
    echo "************************************* SKIPPING STEP 3 ****************************"
fi

if [ $START_AT_STEP -le 4 ] && [ $STOP_AFTER_STEP -ge 4 ]; then
    echo "************************************* STEP 4 *************************************"
    echo "************************************* GS COPY TSVS *******************************"
    # copy tsvs so they can be deleted after next step
    gcloud storage cp "$OSCA_OUTPUT_DIR/*tsv" "$GOOGLE_BUCKET/" && {
        echo "STEP 4 SUCCESSFULLY copied OSCA .tsv outputs to Google Cloud Storage."
    } || {
        echo "STEP 4 FAILED to copy OSCA .tsv outputs to Google Cloud Storage."
        exit 1
    }
else
    echo "************************************* SKIPPING STEP 4 ****************************"
fi

if [ $START_AT_STEP -le 5 ] && [ $STOP_AFTER_STEP -ge 5 ]; then
    echo "************************************* STEP 5 *************************************"
    echo "************************************* SAVE OSCA RDS, MANHATTAN PLOTS *************"
    # process and plot OSCA outputs in parallel. Makes plots and saves huge tsvs as rds
    cat $CC_FILE | \
        xargs -I @ echo Rscript $SCRIPT_DIR/05-process-and-plot-osca-tsv.R \
        --path="$OSCA_OUTPUT_DIR/eqtl_@.tsv" | parallel -j 0 --tmpdir /mnt/accessory/tmp && {
            echo "STEP 5 SUCCESSFULLY processed and plotted OSCA outputs."
        } || {
            echo "STEP 5 FAILED to process and plot OSCA outputs."
            exit 1
        }

    rm $OSCA_OUTPUT_DIR/*tsv
else
    echo "************************************* SKIPPING STEP 5 ****************************"
fi

if [ $START_AT_STEP -le 6 ] && [ $STOP_AFTER_STEP -ge 6 ]; then
    echo "************************************* STEP 6 *************************************"
    echo "************************************* SUBSET TO COMMON SNP-PROBES ****************"
    Rscript $SCRIPT_DIR/06-get-common-snp-probes-osca-rds.R \
        --base="$OSCA_OUTPUT_DIR" && {
            echo "STEP 6 SUCCESSFULLY created common SNP-by-Gene rds."
        } || {
            echo "STEP 6 FAILED to create common SNP-by-Gene rds."
            exit 1
        }
else 
    echo "************************************* SKIPPING STEP 6 ****************************"
fi

if [ $START_AT_STEP -le 7 ] && [ $STOP_AFTER_STEP -ge 7 ]; then
    echo "************************************* STEP 7 *************************************"
    echo "************************************* MASH ***************************************"
    Rscript $SCRIPT_DIR/07-run-mashr.R \
        --path="$OSCA_OUTPUT_DIR/eqtl_present_in_all.rds" \
        --padj-thresh="$MASH_PADJ_THRESH" \
        --num-random="$MASH_NUM_RANDOM" \
        --eps="$MASH_EPS" && {
            echo "STEP 7 SUCCESSFULLY ran mashr."
        } || {
            echo "STEP 7 FAILED to run mashr."
            exit 1
        }
else 
    echo "************************************* SKIPPING STEP 7 ****************************"
fi

if [ $START_AT_STEP -le 8 ] && [ $STOP_AFTER_STEP -ge 8 ]; then
    echo "************************************* STEP 8 *************************************"
    echo "************************************* GS COPY OSCA INPUTS AND OUTPUTS ************"
    gcloud storage cp -r $OSCA_INPUT_DIR $GOOGLE_BUCKET/ && {
        echo "STEP 8 SUCCESSFULLY copied OSCA inputs to Google Cloud Storage."
    } || {
        echo "STEP 8 FAILED to copy OSCA inputs to Google Cloud Storage."
        exit 1
    }
else
    echo "************************************* SKIPPING STEP 8 ****************************"
fi

echo "************************************* PIPELINE COMPLETE **************************"