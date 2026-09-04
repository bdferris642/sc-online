#!/bin/bash
# 00-run-osca-eqtl-pipeline.sh — Orchestrate the full OSCA snRNA-seq cis-eQTL pipeline.
#
# Usage:
#   ./00-run-osca-eqtl-pipeline.sh [OPTIONS]
#
# Required arguments:
#   --input-files FILE        Newline-delimited text file listing absolute paths to h5ad files (step 1)
#   --osca-input-dir DIR      Directory for OSCA inputs and outputs; must contain PLINK bfiles
#   --pb-output-dir DIR       Directory for pseudobulk output CSVs (steps 1 and 2)
#   --pipeline-slogan STR     Short label for this run; used in output paths and GCS bucket
#   --vcf-slogan STR          Prefix of .bed/.bim/.fam and _pca.eigenvec in --osca-input-dir
#   --participants FILE        Newline-delimited list of VCF participant IDs to include
#                             Required unless --id-map is provided, in which case participant IDs
#                             are derived from the vcf_sample_id column of the id-map.
#
# Optional arguments (defaults shown):
#   --ct-id STR               adata.obs column for cell type labels          [cell_class]
#   --gene-log-expr-threshold FLOAT  Min mean log-expression to retain a gene [0.01]
#   --h5ad-cat-covar STR      Space-delimited categorical covariate column names from adata.obs
#                             Used for cov1 (as factors) and SVA formula
#                             [default: "sex case_control study brain_bank dapi_nurr"]
#   --h5ad-quant-covar STR    Space-delimited quantitative covariate column names from adata.obs
#                             Used for cov2 (min-max scaled) and SVA formula  [default: "age pmi"]
#   --id-map FILE             CSV with nucseq_participant_id + vcf_sample_id; remaps IDs (step 1)
#   --mash-eps FLOAT          Floor for SE and beta values passed to mashr    [1e-6]
#   --mash-num-random INT     Random SNP-probe pairs for mashr null correlation [100000]
#   --mash-padj-thresh FLOAT  BH adjusted p-value threshold for mashr strong tests [0.01]
#   --min-num-cells INT       Min cells per sample for pseudobulk inclusion   [10]
#   --sample-id STR           adata.obs column for sample/participant IDs      [participant_id]
#   --start-at-step INT       Resume from this step (inclusive)               [1]
#   --stop-after-step INT     Stop after this step (inclusive)                [1000]
#   --strs-to-skip STR        Comma-separated substrings; h5ad files matching any are skipped
#
# Pipeline steps:
#   1  h5ad → pseudobulk CSVs (ID remap → sum aggregation → CPM → log1p → scale)
#   2  pseudobulk CSVs → OSCA input files (SVA + PCA + covariate assembly)
#   3  OSCA cis-eQTL analysis (parallelized per cell class via GNU parallel)
#   4  (optional) Back up raw .tsv results to Google Cloud Storage
#   5  Two-stage FDR; save RDS; generate Manhattan and p-value histogram plots
#   6  Intersect tested SNP-probe pairs across all cell classes
#   7  mashr effect-size sharing across cell classes
#   8  (optional) Archive all inputs and outputs to Google Cloud Storage
#
# Key outputs  (OUTPUT_ROOT = {osca-input-dir}/eqtl_final_outs/{pipeline-slogan}):
#   Step 1 → {pb-output-dir}/{cell_class}_expression_matrix_ds.csv
#             {pb-output-dir}/{cell_class}_composition_matrix_ds.csv
#   Step 2 → {osca-input-dir}/Phenotype_{cc}_osca.txt       (OSCA expression phenotype)
#             {osca-input-dir}/Upprobe_{cc}.opi              (probe annotation)
#             {osca-input-dir}/cov1_{cc}.txt                 (categorical covariates)
#             {osca-input-dir}/cov2_{cc}_reduced.txt         (quantitative covariates)
#   Step 3 → {OUTPUT_ROOT}/eqtl_{cc}.tsv                    (raw OSCA output; deleted after step 5)
#             {OUTPUT_ROOT}/befile_{cc}.*                    (OSCA BOD binary files)
#             {OUTPUT_ROOT}/eqtl_{cc}_osca_{output,error,progress}.log
#   Step 5 → {OUTPUT_ROOT}/eqtl_{cc}.rds                    (full results + FDR columns)
#             {OUTPUT_ROOT}/eqtl_{cc}_sig.rds                (significant results only)
#             {OUTPUT_ROOT}/plots/eqtl_{cc}_manhattan.png
#             {OUTPUT_ROOT}/plots/eqtl_{cc}_{min_p_gene,padj_gene,pval,padj_snp}_hist.png
#   Step 6 → {OUTPUT_ROOT}/eqtl_present_in_all.rds
#             {OUTPUT_ROOT}/eqtl_present_in_all_sig_in_one.rds
#   Step 7 → {OUTPUT_ROOT}/eqtl_present_in_all__mash_results_{sig,random}.rds
#             {OUTPUT_ROOT}/eqtl_present_in_all__mash_data_{sig,random}.rds
#
# Logs: /mnt/accessory/analysis/eqtl/logs/osca_eqtl_pipeline_{pipeline-slogan}_{datetime}.log

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"
ENV_PREFIX="${SANDBOX}/micromamba_root/envs/osca-venv"
ENV_BIN="${ENV_PREFIX}/bin"

if [ ! -d "${ENV_BIN}" ]; then
    echo "ERROR: osca-venv not found at ${ENV_BIN}."
    echo "       Run: bash ${SCRIPT_DIR}/setup.sh"
    exit 1
fi

# Prepend env bin to PATH so parallel subshells find python, Rscript, osca, parallel
export PATH="${ENV_BIN}:${PATH}"
PYTHON="${ENV_BIN}/python"
RSCRIPT="${ENV_BIN}/Rscript"

# Default values
CT_ID="cell_class"
GENE_EXPR_INPUT_FILES=""
GENE_LOG_EXPR_THRESHOLD=0.01
H5AD_CAT_COVARS="sex case_control study brain_bank dapi_nurr"
H5AD_QUANT_COVARS="age pmi"
ID_MAP_PATH=""
MASH_EPS=1e-6
MASH_NUM_RANDOM=100000
MASH_PADJ_THRESH=0.01
MIN_NUM_CELLS=10
SAMPLE_ID="participant_id"
STRS_TO_SKIP=""
START_AT_STEP=1
STOP_AFTER_STEP=1000

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --ct-id) CT_ID="$2"; shift ;;
        --gene-log-expr-threshold) GENE_LOG_EXPR_THRESHOLD="$2"; shift ;;
        --pb-output-dir) PB_OUTPUT_DIR="$2"; shift ;;
        --h5ad-cat-covar) H5AD_CAT_COVARS="$2"; shift ;;
        --h5ad-quant-covar) H5AD_QUANT_COVARS="$2"; shift ;;
        --input-files) GENE_EXPR_INPUT_FILES="$2"; shift ;;
        --id-map ) ID_MAP_PATH="$2"; shift ;;
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
: "${GENE_EXPR_INPUT_FILES:?--input-files is required}"
: "${OSCA_INPUT_DIR:?--osca-input-dir is required}"
: "${PB_OUTPUT_DIR:?--pb-output-dir is required}"
: "${PIPELINE_SLOGAN:?--pipeline-slogan is required}"
: "${VCF_SLOGAN:?--vcf-slogan is required}"

if [ -z "$PARTICIPANT_FNAME" ] && [ -z "$ID_MAP_PATH" ]; then
    echo "ERROR: --participants is required unless --id-map is provided"
    exit 1
fi


# include date in log filenames
DATE_TIME=$(date +'%Y-%m-%d_%H-%M-%S')
GOOGLE_BUCKET=gs://macosko_data/ferris/eqtl_outs/$PIPELINE_SLOGAN
LOG_DIR="${OSCA_INPUT_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FNAME="${LOG_DIR}/osca_eqtl_pipeline_${PIPELINE_SLOGAN}_${DATE_TIME}.log"
# redirect stdout and stderr to log file and terminal
exec > >(tee -i "$LOG_FNAME") 2>&1
echo "Logging to: $LOG_FNAME"

# store git branch and commit hash to log file so pipeline can be replicated later
branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
commit_hash=$(git rev-parse HEAD 2>/dev/null || echo "unknown")

echo "Branch: $branch"
echo "Commit: $commit_hash"

OSCA_OUTPUT_DIR=$OSCA_INPUT_DIR/eqtl_final_outs/$PIPELINE_SLOGAN

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

# If --id-map provided but no --participants, derive participant list from id-map vcf_sample_id column
if [ -z "$PARTICIPANT_FNAME" ] && [ -n "$ID_MAP_PATH" ]; then
    if [ ! -f "$ID_MAP_PATH" ]; then
        echo "ERROR: id-map file not found: $ID_MAP_PATH"
        exit 1
    fi
    PARTICIPANT_FNAME=$(mktemp /tmp/participants_from_idmap_XXXXXX.txt)
    "$PYTHON" -c "
import csv
with open('$ID_MAP_PATH') as f:
    for row in csv.DictReader(f):
        print(row['vcf_sample_id'])
" > "$PARTICIPANT_FNAME"
    echo "Derived $(wc -l < "$PARTICIPANT_FNAME") participants from id-map: $ID_MAP_PATH"
elif [ ! -f "$PARTICIPANT_FNAME" ]; then
    echo "ERROR: participants file not found: $PARTICIPANT_FNAME"
    exit 1
fi

# write all function arguments to log file
echo "Running OSCA eQTL pipeline with the following arguments:"
echo "CT_ID: $CT_ID"
echo "GENE_EXPR_INPUT_FILES: $GENE_EXPR_INPUT_FILES"
echo "H5AD_CAT_COVARS: $H5AD_CAT_COVARS"
echo "H5AD_QUANT_COVARS: $H5AD_QUANT_COVARS"
echo "ID_MAP_PATH: $ID_MAP_PATH"
echo "PB_OUTPUT_DIR: $PB_OUTPUT_DIR"
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

if [ $START_AT_STEP -le 1 ]; then
    echo "************************************* STEP 1 *************************************"
    echo "************************************* PSEUDOBULK *********************************"
    "$PYTHON" "$SCRIPT_DIR/01-make-eqtl-pseudobulk.py" \
        --input-files "$GENE_EXPR_INPUT_FILES" \
        --output-dir "$PB_OUTPUT_DIR" \
        --min-num-cells "$MIN_NUM_CELLS" \
        --gene-log-expr-threshold "$GENE_LOG_EXPR_THRESHOLD" \
        --ct-id "$CT_ID" \
        --sample-id "$SAMPLE_ID" \
        --strs-to-skip "$STRS_TO_SKIP" \
        ${ID_MAP_PATH:+--id-map "$ID_MAP_PATH"} && {
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
    # Build step 2 args as an array so space-delimited covariate lists pass safely
    STEP2_ARGS=(
        "--expression-dir=$PB_OUTPUT_DIR"
        "--output-dir=$OSCA_INPUT_DIR"
        "--vcf-slogan=$VCF_SLOGAN"
        "--participants=$PARTICIPANT_FNAME"
    )
    [ -n "$H5AD_CAT_COVARS" ]   && STEP2_ARGS+=("--h5ad-cat-covars=$H5AD_CAT_COVARS")
    [ -n "$H5AD_QUANT_COVARS" ] && STEP2_ARGS+=("--h5ad-quant-covars=$H5AD_QUANT_COVARS")
    # format the OSCA inputs
    "$RSCRIPT" "$SCRIPT_DIR/02-run-osca-formatting-scanpy.R" "${STEP2_ARGS[@]}" && {
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
    # Discover cell classes from Phenotype files written by step 2
    mapfile -t CCS < <(
        ls "$OSCA_INPUT_DIR"/Phenotype_*_osca.txt 2>/dev/null \
        | sed 's|.*/Phenotype_||; s|_osca\.txt$||'
    )
    if [ ${#CCS[@]} -eq 0 ]; then
        echo "ERROR: No Phenotype_*_osca.txt files found in $OSCA_INPUT_DIR. Run step 2 first."
        exit 1
    fi
    echo "Cell classes: ${CCS[*]}"

    # run OSCA on all cell classes simultaneously
    printf '%s\n' "${CCS[@]}" | \
        xargs -I @ echo $SCRIPT_DIR/03-build-eqtl.sh \
            "$OSCA_INPUT_DIR/Phenotype_@_osca.txt" \
            "$OSCA_INPUT_DIR/befile_@" \
            "$OSCA_INPUT_DIR/$VCF_SLOGAN" \
            "$OSCA_INPUT_DIR/Upprobe_@.opi" \
            "$OSCA_INPUT_DIR/cov1_@.txt" \
            "$OSCA_INPUT_DIR/cov2_@_reduced.txt" \
            5 \
            "$OSCA_OUTPUT_DIR/eqtl_@.tsv" | \
        parallel -j 0 --tmpdir "$OSCA_OUTPUT_DIR" && {
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
    # Discover TSVs directly — no cc-file needed
    ls "$OSCA_OUTPUT_DIR"/eqtl_*.tsv 2>/dev/null | \
        xargs -I {} echo "$RSCRIPT" "$SCRIPT_DIR/05-process-and-plot-osca-tsv.R" \
        --path="{}" | parallel -j 0 --tmpdir "$OSCA_OUTPUT_DIR" && {
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
    "$RSCRIPT" "$SCRIPT_DIR/06-get-common-snp-probes-osca-rds.R" \
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
    "$RSCRIPT" "$SCRIPT_DIR/07-run-mashr.R" \
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
