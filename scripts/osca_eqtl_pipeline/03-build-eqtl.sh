#!/usr/bin/env bash
# 03-build-eqtl.sh — Run OSCA eQTL analysis for a single cell class.
#
# Positional arguments:
#   $1  efile          — expression phenotype file (OSCA txt format; from step 2 Phenotype_{cc}_osca.txt)
#   $2  befile_prefix  — output prefix for OSCA BOD binary files (e.g. befile_{cc})
#   $3  bfile          — PLINK bfile prefix; expects {bfile}.bed, {bfile}.bim, {bfile}.fam
#   $4  update_opi     — updated probe info file (from step 2 Upprobe_{cc}.opi)
#   $5  covar_file     — categorical covariate file (from step 2 cov1_{cc}.txt)
#   $6  qcovar_file    — quantitative covariate file (from step 2 cov2_{cc}_reduced.txt)
#   $7  cores          — number of threads for OSCA
#   $8  final_output   — output prefix for final TSV (e.g. {osca_output_dir}/eqtl_{cc})
#
# OSCA sub-steps:
#   1. --make-bod         Encode expression phenotype as binary BOD format
#   2. --update-opi       Update probe info (gene positions) in BOD
#   3. --eqtl --cis       Cis-eQTL linear model test within CIS_WIND kb
#   4. --query 1          Export all cis SNP-probe pairs to TSV
#
# Outputs:
#   {output_dir}/eqtl_{cc}.tsv           — raw OSCA cis-eQTL results
#   {output_dir}/befile_{cc}.*           — OSCA BOD binary expression files
#   {output_dir}/eqtl_{cc}_osca_*.log   — per-step OSCA logs

# Check if correct number of arguments is provided
if [[ $# -lt 8 ]]; then
    echo "Usage: $0 <efile> <befile_prefix> <bfile> <update_opi> <covar_file> <qcovar_file> <cores> <final_output>"
    exit 1
fi

# Parse input arguments
efile=$1 # constructed from OSCA_formatting-scanpy.ipynb  ( has _*osca.txt) --expression file formatted
befile_prefix=$2 #just to write a prefix for outputs
bfile=$3 # your plink bfiles. should have <bfile>.fam, <bfile>.bed, <bfile>.bim
update_opi=$4 # constructed from OSCA_formatting-scanpy.ipynb . ( has  *opi ) . Has updated probe info
covar_file=$5 # cov1*.txt. Categorical covariate file from OSCA_formatting-scanpy.ipynb .
qcovar_file=$6 # cov2*.txt. Quantitative covariate file from OSCA_formatting-scanpy.ipynb .
cores=$7 # number of cores to use
final_output=$8 # output prefix

# CIS_WIND unit is Kb. 1000 Kb = 1 Mb (standard cis-eQTL window).
# OSCA's default is 2000 Kb (2 Mb). We use 1 Mb here.
CIS_WIND=1000

# Set output directory from final_output path
output_dir=$(dirname "$final_output")
mkdir -p "$output_dir"

# Log file names with prefix from final output. Remove the file extension
log_prefix="$output_dir/$(basename "${final_output%.*}")_"
output_log="${log_prefix}osca_output.log"
error_log="${log_prefix}osca_error.log"
progress_log="${log_prefix}osca_progress.log"

# Validate required input files exist before any OSCA call
for f in "$efile" "$bfile.bed" "$bfile.bim" "$bfile.fam" "$update_opi" "$covar_file" "$qcovar_file"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: required file not found: $f"
        exit 1
    fi
done

# Step 1: Create BOD file
echo "Step 1: Running gene-expression --make-bod" >> ${progress_log}
osca --efile "$efile" --gene-expression --make-bod --out "$befile_prefix" >> ${progress_log} 2>&1 && {
    echo "Step 1, Creating BOD File: COMPLETED" >> ${progress_log}
} || {
    echo "Step 1, Creating BOD File: FAILED" >> ${progress_log}
    exit 1
}

# Step 2: Update OPI with user-specified input
echo "Step 2: Updating OPI using: $update_opi" >> ${progress_log}
osca --befile "$befile_prefix" --update-opi "$update_opi" >> ${progress_log} 2>&1 && {
    echo "Step 2, Updating OPI: COMPLETED" >> ${progress_log}
} || {
    echo "Step 2, Updating OPI: FAILED" >> ${progress_log}
    exit 1
}

# Step 3: Run eQTL analysis
temp_eqtl_output="$output_dir/tempeqtl_$(basename "$final_output")"
echo "Step 3: Running eQTL analysis" >> ${progress_log}
osca --eqtl --bfile "$bfile" --befile "$befile_prefix" --cis --cis-wind "$CIS_WIND" \
    --covar "$covar_file" --qcovar "$qcovar_file" --to-smr --thread-num "$cores" \
    --out "$temp_eqtl_output" >> ${progress_log} 2>&1 && {
    echo "Step 3, Running eQTL analysis: COMPLETED" >> ${progress_log}
} || {
    echo "Step 3, Running eQTL analysis: FAILED" >> ${progress_log}
    exit 1
}

# Step 4: Query eQTL summary with user-specified output
echo "Step 4: Querying eQTL summary with output: $final_output" >> ${progress_log}
osca --beqtl-summary "$temp_eqtl_output" --query 1 --out "$final_output" >> ${progress_log} 2>&1 && {
    echo "Step 4, Querying eQTL summary: COMPLETED" >> ${progress_log}
} || {
    echo "Step 4, Querying eQTL summary: FAILED" >> ${progress_log}
    exit 1
}

# Completion
echo "OSCA job completed at $(date)" >> ${progress_log}
