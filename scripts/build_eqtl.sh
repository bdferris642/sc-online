# Check if correct number of arguments is provided
if [[ $# -lt 7 ]]; then
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



# Log file names with prefix from final output. Remove the file extension
log_prefix="${final_output%.*}_"
output_log="${log_prefix}osca_output.log"
error_log="${log_prefix}osca_error.log"
progress_log="${log_prefix}osca_progress.log"

# Step 1: Create BOD file
echo "Step 1: Running gene-expression --make-bod" >> ${progress_log}
osca --efile $efile --gene-expression --make-bod --out $befile_prefix >> ${progress_log} 2>&1 && {
    echo "Step 1, Creating BOD File: COMPLETED" >> ${progress_log}
} || {
    echo "Step 1, Creating BOD File: FAILED" >> ${progress_log}
    exit 1
}

# Step 2: Update OPI with user-specified input
echo "Step 2: Updating OPI using: $update_opi" >> ${progress_log}
osca --befile $befile_prefix --update-opi $update_opi >> ${progress_log} 2>&1 && {
    echo "Step 2, Updating OPI: COMPLETED" >> ${progress_log}
} || {
    echo "Step 2, Updating OPI: FAILED" >> ${progress_log}
    exit 1
}

# Step 3: Run eQTL analysis
echo "Step 3: Running eQTL analysis" >> ${progress_log}
osca --eqtl --bfile $bfile --befile $befile_prefix --cis --cis-wind 2000 \
    --covar $covar_file --qcovar $qcovar_file --to-smr --thread-num $cores \
    --out tempeqtl_$final_output >> ${progress_log} 2>&1 && {
    echo "Step 3, Running eQTL analysis: COMPLETED" >> ${progress_log}
} || {
    echo "Step 3, Running eQTL analysis: FAILED" >> ${progress_log}
    exit 1
}


# Step 4: Query eQTL summary with user-specified output
echo "Step 4: Querying eQTL summary with output: $final_output" >> ${progress_log}
osca --beqtl-summary tempeqtl_$final_output  --query 1 --out $final_output >> ${progress_log} 2>&1 && {
    echo "Step 4, Querying eQTL summary: COMPLETED" >> ${progress_log}
} || {
    echo "Step 4, Querying eQTL summary: FAILED" >> ${progress_log}
    exit 1
}


# Completion
echo "OSCA job completed at \$(date)" >> ${progress_log}
