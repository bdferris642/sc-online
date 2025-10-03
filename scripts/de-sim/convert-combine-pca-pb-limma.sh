#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

START_AT_STEP=$1
STOP_AFTER_STEP=$2
geneset=$3

DRY_RUN=0  # Set to 1 to enable dry run mode (echo commands instead of executing)

# ---- arg parsing ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-at-step)
      [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || { echo "Error: --start-at-step requires an integer."; exit 2; }
      START_AT_STEP="$2"; shift 2 ;;
    --stop-after-step)
      [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || { echo "Error: --stop-after-step requires an integer."; exit 2; }
      STOP_AFTER_STEP="$2"; shift 2 ;;
    *)
      echo "Unknown parameter: $1"; exit 1 ;;
  esac
done

slogan="zonated_objs_combined_with_md"
GS_LOC="gs://macosko_data/ferris/bican/simulations/spn_type/split_001"
split_dir="/mnt/accessory/seq_data/bican_str/zonated_objs/simulations/${slogan}/${geneset}"
geneset_dir="${split_dir}/${geneset}"
ctr_qs="${split_dir}/${slogan}__control.qs"


echo "Starting processing for ${slogan}..."
echo "Start at step: ${START_AT_STEP}"
echo "Stop after step: ${STOP_AFTER_STEP}"
echo "Split directory: ${split_dir}"
echo "Geneset directory: ${geneset_dir}"
echo "Control qs path: ${ctr_qs}"

# ---- step 1 ----
if (( START_AT_STEP <= 1 )); then
  echo "Converting h5ad to seurat and qs objects for ${slogan}..."
  # Run convert script for each .h5ad in split_dir
  # - Using GNU parallel directly on file list
  find "$geneset_dir" -type f -name *.h5ad -exec realpath {} + | \
    xargs -I {} echo Rscript ~/sc-online/scripts/convert-anndata-h5ad-to-seurat-qs.R --path={} | bash

  # Build cases.txt: full paths to *_case_*.qs files
  rm -f "${geneset_dir}/cases.txt"
  touch "${geneset_dir}/cases.txt"
  find "$geneset_dir" -type f -name *.qs -exec realpath {} + | \
    >> "${geneset_dir}/cases.txt"

  if (( STOP_AFTER_STEP == 1 )); then
    echo "Stopping after step 1 as requested."
    exit 0
  fi
fi

# # ---- step 2 ----
# if (( START_AT_STEP <= 2 )); then
#   Rscript ~/sc-online/scripts/normalizeScalePcaClusterUmap.R \
#     --path="${ctr_qs}" \
#     --var-feature-subset-col=donor_id \
#     --resolutions=0.2,0.5
#   if (( STOP_AFTER_STEP == 2 )); then
#     echo "Stopping after step 2 as requested."
#     exit 0
#   fi
# fi 

# ---- step 3 ----
if (( START_AT_STEP <= 3 )); then

  cat "${geneset_dir}/cases.txt" | \
  xargs -I @ echo \
    Rscript ~/sc-online/scripts/add-ctr-pcs.R \
      --ctr-path="${ctr_qs}" \
      --sobj-path-list=@ | bash

  if (( STOP_AFTER_STEP == 3 )); then
    echo "Stopping after step 3 as requested."
    exit 0
  fi
fi

# ---- step 4 ----
if (( START_AT_STEP <= 4 )); then
  echo "Combining seurat objects for ${slogan}..."
  # Combine control with each case .qs
  find "$geneset_dir" -maxdepth 1 -type f -name '*__case__*' -name '*.qs' \
  -exec realpath -z {} + \
  | xargs -0 -I{} Rscript ~/sc-online/scripts/combine-sobj-qs.R --p1="$ctr_qs" --p2="{}" | bash
  ls "$geneset_dir" | grep __case__ | grep qs | xargs -I @ echo Rscript ~/sc-online/scripts/combine-sobj-qs.R --p1=$ctr_qs --p2=$(readlink -f @) | bash
  if (( STOP_AFTER_STEP == 4 )); then
    echo "Stopping after step 4 as requested."
    exit 0
  fi
fi

# ---- step 5 ----
if (( START_AT_STEP <= 5 )); then
  echo "Making PseudoBulk objects for ${slogan}..."
    ls "$geneset_dir" | grep '__combined__' | grep qs | \
        xargs -I {} echo Rscript ~/sc-online/scripts/make-pseudobulks.R \
        --path="${geneset_dir}"/{} --contrast-col=case_control --sample-col=donor_id  \
        --grouping-cols=donor_id,library,case_control,Biobank,Age.at.Death,Sex,Race,Ethnicity,PMI \
        --cols-to-weighted-avg=pct_mt,pct_intronic,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
        --n-svs=0 | \
        parallel -j 3 && find "${geneset_dir}/pseudobulk" | grep qs$ | xargs -I @ mv @ ${geneset_dir}/pseudobulk
    if (( STOP_AFTER_STEP == 5 )); then
        echo "Stopping after step 5 as requested."
        exit 0
    fi
fi

# ---- step 6 ----
if (( START_AT_STEP <= 6 )); then
  ls "${geneset_dir}/pseudobulk" | grep -E '\.qs$' | \
    xargs -I @ echo Rscript ~/sc-online/scripts/run-limma-v2.R \
      --path=$(readlink -f @) --contrast-col=case_control --rand-var=donor_id --num-threads=1 \
      --cov-list=case_control,Age.at.Death,Sex,PMI,pct_mt,pct_intronic,log10_nUMI \
      --suffix=vanilla --calc-purity=F | parallel -j 4

  ls "${geneset_dir}/pseudobulk" | grep -E '\.qs$' | \
    xargs -I @ echo Rscript ~/sc-online/scripts/run-limma-v2.R \
      --path=$(readlink -f @) --contrast-col=case_control --rand-var=donor_id --num-threads=1 \
      --cov-list=case_control,Age.at.Death,Sex,PMI,pct_mt,pct_intronic,log10_nUMI,PC_1,PC_2,PC_3,PC_4 \
      --suffix=4pcs --calc-purity=F | parallel -j 4

  ls "${geneset_dir}/pseudobulk" | grep -E '\.qs$' | \
    xargs -I @ echo Rscript ~/sc-online/scripts/run-limma-v2.R \
      --path=$(readlink -f @) --contrast-col=case_control --rand-var=donor_id --num-threads=1 \
      --cov-list=case_control,Age.at.Death,Sex,PMI,pct_mt,pct_intronic,log10_nUMI,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
      --suffix=20pcs --calc-purity=F | parallel -j 4

  mkdir -p "${split_dir}/de_results/${geneset}"
  find "${geneset_dir}/pseudobulk" -name *case_control.csv | xargs -I @ mv @ "${split_dir}/de_results/${geneset}"

  gcloud storage mv ${geneset_dir} $GS_LOC

  if (( STOP_AFTER_STEP == 6 )); then
    echo "Stopping after step 6 as requested."
    exit 0
  fi
fi