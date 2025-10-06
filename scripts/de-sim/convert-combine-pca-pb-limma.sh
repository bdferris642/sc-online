#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ----------------- args -----------------
START_AT_STEP=1
STOP_AFTER_STEP=999
geneset=""

pos=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-at-step)
      [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || { echo "Error: --start-at-step requires an integer."; exit 2; }
      START_AT_STEP="$2"; shift 2;;
    --stop-after-step)
      [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || { echo "Error: --stop-after-step requires an integer."; exit 2; }
      STOP_AFTER_STEP="$2"; shift 2;;
    --geneset)
      [[ $# -ge 2 ]] || { echo "Error: --geneset requires a value."; exit 2; }
      geneset="$2"; shift 2;;
    --) shift; break;;
    -*) echo "Unknown parameter: $1"; exit 1;;
    *) pos+=("$1"); shift;;
  esac
done
if (( ${#pos[@]} >= 1 )); then START_AT_STEP="${pos[0]}"; fi
if (( ${#pos[@]} >= 2 )); then STOP_AFTER_STEP="${pos[1]}"; fi
if (( ${#pos[@]} >= 3 )); then geneset="${pos[2]}"; fi
[[ -n "${geneset}" ]] || { echo "Error: geneset is required (positional or --geneset)."; exit 2; }

# ----------------- paths -----------------
slogan="zonated_objs_combined_with_md"
GS_LOC="gs://macosko_data/ferris/bican/simulations/spn_type/split_001"

split_dir="/mnt/accessory/seq_data/bican_str/zonated_objs/simulations/${slogan}/split_001"
geneset_dir="${split_dir}/${geneset}"
ctr_qs="${split_dir}/${slogan}__control.qs"

mkdir -p "${geneset_dir}"

echo "Starting processing for ${slogan}..."
echo "Start at step: ${START_AT_STEP}"
echo "Stop after step: ${STOP_AFTER_STEP}"
echo "Split directory: ${split_dir}"
echo "Geneset directory: ${geneset_dir}"
echo "Control qs path: ${ctr_qs}"

# ----------------- step 1 -----------------
if (( START_AT_STEP <= 1 )); then
  echo "Step 1: Convert .h5ad -> Seurat/qs …"
  find "$geneset_dir" -maxdepth 1 -type f -name '*.h5ad' -print0 \
    | xargs -0 -I @ echo Rscript "$HOME/sc-online/scripts/convert-anndata-h5ad-to-seurat-qs.R" --path=@ | bash

  # Build cases.txt (full paths to *__case__*.qs)
  find "$geneset_dir" -maxdepth 1 -type f -name '*__case__*.qs' -print > "${geneset_dir}/cases.txt"

  if (( STOP_AFTER_STEP == 1 )); then
    echo "Stopping after step 1 as requested."
    exit 0
  fi
fi

# ----------------- step 3 -----------------
if (( START_AT_STEP <= 3 )); then
  echo "Step 3: Add control PCs to each case sobj …"
  if [[ ! -s "${geneset_dir}/cases.txt" ]]; then
    echo "cases.txt missing or empty; rebuilding."
    find "$geneset_dir" -maxdepth 1 -type f -name '*__case__*.qs' -print > "${geneset_dir}/cases.txt"
  fi

  # simple sequential loop (can parallelize similarly if desired)
  while IFS= read -r case_qs; do
    [[ -z "$case_qs" ]] && continue
    Rscript "$HOME/sc-online/scripts/add-ctr-pcs.R" \
      --ctr-path="${ctr_qs}" \
      --sobj-path="${case_qs}"
  done < "${geneset_dir}/cases.txt"

  if (( STOP_AFTER_STEP == 3 )); then
    echo "Stopping after step 3 as requested."
    exit 0
  fi
fi

# ----------------- step 4 -----------------
if (( START_AT_STEP <= 4 )); then
  echo "Step 4: Combine control with each case qs …"
  find "$geneset_dir" -maxdepth 1 -type f -name '*__case__*.qs' -print0 \
    | xargs -0 -I @ echo Rscript "$HOME/sc-online/scripts/combine-sobj-qs.R" \
      --p1="${ctr_qs}" --p2=@ --out-dir="${geneset_dir}" | parallel -j 2

  if (( STOP_AFTER_STEP == 4 )); then
    echo "Stopping after step 4 as requested."
    exit 0
  fi
fi

# ----------------- step 5 -----------------
if (( START_AT_STEP <= 5 )); then
  echo "Step 5: Make pseudobulk objects …"
  mkdir -p "${geneset_dir}/pseudobulk"

  find "$geneset_dir" -maxdepth 1 -type f -name '*__combined__*.qs' -print0 \
    | xargs -0 -I @ echo Rscript "$HOME/sc-online/scripts/make-pseudobulks.R" \
        --path=@ --contrast-col=case_control --sample-col=donor_id \
        --grouping-cols=donor_id,library,case_control,Biobank,Age.at.Death,Sex,Race,Ethnicity,PMI \
        --cols-to-weighted-avg=pct_mt,pct_intronic,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
        --n-svs=0 | parallel -j 3

  # Consolidate .qs into the pseudobulk dir (if the R script wrote elsewhere)
  find "${geneset_dir}/pseudobulk" -type f -name '*.qs' -print0 \
    | xargs -0 -I{} mv -t "${geneset_dir}/pseudobulk" "{}" || true

  if (( STOP_AFTER_STEP == 5 )); then
    echo "Stopping after step 5 as requested."
    exit 0
  fi
fi

# ----------------- step 6 -----------------
if (( START_AT_STEP <= 6 )); then
  echo "Step 6: Run limma (three covar sets) …"
  pb_dir="${geneset_dir}/pseudobulk"
  mkdir -p "${split_dir}/de_results/${geneset}"

  # 6a) vanilla
  find "$pb_dir" -maxdepth 1 -type f -name '*.qs' -print0 \
    | xargs -0 -I @ echo Rscript "$HOME/sc-online/scripts/run-limma-v2.R" \
        --path=@ --contrast-col=case_control --rand-var=donor_id --num-threads=1 \
        --cov-list=case_control,Age.at.Death,Sex,PMI,pct_mt,pct_intronic,log10_nUMI \
        --suffix=vanilla --calc-purity=F | parallel -j 3

  # 6b) + 4 PCs
  find "$pb_dir" -maxdepth 1 -type f -name '*.qs' -print0 \
    | xargs -0 -I @ echo Rscript "$HOME/sc-online/scripts/run-limma-v2.R" \
        --path=@ --contrast-col=case_control --rand-var=donor_id --num-threads=1 \
        --cov-list=case_control,Age.at.Death,Sex,PMI,pct_mt,pct_intronic,log10_nUMI,PC_1,PC_2,PC_3,PC_4 \
        --suffix=4pcs --calc-purity=F | parallel -j 3

  # 6c) + 20 PCs
  find "$pb_dir" -maxdepth 1 -type f -name '*.qs' -print0 \
    | xargs -0 -I @ echo Rscript "$HOME/sc-online/scripts/run-limma-v2.R" \
        --path=@ --contrast-col=case_control --rand-var=donor_id --num-threads=1 \
        --cov-list=case_control,Age.at.Death,Sex,PMI,pct_mt,pct_intronic,log10_nUMI,PC_1,PC_2,PC_3,PC_4,PC_5,PC_6,PC_7,PC_8,PC_9,PC_10,PC_11,PC_12,PC_13,PC_14,PC_15,PC_16,PC_17,PC_18,PC_19,PC_20 \
        --suffix=20pcs --calc-purity=F | parallel -j 3

  # Ship to GCS
  gcloud storage mv "${geneset_dir}" "${GS_LOC}/"

  if (( STOP_AFTER_STEP == 6 )); then
    echo "Stopping after step 6 as requested."
    exit 0
  fi
fi
