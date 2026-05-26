#!/usr/bin/env bash
set -euo pipefail

# Orchestrates memory-friendly, gene-wise NEBULA runs.

# -----------------------------
# Args
# -----------------------------
PATH_QS=""
ID_COL="donor_id"
COVS=""
OFFSET_COL="nCount_RNA"
N_FOLDS=""
N_CORES=""
SAVE_TMP="0"
SUFFIX=""
MIN_COUNTS_GENE="10"
MIN_FRAC_GENE="0.01"

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") \\
    --path /path/to/object.qs \\
    --covs sex,batch,age \\
    --offset-col nCount_RNA \\
    --n-folds 8 \\
    --n-cores 8 \\
    --min-counts-gene 10 \\
    --min-frac-gene 0.01 \\
    --save-tmp 0|1 \\
    --suffix <suffix>


Notes:
- chunk-sobj-by-genes.R (creates chunks with ALL CELLS and a CHUNK OF GENES, and meta {id, covs, nUMI}).
- run-mast.R on each chunk with in parallel using GNU parallel.
- combine-mast-results.R (row-binds fields and writes combined .csv for each coef).
- optionally removes tmp.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --path) PATH_QS="$2"; shift 2 ;;
    --id-col) ID_COL="$2"; shift 2 ;;
    --covs) COVS="$2"; shift 2 ;;
    --offset-col) OFFSET_COL="$2"; shift 2 ;;
    --n-folds) N_FOLDS="$2"; shift 2 ;;
    --n-cores) N_CORES="$2"; shift 2 ;;
    --save-tmp) SAVE_TMP="$2"; shift 2 ;;
    --suffix ) SUFFIX="$2"; shift 2 ;;
    --min-counts-gene ) MIN_COUNTS_GENE="$2"; shift 2 ;;
    --min-frac-gene ) MIN_FRAC_GENE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done


[[ -z "${PATH_QS}" || -z "${N_FOLDS}" || -z "${N_CORES}" ]] && { echo "ERROR: --path, --n-folds, --n-cores required"; usage; exit 1; }
[[ ! -f "$PATH_QS" ]] && { echo "ERROR: File not found: $PATH_QS"; exit 1; }


SCRIPT_DIR="$(dirname "$(realpath "$0")")" # ~/sc-online/scripts
base="$(basename "$PATH_QS")"
d="$(dirname "$PATH_QS")"

# slogan like .qs -> __mast.qs
slogan="${base/%.qs/__mast.qs}"
tmp_base="${slogan%.qs}"                         # drop .qs
tmp_dir="${d}/${tmp_base}"                      # tmp dir next to input file
chunks_dir="${tmp_dir}/chunks"                   # where chunk .qs files go
tmp_de_dir="${tmp_dir}/de_results/"  # where per-chunk NEBULA outputs go

final_dir="${d}/de_results/mast"

if [[ -n "${SUFFIX:-}" ]]; then
  out_basename="${base/%.qs/__mast_combined__$SUFFIX.qs}"
  suffix_str="--suffix $SUFFIX"
else
  out_basename="${base/%.qs/__mast_combined.qs}"
  suffix_str=""
fi

final_qs="${final_dir}/${out_basename}"
mkdir -p "$final_dir" "$chunks_dir" "$tmp_de_dir"

echo "[INFO] Input: $PATH_QS"
echo "[INFO] Script dir: $SCRIPT_DIR"
echo "[INFO] Base name: $base"
echo "[INFO] n-folds: $N_FOLDS"
echo "[INFO] parallel jobs: $N_CORES"

# comma split covs, put each on a newline
echo "[INFO] covs:"
if [[ -n "${COVS}" ]]; then
  IFS=',' read -ra COVS_ARR <<< "$COVS"
  for cov in "${COVS_ARR[@]}"; do
    echo "       - $cov"
  done
else
  echo "       - <none>"
fi
#echo "[INFO] covs: ${COVS:-<none>}"

echo "[INFO] Tmp dir: $tmp_dir"
echo "[INFO] Chunks dir: $chunks_dir"
echo "[INFO] Tmp de_dir: $tmp_de_dir"

echo "[INFO] final_dir: $final_dir"
echo "[INFO] nebula_out_basename: $out_basename"
echo "[INFO] Final combined results: $final_qs"

echo


# -----------------------------
# Chunk (separate script)
# -----------------------------
Rscript "$(dirname "$0")/chunk-sobj-by-genes.R" \
  --path "$PATH_QS" \
  --covs "$COVS" \
  --n-folds "$N_FOLDS" \
  --outdir "$chunks_dir" \
  --new-offset-colname "nUMI" || {
    echo "[ERROR] chunk-sobj-by-genes.R failed"; exit 1;
  }

echo

# -----------------------------
# Run run-mast.R per chunk (GNU parallel)
#   - Each job uses: -n 1
#   - Add a unique -s <chunk_name> for disambiguation
# -----------------------------
echo "[INFO] Launching chunk jobs with GNU parallel..."

# Prepare command stream for parallel
find "$chunks_dir" -maxdepth 1 -type f -name "chunk_*.qs" -print0 \
    | xargs -0 -I @ echo Rscript "$SCRIPT_DIR/run-mast.R" \
        --path @ \
        --covs "$COVS" \
        --min-frac-gene "$MIN_FRAC_GENE" \
        --min-counts-gene "$MIN_COUNTS_GENE" \
        --n-cores 1 --out-dir "$tmp_de_dir" "$suffix_str" \
    | parallel -j "${N_CORES}"


echo "[INFO] All chunk jobs completed."
echo

# -----------------------------
# Step 5: Combine (separate script)
# -----------------------------

Rscript "$SCRIPT_DIR/combine-mast-results.R" \
  --in-dir "$tmp_de_dir" \
  --out "$final_qs"

echo "[INFO] Final combined results: ${final_qs}"

# -----------------------------
# Step 6: Cleanup
# -----------------------------
if [[ "${SAVE_TMP}" == "0" ]]; then
  echo "[INFO] Removing tmp dir: ${tmp_dir}"
  rm -rf "${tmp_dir}"
else
  echo "[INFO] Leaving tmp dir (save-tmp=1): ${tmp_dir}"
fi

echo "[ALL DONE]"

