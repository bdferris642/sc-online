#!/usr/bin/env bash
set -euo pipefail

# Isolate R from personal startup files and libraries: a ~/.Rprofile that loads packages absent
# from this env (e.g. library(repr)) aborts every R start, and R_LIBS* could shadow the env's packages.
export R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null
unset R_LIBS R_LIBS_USER

# run_nebula_parallel: gene-chunked NEBULA (LN) straight from an h5ad, optionally preceded by the
# control-PC step, with rank purity from the h5ad. Portable and self-activating: every script it
# calls lives in this directory and runs under the run-nebula-parallel-env micromamba env
# (create it with setup.sh). See README.md.

WF_DIR="$(cd "$(dirname "$(realpath "$0")")" && pwd)"
SANDBOX="$(cd "${WF_DIR}/.." && pwd)"
ENV_PREFIX="${SANDBOX}/micromamba_root/envs/run-nebula-parallel-env"
ENV_BIN="${ENV_PREFIX}/bin"
PY="${ENV_BIN}/python"
RS="${ENV_BIN}/Rscript"
[[ -x "$RS" && -x "$PY" ]] || { echo "ERROR: workflow env not found at ${ENV_PREFIX}; run ${WF_DIR}/setup.sh"; exit 1; }
export PATH="${ENV_BIN}:${PATH}"

# -----------------------------
# Args
# -----------------------------
H5AD=""
LAYER="counts"
OUTDIR=""
ID_COL="participant_id"
COVS=""
OFFSET_COL="nCount_RNA"
CONTRAST_COL="case_control"
CASE_VAL="pd"
RUN_WITHOUT_CONTRAST=0
N_FOLDS=""
N_CORES=""
CHUNK_WORKERS=8
IO_WORKERS=24
SAVE_TMP="0"
SUFFIX=""
GET_RANK_PURITY=0
PURITY_SCALE="log-rate"
PURITY_LAYER="norm_log"
DERIVE_REGION=0
COMPUTE_PCS=0
CONTROL_VAL="ctr"
N_HVGS=2000
N_PCS=30
PC_TEST_COVS="age,sex,log10_nUMI,frac_mito,frac_intronic,chemistry"
PC_FDR=0.05
PCA_MODEL_DIR=""
PCS_TO_ADD=()
PC_MAX=20

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") \\
    --path /path/to/object.h5ad --outdir /path/to/outdir \\
    --covs case_control,age,sex,... --contrast-col case_control --case-val pd \\
    --n-folds 40 --n-cores 40 \\
    [--compute-control-pcs 0|1] [--pcs-to-add elbow significant] [--pc-max 20] \\
    [--get-rank-purity 0|1] [--suffix <suffix>] [other options below]

Core:
  --path, --outdir, --n-folds, --n-cores (required)
  --layer counts            CSR raw-count layer (cells x genes); 'X' = adata.X
  --id-col participant_id
  --covs a,b,c              NEBULA covariates (interactions a:b / a*b allowed)
  --offset-col nCount_RNA   nCount_RNA = per-cell total of --layer over all genes; else an obs column
  --contrast-col / --case-val   the tested contrast; MUST be in --covs (see --run-without-contrast)
  --run-without-contrast 0|1    allow a model without the contrast column (rank purity then disallowed)
  --chunk-workers 8 --io-workers 24 --save-tmp 0|1 --suffix S --derive-region-dummies 0|1

Control PCs (step 0):
  --compute-control-pcs 0|1  fit a PCA on control cells (2000 HVGs), project all cells, test each PC
                             case vs control on donor means (OLS), and replace obs control_PC* +
                             uns/control_pca IN PLACE in the h5ad (sidecar backup of obs/uns first)
  --control-val ctr          control level of --contrast-col
  --n-hvgs 2000 --n-pcs 30 --pc-fdr 0.05
  --pc-test-covs ${PC_TEST_COVS}
  --pca-model-dir DIR        default <h5ad dir>/pca_models
  --pcs-to-add elbow [significant] [all]   append PCs (union of the sets) to --covs; one or more
                             words; 'all' = every computed PC (control_PC_01, 02, ...)
  --pc-max 20                only PCs with index <= pc-max are added (with 'all': PCs 01..pc-max)

Rank purity:
  --get-rank-purity 0|1  --purity-scale log-rate|mean  --purity-layer norm_log
USAGE
}

ARGS_ORIG=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --path) H5AD="$2"; shift 2 ;;
    --layer) LAYER="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --id-col) ID_COL="$2"; shift 2 ;;
    --covs) COVS="$2"; shift 2 ;;
    --offset-col) OFFSET_COL="$2"; shift 2 ;;
    --contrast-col) CONTRAST_COL="$2"; shift 2 ;;
    --case-val) CASE_VAL="$2"; shift 2 ;;
    --run-without-contrast) RUN_WITHOUT_CONTRAST="$2"; shift 2 ;;
    --n-folds) N_FOLDS="$2"; shift 2 ;;
    --n-cores) N_CORES="$2"; shift 2 ;;
    --chunk-workers) CHUNK_WORKERS="$2"; shift 2 ;;
    --io-workers) IO_WORKERS="$2"; shift 2 ;;
    --save-tmp) SAVE_TMP="$2"; shift 2 ;;
    --suffix) SUFFIX="$2"; shift 2 ;;
    --get-rank-purity) GET_RANK_PURITY="$2"; shift 2 ;;
    --purity-scale) PURITY_SCALE="$2"; shift 2 ;;
    --purity-layer) PURITY_LAYER="$2"; shift 2 ;;
    --derive-region-dummies) DERIVE_REGION="$2"; shift 2 ;;
    --compute-control-pcs) COMPUTE_PCS="$2"; shift 2 ;;
    --control-val) CONTROL_VAL="$2"; shift 2 ;;
    --n-hvgs) N_HVGS="$2"; shift 2 ;;
    --n-pcs) N_PCS="$2"; shift 2 ;;
    --pc-test-covs) PC_TEST_COVS="$2"; shift 2 ;;
    --pc-fdr) PC_FDR="$2"; shift 2 ;;
    --pca-model-dir) PCA_MODEL_DIR="$2"; shift 2 ;;
    --pc-max) PC_MAX="$2"; shift 2 ;;
    --pcs-to-add)
      shift
      while [[ $# -gt 0 && "$1" != --* ]]; do PCS_TO_ADD+=("${1//,/ }"); shift; done
      read -ra PCS_TO_ADD <<< "${PCS_TO_ADD[*]}"
      [[ ${#PCS_TO_ADD[@]} -eq 0 ]] && { echo "ERROR: --pcs-to-add needs one or more of: elbow significant all"; exit 1; } ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

for flag in GET_RANK_PURITY DERIVE_REGION COMPUTE_PCS RUN_WITHOUT_CONTRAST SAVE_TMP; do
  [[ "${!flag}" == "0" || "${!flag}" == "1" ]] || { echo "ERROR: ${flag} must be 0 or 1 (got '${!flag}')"; exit 1; }
done
[[ -z "${H5AD}" || -z "${OUTDIR}" || -z "${N_FOLDS}" || -z "${N_CORES}" ]] && { echo "ERROR: --path, --outdir, --n-folds, --n-cores required"; usage; exit 1; }
[[ ! -f "$H5AD" ]] && { echo "ERROR: File not found: $H5AD"; exit 1; }
[[ "$PURITY_SCALE" != "log-rate" && "$PURITY_SCALE" != "mean" ]] && { echo "ERROR: --purity-scale must be log-rate or mean"; exit 1; }
for k in "${PCS_TO_ADD[@]}"; do
  [[ "$k" == "elbow" || "$k" == "significant" || "$k" == "all" ]] || { echo "ERROR: --pcs-to-add accepts elbow / significant / all, got '$k'"; exit 1; }
done
if [[ ${#PCS_TO_ADD[@]} -gt 0 ]] && grep -qE '(^|,|:|\*)control_PC' <<< "$COVS"; then
  echo "ERROR: --covs already contains control_PC* terms; drop them or drop --pcs-to-add"; exit 1
fi
derive_flag=()
[[ "$DERIVE_REGION" == "1" ]] && derive_flag=(--derive-region-dummies)
command -v parallel > /dev/null || { echo "ERROR: GNU parallel not found"; exit 1; }

# -----------------------------
# Paths
# -----------------------------
H5AD="$(realpath "$H5AD")"
mkdir -p "$OUTDIR"
OUTDIR="$(realpath "$OUTDIR")"
base="$(basename "$H5AD")"
slogan="${base%.h5ad}"
tmp_dir="${OUTDIR}/${slogan}__nebula_ln"
raw_dir="${tmp_dir}/raw"
chunks_dir="${tmp_dir}/chunks"
tmp_de_dir="${tmp_dir}/de_results/"
if [[ "$PURITY_SCALE" == "log-rate" ]]; then
  pb_layer="$LAYER"; pb_dir="${tmp_dir}/pseudobulk_${LAYER}_lograte"
else
  pb_layer="$PURITY_LAYER"; pb_dir="${tmp_dir}/pseudobulk_${PURITY_LAYER}"
fi
final_dir="${OUTDIR}/de_results/nebula_ln"
log_dir="${OUTDIR}/logs"
mkdir -p "$final_dir" "$chunks_dir" "$tmp_de_dir" "$log_dir"

if [[ -n "${SUFFIX}" ]]; then
  final_qs="${final_dir}/${slogan}__nebula_ln_combined__${SUFFIX}.qs"
else
  final_qs="${final_dir}/${slogan}__nebula_ln_combined.qs"
fi

stamp="$(date +%Y%m%d_%H%M%S)"
log_file="${log_dir}/run_${slogan}${SUFFIX:+__${SUFFIX}}_${stamp}.log"
exec > >(tee -a "$log_file") 2>&1

echo "[INFO] Workflow: $WF_DIR (env: $ENV_PREFIX)"
echo "[INFO] Input h5ad: $H5AD (layer: $LAYER)"
echo "[INFO] Outdir: $OUTDIR"

# -----------------------------
# Step 0: control PCs (optional), then covariate list
# -----------------------------
if [[ "$COMPUTE_PCS" == "1" ]]; then
  echo "[INFO] Step 0: computing control PCs (contrast ${CONTRAST_COL}, control '${CONTROL_VAL}') -> h5ad updated in place"
  pc_args=(compute --h5ad "$H5AD" --layer "$LAYER" --id-col "$ID_COL" --contrast-col "$CONTRAST_COL"
           --control-val "$CONTROL_VAL" --case-val "$CASE_VAL" --n-hvgs "$N_HVGS" --n-pcs "$N_PCS"
           --pc-test-covs "$PC_TEST_COVS" --pc-fdr "$PC_FDR" --workers "$IO_WORKERS")
  [[ -n "$PCA_MODEL_DIR" ]] && pc_args+=(--pca-model-dir "$PCA_MODEL_DIR")
  "$PY" "${WF_DIR}/python/control_pcs.py" "${pc_args[@]}" || { echo "[ERROR] control_pcs.py failed"; exit 1; }
fi

COVS_FINAL="$COVS"
PCS_ADDED=""
PCA_CREATED=""
if [[ ${#PCS_TO_ADD[@]} -gt 0 ]]; then
  sel="$("$PY" "${WF_DIR}/python/control_pcs.py" select --h5ad "$H5AD" --pcs-to-add "${PCS_TO_ADD[@]}" --pc-max "$PC_MAX")" \
    || { echo "[ERROR] control PC selection failed"; exit 1; }
  PCS_ADDED="$(sed -n 1p <<< "$sel")"
  PCA_CREATED="$(sed -n 2p <<< "$sel")"
  echo "[INFO] PCs added (${PCS_TO_ADD[*]}, index <= ${PC_MAX}; control_pca created ${PCA_CREATED}): ${PCS_ADDED:-<none>}"
  [[ -n "$PCS_ADDED" ]] && COVS_FINAL="${COVS:+${COVS},}${PCS_ADDED}"
elif grep -qE '(^|,|:|\*)control_PC' <<< "$COVS"; then
  # PCs typed by hand: still tie the chunks to the PC values they were made from
  PCA_CREATED="$("$PY" -c 'import sys, h5py, anndata as ad
with h5py.File(sys.argv[1]) as h:
    cp = ad.io.read_elem(h["uns/control_pca"]) if "uns" in h and "control_pca" in h["uns"] else {}
print(cp.get("created", "unknown"))' "$H5AD")"
fi

# contrast must be a model term unless explicitly waived
IFS=',' read -ra COVS_ARR <<< "$COVS_FINAL"
has_contrast=0
for cov in "${COVS_ARR[@]}"; do [[ "$cov" == "$CONTRAST_COL" ]] && has_contrast=1; done
if [[ "$has_contrast" == "0" ]]; then
  if [[ "$RUN_WITHOUT_CONTRAST" == "1" ]]; then
    [[ "$GET_RANK_PURITY" == "1" ]] && { echo "ERROR: --get-rank-purity 1 needs the contrast coefficient; '${CONTRAST_COL}' is not in --covs"; exit 1; }
    echo "[WARN] '${CONTRAST_COL}' is not a covariate: no contrast coefficient will be fit (--run-without-contrast 1)"
  else
    echo "ERROR: --contrast-col '${CONTRAST_COL}' is not in --covs, so NEBULA would fit no coefficient for it."
    echo "       Add it to --covs, or pass --run-without-contrast 1 to fit a model without the contrast."
    exit 1
  fi
fi

# persist run args
"$PY" -c 'import json,sys; json.dump(dict(zip(sys.argv[2::2], sys.argv[3::2])), open(sys.argv[1], "w"), indent=2)' \
  "${log_dir}/run_args_${stamp}.json" path "$H5AD" layer "$LAYER" outdir "$OUTDIR" id_col "$ID_COL" covs "$COVS" \
  covs_final "$COVS_FINAL" pcs_to_add "${PCS_TO_ADD[*]}" pc_max "$PC_MAX" pcs_added "$PCS_ADDED" control_pca_created "$PCA_CREATED" \
  compute_control_pcs "$COMPUTE_PCS" control_val "$CONTROL_VAL" n_hvgs "$N_HVGS" n_pcs "$N_PCS" pc_test_covs "$PC_TEST_COVS" \
  pc_fdr "$PC_FDR" pca_model_dir "$PCA_MODEL_DIR" offset_col "$OFFSET_COL" \
  contrast_col "$CONTRAST_COL" case_val "$CASE_VAL" run_without_contrast "$RUN_WITHOUT_CONTRAST" n_folds "$N_FOLDS" n_cores "$N_CORES" \
  chunk_workers "$CHUNK_WORKERS" io_workers "$IO_WORKERS" save_tmp "$SAVE_TMP" suffix "$SUFFIX" \
  get_rank_purity "$GET_RANK_PURITY" purity_scale "$PURITY_SCALE" purity_layer "$PURITY_LAYER" derive_region_dummies "$DERIVE_REGION" \
  command "$0 ${ARGS_ORIG[*]}"

echo "[INFO] Rscript: $RS"
echo "[INFO] n-folds: $N_FOLDS   nebula jobs: $N_CORES   chunk writers: $CHUNK_WORKERS   io workers: $IO_WORKERS"
echo "[INFO] id-col: $ID_COL"
echo "[INFO] covs (final):"
for cov in "${COVS_ARR[@]}"; do echo "       - $cov"; done
echo "[INFO] offset-col: $OFFSET_COL  -> meta 'nUMI_nebula'"
echo "[INFO] rank purity: $GET_RANK_PURITY (scale: $PURITY_SCALE, layer: $pb_layer)"
echo "[INFO] derive region dummies from region: $DERIVE_REGION"
echo "[INFO] Tmp dir: $tmp_dir"
echo "[INFO] Final combined results: $final_qs"
echo "[INFO] Log: $log_file"
echo

# -----------------------------
# Step 1-3: Chunk straight from the h5ad
# -----------------------------
chunk_params="layer=${LAYER};n_folds=${N_FOLDS};id_col=${ID_COL};covs=${COVS_FINAL};offset_col=${OFFSET_COL};derive_region=${DERIVE_REGION};control_pca_created=${PCA_CREATED}"
params_file="${chunks_dir}/chunk_params.txt"
n_existing=$(find "$chunks_dir" -maxdepth 1 -name "chunk_*_of_*.qs" | wc -l)

if [[ -f "$params_file" && "$(cat "$params_file")" != "$chunk_params" ]]; then
  echo "[ERROR] Existing chunks in ${chunks_dir} were made with different parameters:"
  echo "        existing: $(cat "$params_file")"
  echo "        current:  ${chunk_params}"
  echo "        Remove ${tmp_dir} or use another --outdir."
  exit 1
fi

if [[ "$n_existing" -eq "$N_FOLDS" && -f "$params_file" ]]; then
  echo "[INFO] All ${N_FOLDS} chunk files already exist in ${chunks_dir}, skipping chunking."
else
  echo "$chunk_params" > "$params_file"
  if [[ ! -f "${raw_dir}/summary.json" ]]; then
    rm -rf "$raw_dir"
    echo "[INFO] Streaming ${LAYER} into ${N_FOLDS} gene chunks..."
    "$PY" "${WF_DIR}/python/h5ad_to_gene_chunks.py" \
      --h5ad "$H5AD" --layer "$LAYER" --raw-dir "$raw_dir" --n-folds "$N_FOLDS" \
      --id-col "$ID_COL" --covs "$COVS_FINAL" --offset-col "$OFFSET_COL" --workers "$IO_WORKERS" "${derive_flag[@]}" \
      || { echo "[ERROR] h5ad_to_gene_chunks.py failed"; exit 1; }
  fi
  echo "[INFO] Writing chunk .qs files (${CHUNK_WORKERS} at a time)..."
  seq 1 "$N_FOLDS" | parallel --halt now,fail=1 -j "$CHUNK_WORKERS" \
    "$RS" "${WF_DIR}/R/write_chunk_qs.R" --raw-dir "$raw_dir" --chunk {} --outdir "$chunks_dir" \
    || { echo "[ERROR] write_chunk_qs.R failed"; exit 1; }
  n_existing=$(find "$chunks_dir" -maxdepth 1 -name "chunk_*_of_*.qs" | wc -l)
  [[ "$n_existing" -ne "$N_FOLDS" ]] && { echo "[ERROR] expected ${N_FOLDS} chunk files, found ${n_existing}"; exit 1; }
  cp "${raw_dir}/chunk_plan.tsv" "${raw_dir}/summary.json" "$chunks_dir/"
  rm -rf "$raw_dir"
fi
echo

# -----------------------------
# Rank purity pseudobulk (runs alongside NEBULA)
# -----------------------------
pb_pid=""
if [[ "$GET_RANK_PURITY" == "1" && ! -f "${pb_dir}/pb_mean.npy" ]]; then
  echo "[INFO] Pseudobulking ${pb_layer} (${PURITY_SCALE}) in the background -> ${pb_dir}"
  "$PY" "${WF_DIR}/python/pseudobulk_h5ad.py" --h5ad "$H5AD" --layer "$pb_layer" --scale "$PURITY_SCALE" \
    --offset-col "$OFFSET_COL" --id-col "$ID_COL" \
    --covs "$COVS_FINAL" --contrast-col "$CONTRAST_COL" --out-dir "$pb_dir" --workers "$IO_WORKERS" "${derive_flag[@]}" &
  pb_pid=$!
fi

# -----------------------------
# Step 4: run_nebula.R per chunk (GNU parallel)
# -----------------------------
echo "[INFO] Launching NEBULA chunk jobs with GNU parallel..."
jobs_file="${tmp_dir}/nebula_jobs.txt"
: > "$jobs_file"
for chunk in "$chunks_dir"/chunk_*_of_*.qs; do
  cmd=("$RS" "${WF_DIR}/R/run_nebula.R" --path "$chunk" --id-col "$ID_COL" --offset-col nUMI_nebula
       --n-cores 1 --out-dir "$tmp_de_dir")
  [[ -n "$COVS_FINAL" ]] && cmd+=(--covs "$COVS_FINAL")
  [[ -n "$SUFFIX" ]] && cmd+=(--suffix "$SUFFIX")
  printf '%q ' "${cmd[@]}" >> "$jobs_file"
  echo >> "$jobs_file"
done
parallel -j "$N_CORES" < "$jobs_file" || { echo "[ERROR] one or more run_nebula.R jobs failed"; exit 1; }
echo "[INFO] All chunk jobs completed."
echo

# -----------------------------
# Step 5: combine, then rank purity from the h5ad
# -----------------------------
combine_args=(--in-dir "$tmp_de_dir" --out "$final_qs")
[[ -n "$SUFFIX" ]] && combine_args+=(--suffix "$SUFFIX")
"$RS" "${WF_DIR}/R/combine_nebula_results.R" "${combine_args[@]}" \
  || { echo "[ERROR] combine_nebula_results.R failed"; exit 1; }

if [[ "$GET_RANK_PURITY" == "1" ]]; then
  if [[ -n "$pb_pid" ]]; then
    wait "$pb_pid" || { echo "[ERROR] pseudobulk_h5ad.py failed"; exit 1; }
  fi
  purity_tsv="${final_qs%.qs}__rank_purity.tsv"
  "$PY" "${WF_DIR}/python/rank_purity.py" --combined-csv "${final_qs%.qs}.csv" --pb-dir "$pb_dir" \
    --covs "$COVS_FINAL" --contrast-col "$CONTRAST_COL" --case-val "$CASE_VAL" --out "$purity_tsv" \
    || { echo "[ERROR] rank_purity.py failed"; exit 1; }
  "$RS" "${WF_DIR}/R/add_rank_purity_to_combined.R" --combined-qs "$final_qs" --purity-tsv "$purity_tsv" \
    || { echo "[ERROR] add_rank_purity_to_combined.R failed"; exit 1; }
fi
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
