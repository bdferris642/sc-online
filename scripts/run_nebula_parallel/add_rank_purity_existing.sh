#!/usr/bin/env bash
# Add rank_purity to an EXISTING combined NEBULA result (e.g. from upstream run-nebula-parallel.sh),
# computed from the h5ad the model was fit on. The source .qs/.csv are never modified: copies are
# written to <outdir>/de_results/nebula_ln/ with the rank_purity column joined, plus the purity TSV.
#
# Usage:
#   add_rank_purity_existing.sh --h5ad X.h5ad --combined-qs Y__nebula_ln_combined__S.qs \
#     --covs case_control,age,... --outdir DIR [--id-col participant_id] [--contrast-col case_control] [--case-val pd] \
#     [--layer counts] [--offset-col nCount_RNA] [--purity-scale log-rate|mean] [--purity-layer norm_log] \
#     [--min-cells-per-donor 5] [--io-workers 24]
# --covs must be the exact --covs the NEBULA model used: rank_purity.py errors if any model
# coefficient is not claimed by a covariate.
set -euo pipefail

WF_DIR="$(cd "$(dirname "$(realpath "$0")")" && pwd)"
ENV_PREFIX="$(cd "${WF_DIR}/.." && pwd)/micromamba_root/envs/run-nebula-parallel-env"
[[ -x "${ENV_PREFIX}/bin/python" ]] || { echo "ERROR: workflow env not found at ${ENV_PREFIX}; run ${WF_DIR}/setup.sh"; exit 1; }
export PATH="${ENV_PREFIX}/bin:${PATH}"
PY="${ENV_PREFIX}/bin/python"

H5AD=""; COMBINED_QS=""; COVS=""; OUTDIR=""
ID_COL="participant_id"; CONTRAST_COL="case_control"; CASE_VAL="pd"; LAYER="counts"; OFFSET_COL="nCount_RNA"
PURITY_SCALE="log-rate"; PURITY_LAYER="norm_log"; MIN_CELLS=5; IO_WORKERS=24

while [[ $# -gt 0 ]]; do
  case "$1" in
    --h5ad) H5AD="$2"; shift 2 ;;
    --combined-qs) COMBINED_QS="$2"; shift 2 ;;
    --covs) COVS="$2"; shift 2 ;;
    --id-col) ID_COL="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --contrast-col) CONTRAST_COL="$2"; shift 2 ;;
    --case-val) CASE_VAL="$2"; shift 2 ;;
    --layer) LAYER="$2"; shift 2 ;;
    --offset-col) OFFSET_COL="$2"; shift 2 ;;
    --purity-scale) PURITY_SCALE="$2"; shift 2 ;;
    --purity-layer) PURITY_LAYER="$2"; shift 2 ;;
    --min-cells-per-donor) MIN_CELLS="$2"; shift 2 ;;
    --io-workers) IO_WORKERS="$2"; shift 2 ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done
[[ -z "$H5AD" || -z "$COMBINED_QS" || -z "$COVS" || -z "$OUTDIR" ]] && { echo "ERROR: --h5ad, --combined-qs, --covs, --outdir required"; exit 1; }
[[ -f "$H5AD" ]] || { echo "ERROR: not found: $H5AD"; exit 1; }
[[ -f "$COMBINED_QS" && -f "${COMBINED_QS%.qs}.csv" ]] || { echo "ERROR: need $COMBINED_QS and its .csv"; exit 1; }

final_dir="${OUTDIR}/de_results/nebula_ln"; log_dir="${OUTDIR}/logs"
mkdir -p "$final_dir" "$log_dir"
stem="$(basename "${COMBINED_QS%.qs}")"
out_qs="${final_dir}/${stem}.qs"
purity_tsv="${final_dir}/${stem}__rank_purity.tsv"
pb_dir="${OUTDIR}/${stem}__pseudobulk"
[[ -e "$out_qs" ]] && { echo "ERROR: $out_qs exists; remove it to recompute"; exit 1; }

exec > >(tee -a "${log_dir}/add_rank_purity_${stem}_$(date +%Y%m%d_%H%M%S).log") 2>&1
echo "[INFO] h5ad: $H5AD"
echo "[INFO] source: $COMBINED_QS"
echo "[INFO] covs: $COVS"
echo "[INFO] scale: $PURITY_SCALE  offset: $OFFSET_COL  min cells/donor: $MIN_CELLS"

if [[ "$PURITY_SCALE" == "log-rate" ]]; then pb_layer="$LAYER"; else pb_layer="$PURITY_LAYER"; fi
"$PY" "${WF_DIR}/python/pseudobulk_h5ad.py" --h5ad "$H5AD" --layer "$pb_layer" --scale "$PURITY_SCALE" \
  --offset-col "$OFFSET_COL" --id-col "$ID_COL" --covs "$COVS" --contrast-col "$CONTRAST_COL" \
  --out-dir "$pb_dir" --workers "$IO_WORKERS"
"$PY" "${WF_DIR}/python/rank_purity.py" --combined-csv "${COMBINED_QS%.qs}.csv" --pb-dir "$pb_dir" --covs "$COVS" \
  --contrast-col "$CONTRAST_COL" --case-val "$CASE_VAL" --min-cells-per-donor "$MIN_CELLS" --out "$purity_tsv"
cp "$COMBINED_QS" "$out_qs"
cp "${COMBINED_QS%.qs}.csv" "${out_qs%.qs}.csv"
"${ENV_PREFIX}/bin/Rscript" "${WF_DIR}/R/add_rank_purity_to_combined.R" --combined-qs "$out_qs" --purity-tsv "$purity_tsv"
rm -rf "$pb_dir"
echo "[INFO] done: $out_qs"
