#!/usr/bin/env bash
set -euo pipefail

BCL=$1
INDEX=$2

SCRIPT_DIR=$(dirname "$(realpath "$0")")
GS_BASE=gs://fc-secure-d99fbd65-eb27-4989-95b4-4cf559aa7d36/gene-expression
BASE_PATH=/mnt/accessory/seq_data/pd-freeze/gene-expression 
CB_DIR="${BASE_PATH}/${BCL}/${INDEX}/cellbender"
DS_DIR="${BASE_PATH}/${BCL}/${INDEX}/dropsift_outputs"
CR_DIR="${BASE_PATH}/${BCL}/${INDEX}/outs"
VIREO_DIR="${BASE_PATH}/${BCL}/${INDEX}/vireo"


LOG_FILE="${BASE_PATH}/run-load-data-errors-2"

# add date to end of log file name
LOG_FILE="${LOG_FILE}-$(date '+%Y%m%d').log"

if [ ! -e "${LOG_FILE}" ]; then
    touch "${LOG_FILE}"
fi

mkdir -p "$CB_DIR" "$DS_DIR" "$CR_DIR" "$VIREO_DIR"

echo "Loading data for BCL: ${BCL}, INDEX: ${INDEX}"
echo "Downloading CellBender output..."
gcloud storage cp "${GS_BASE}/${BCL}/${INDEX}/cellbender/*_out_filtered.h5" "${CB_DIR}/" || {
# gcloud storage ls "${GS_BASE}/${BCL}/${INDEX}/cellbender/*_out_filtered.h5"|| {
    echo "{$(date '+%Y-%m-%d %H:%M:%S')}: Failed to download CellBender Output for BCL: ${BCL}, INDEX: ${INDEX}" >> "${LOG_FILE}"
    exit 1
}

echo "Downloading dropsift output..."
gcloud storage cp "${GS_BASE}/${BCL}/${INDEX}/dropsift_outputs/dropsift_output.csv" "${DS_DIR}/" || {
# gcloud storage ls "${GS_BASE}/${BCL}/${INDEX}/dropsift_outputs/dropsift_output.csv" || {
    echo "{$(date '+%Y-%m-%d %H:%M:%S')}: Failed to download dropsift output for BCL: ${BCL}, INDEX: ${INDEX}" >> "${LOG_FILE}"
    exit 1
}

echo "Downloading dropsift CellRanger Molecule Info..."
gcloud storage cp "${GS_BASE}/${BCL}/${INDEX}/outs/molecule_info.h5" "${CR_DIR}/" || {
# gcloud storage ls "${GS_BASE}/${BCL}/${INDEX}/outs/molecule_info.h5" || {
    echo "{$(date '+%Y-%m-%d %H:%M:%S')}: Failed to download CellRanger Molecule Info for BCL: ${BCL}, INDEX: ${INDEX}" >> "${LOG_FILE}"
    exit 1
}

echo "Downloading Vireo donor IDs..."
gcloud storage cp "${GS_BASE}/${BCL}/${INDEX}/vireo/05_VIREO/vireo_final_out/donor_ids.tsv" "${VIREO_DIR}/" || {
#gcloud storage ls "${GS_BASE}/${BCL}/${INDEX}/vireo/05_VIREO/vireo_final_out/donor_ids.tsv" || {
    echo "{$(date '+%Y-%m-%d %H:%M:%S')}: Failed to download Vireo donor IDs for BCL: ${BCL}, INDEX: ${INDEX}" >> "${LOG_FILE}"
    exit 1
}

echo "building anndata with load_data.py..."
python "${SCRIPT_DIR}"/load_data.py --basepath "${BASE_PATH}" --bcl-index "${BCL}/${INDEX}"  && {
    echo "Successfully loaded data for BCL: ${BCL}, INDEX: ${INDEX}"
    echo "Successfully loaded data for BCL: ${BCL}, INDEX: ${INDEX}" >> "${LOG_FILE}"
    rm -r "${CB_DIR}" "${DS_DIR}" "${CR_DIR}" "${VIREO_DIR}"
} || {
    echo "{$(date '+%Y-%m-%d %H:%M:%S')}: Failed to load data for BCL: ${BCL}, INDEX: ${INDEX}" >> "${LOG_FILE}"
    exit 1
}