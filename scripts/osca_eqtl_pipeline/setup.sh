#!/usr/bin/env bash
# setup.sh — One-time setup: installs micromamba and creates the osca-venv conda
# environment with all Python and R dependencies for the OSCA eQTL pipeline.
# Everything is confined to cc-sandbox.
#
# Usage (from any directory):
#   bash /path/to/cc-sandbox/osca_eqtl_pipeline_patches/setup.sh

set -euo pipefail

SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"   # cc-sandbox root
BIN_DIR="${SANDBOX}/bin"
MAMBA_ROOT="${SANDBOX}/micromamba_root"
MAMBA="${BIN_DIR}/micromamba"
ENV_NAME="osca-venv"
ENV_PREFIX="${MAMBA_ROOT}/envs/${ENV_NAME}"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT}"

echo "════════════════════════════════════════════════════════"
echo " OSCA eQTL pipeline — environment setup"
echo " Sandbox : ${SANDBOX}"
echo " Env     : ${ENV_NAME}"
echo "════════════════════════════════════════════════════════"

# ── 1. micromamba binary ──────────────────────────────────────────────────────
mkdir -p "${BIN_DIR}"
if [ ! -f "${MAMBA}" ]; then
  echo "[micromamba] Downloading binary …"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvjf - -C "${BIN_DIR}" --strip-components=1 bin/micromamba
else
  echo "[micromamba] Already present: $("${MAMBA}" --version)"
fi

# ── 2. Create env (skip if already exists at any known location) ──────────────
_EXISTING="/home/ferris/sc-online/scripts/micromamba_root/envs/${ENV_NAME}"
if "${MAMBA}" env list | grep -q "^${ENV_NAME}" || [ -d "${_EXISTING}/bin" ]; then
  echo "[conda] '${ENV_NAME}' already exists — skipping creation."
  # Point ENV_PREFIX at whichever location exists
  [ -d "${ENV_PREFIX}/bin" ] || ENV_PREFIX="${_EXISTING}"
else
  echo "[conda] Creating environment '${ENV_NAME}' …"
  # CONDA_PKGS_DIRS stored within MAMBA_ROOT (persistent, no /tmp dependency).
  # r-mashr and r-ashr are on conda-forge; bioconductor-sva is on bioconda.
  # If the solve fails because r-mashr is unavailable for the resolved R version,
  # add r-base=4.3 to the package list below.
  CONDA_PKGS_DIRS="${MAMBA_ROOT}/conda-pkgs" \
  "${MAMBA}" create -n "${ENV_NAME}" -y \
    -c conda-forge -c bioconda \
    "python=3.12" \
    numpy pandas anndata scanpy "sccoda=0.1.9" decoupler-py \
    scipy matplotlib seaborn gseapy pyreadr ncls \
    r-base r-getopt r-mashr r-ashr r-matrix \
    r-tidyverse r-ggrepel r-repr r-cairo \
    bioconductor-sva bioconductor-qvalue \
    parallel plink2
fi

# ── 2b. Packages to ensure are present in existing envs ──────────────────────
# Runs unconditionally so re-running setup patches an already-created env.
CONDA_PKGS_DIRS="${MAMBA_ROOT}/conda-pkgs" \
"${MAMBA}" install -n "${ENV_NAME}" -y -c conda-forge -c bioconda \
    r-repr r-cairo plink2 \
    scipy matplotlib seaborn gseapy pyreadr ncls \
    bioconductor-qvalue

# ── 3. Symlink OSCA binary into env bin (so subshells find it on PATH) ────────
OSCA_SRC="/mnt/accessory/analysis/eqtl/osca/osca"
OSCA_DST="${ENV_PREFIX}/bin/osca"
if [ ! -e "${OSCA_DST}" ]; then
  if [ ! -f "${OSCA_SRC}" ]; then
    echo "WARNING: OSCA binary not found at ${OSCA_SRC}. Add to PATH manually."
  else
    ln -s "${OSCA_SRC}" "${OSCA_DST}"
    echo "[osca] Symlinked ${OSCA_SRC} → ${OSCA_DST}"
  fi
fi

# ── 4. R package verification ─────────────────────────────────────────────────
echo "[R] Verifying installed R packages …"
"${MAMBA}" run -n "${ENV_NAME}" \
  Rscript --no-save --no-restore --no-site-file --no-init-file -e "
  cat('\n=== R package verification ===\n')
  pkgs <- c('mashr', 'ashr', 'sva', 'Matrix', 'getopt',
            'tidyverse', 'dplyr', 'ggplot2', 'ggrepel', 'tidyr', 'repr', 'Cairo',
            'qvalue')
  missing <- c()
  for (pkg in pkgs) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf('  %-30s %s\n', pkg, if (ok) 'OK' else 'MISSING'))
    if (!ok) missing <- c(missing, pkg)
  }
  if (length(missing) > 0) stop(paste('Missing R packages:', paste(missing, collapse = ', ')))
  cat('All R packages OK.\n')
"

# ── 5. Python package verification ───────────────────────────────────────────
echo "[Python] Verifying installed Python packages …"
"${MAMBA}" run -n "${ENV_NAME}" \
  python -c "
import sys
pkgs = ['numpy', 'pandas', 'anndata', 'scanpy', 'sccoda', 'decoupler',
        'scipy', 'matplotlib', 'seaborn', 'gseapy', 'pyreadr', 'ncls']
missing = []
for p in pkgs:
    try:
        __import__(p)
        print(f'  {p:<30} OK')
    except ImportError:
        print(f'  {p:<30} MISSING')
        missing.append(p)
if missing:
    sys.exit(f'Missing Python packages: {missing}')
print('All Python packages OK.')
"

# ── 6. Download external reference resources ─────────────────────────────────
RESOURCES="${SANDBOX}/resources"
mkdir -p "${RESOURCES}"

# GTEx v8 Brain_Substantia_nigra significant eQTL pairs (~4 MB gzipped)
# File lives inside the full eQTL tar at the adult-gtex GCS bucket.
# We stream the tar and extract only the SN file to avoid downloading 1.5 GB.
GTEX_SN="${RESOURCES}/gtex_sn_signif_pairs.txt.gz"
GTEX_TAR_URL="https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar"
GTEX_TAR_PATH="GTEx_Analysis_v8_eQTL/Brain_Substantia_nigra.v8.signif_variant_gene_pairs.txt.gz"
if [ ! -f "${GTEX_SN}" ]; then
    echo "[resources] Streaming GTEx v8 eQTL tar, extracting SN significant pairs …"
    GTEX_TMP="$(mktemp -d "${RESOURCES}/.gtex_tmp_XXXXXX")"
    curl -fL "${GTEX_TAR_URL}" \
        | tar -x "${GTEX_TAR_PATH}" --strip-components=1 -C "${GTEX_TMP}" \
        && mv "${GTEX_TMP}/Brain_Substantia_nigra.v8.signif_variant_gene_pairs.txt.gz" "${GTEX_SN}" \
        && rm -rf "${GTEX_TMP}" \
        && echo "[resources] GTEx SN download complete: ${GTEX_SN}" \
        || { echo "WARNING: GTEx SN download failed. Provide manually: ${GTEX_SN}"; rm -rf "${GTEX_TMP}"; }
else
    echo "[resources] GTEx SN file already present: ${GTEX_SN}"
fi

# Corces 2020 snATAC-seq peaks from dopaminergic neuron dataset (GEO GSE147672)
# IDR peaks are per-cluster; we merge all 24 clusters into a single BED with
# col4 = ClusterN for use in cell-type-level enrichment analysis.
CORCES_ATAC="${RESOURCES}/corces_2020_da_atac_peaks.bed.gz"
CORCES_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147672/suppl/GSE147672_scATAC_idr_peaks.tar.gz"
if [ ! -f "${CORCES_ATAC}" ]; then
    echo "[resources] Downloading and merging Corces 2020 snATAC IDR peaks …"
    CORCES_TMP="$(mktemp -d "${RESOURCES}/.corces_tmp_XXXXXX")"
    curl -fL "${CORCES_URL}" | tar -xz -C "${CORCES_TMP}" \
        && "${MAMBA}" run -n "${ENV_NAME}" python -c "
import gzip, re, glob
files = sorted(glob.glob('${CORCES_TMP}/*.gz'),
               key=lambda x: int(re.search(r'Cluster(\d+)', x).group(1)))
with gzip.open('${CORCES_ATAC}', 'wt') as fout:
    for f in files:
        cluster = re.search(r'Cluster(\d+)', f).group(1)
        with gzip.open(f, 'rt') as fin:
            for line in fin:
                p = line.strip().split('\t')
                fout.write(f'{p[0]}\t{p[1]}\t{p[2]}\tCluster{cluster}\n')
print('Done')
" \
        && rm -rf "${CORCES_TMP}" \
        && echo "[resources] Corces ATAC download complete: ${CORCES_ATAC}" \
        || { echo "WARNING: Corces ATAC download failed. Provide manually: ${CORCES_ATAC}"; rm -rf "${CORCES_TMP}"; }
else
    echo "[resources] Corces ATAC file already present: ${CORCES_ATAC}"
fi

echo ""
echo "════════════════════════════════════════════════════════"
echo " Setup complete."
echo " Run the pipeline (no env activation needed):"
echo "   bash $(dirname "$0")/00-run-osca-eqtl-pipeline.sh \\"
echo "     --input-files /path/to/h5ad-paths.txt \\"
echo "     --osca-input-dir /path/to/osca-inputs \\"
echo "     --pb-output-dir /path/to/pseudobulk \\"
echo "     --pipeline-slogan my-run \\"
echo "     --vcf-slogan my-vcf \\"
echo "     --participants /path/to/participants.txt"
echo "════════════════════════════════════════════════════════"
