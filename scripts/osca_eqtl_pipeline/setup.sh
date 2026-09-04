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

# ── 2. Create env (skip if already exists) ────────────────────────────────────
if "${MAMBA}" env list | grep -q "^${ENV_NAME}"; then
  echo "[conda] '${ENV_NAME}' already exists — skipping creation."
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
    r-base r-getopt r-mashr r-ashr r-matrix \
    r-tidyverse r-ggrepel r-repr r-cairo \
    bioconductor-sva \
    parallel
fi

# ── 2b. Packages to ensure are present in existing envs ──────────────────────
# Runs unconditionally so re-running setup patches an already-created env.
CONDA_PKGS_DIRS="${MAMBA_ROOT}/conda-pkgs" \
"${MAMBA}" install -n "${ENV_NAME}" -y -c conda-forge r-repr r-cairo

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
            'tidyverse', 'dplyr', 'ggplot2', 'ggrepel', 'tidyr', 'repr', 'Cairo')
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
pkgs = ['numpy', 'pandas', 'anndata', 'scanpy', 'sccoda', 'decoupler']
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
