#!/usr/bin/env bash
# setup_env.sh
# One-time setup: installs micromamba, creates a conda environment with
# R 4.4 + all crumblr dependencies, and sets up a Python venv.
# Everything is confined to cc-sandbox.
#
# Usage (from any directory):
#   bash /path/to/crumblr/setup_env.sh

set -euo pipefail

SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT_DIR="${SANDBOX}/crumblr"
BIN_DIR="${SANDBOX}/bin"
MAMBA_ROOT="${SANDBOX}/micromamba_root"
VENV_DIR="${SANDBOX}/crumblr-venv"
ENV_NAME="crumblr"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT}"
MAMBA="${BIN_DIR}/micromamba"

echo "════════════════════════════════════════════════════════"
echo " crumblr workflow — environment setup"
echo " Sandbox : ${SANDBOX}"
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

# ── 2. conda environment with R 4.4 + Bioconductor packages ──────────────────
# If environment.yml (pinned lock file) exists, use it for a reproducible build.
# Otherwise fall back to installing by package name (will resolve latest versions).
if "${MAMBA}" env list | grep -q "^${ENV_NAME}"; then
  echo "[conda] Environment '${ENV_NAME}' already exists — skipping creation."
elif [ -f "${SCRIPT_DIR}/environment.yml" ]; then
  echo "[conda] Creating environment '${ENV_NAME}' from pinned environment.yml …"
  "${MAMBA}" env create -n "${ENV_NAME}" -y \
    -f "${SCRIPT_DIR}/environment.yml"
else
  echo "[conda] Creating environment '${ENV_NAME}' with R 4.4 (latest packages) …"
  "${MAMBA}" create -n "${ENV_NAME}" -y \
    -c conda-forge -c bioconda \
    r-base=4.4 \
    r-lme4 \
    r-ggplot2 \
    r-remotes \
    bioconductor-variancepartition \
    bioconductor-edger \
    bioconductor-limma \
    bioconductor-biocparallel \
    bioconductor-summarizedexperiment \
    bioconductor-singlecellexperiment \
    bioconductor-genomicranges \
    bioconductor-delayedarray \
    bioconductor-genomeinfoDb \
    bioconductor-ggtree
fi

# ── 3. Install BiocManager, upgrade variancePartition, install crumblr ───────
# conda-forge only ships variancePartition 1.36.2; crumblr requires >=1.36.3.
# Install BiocManager from CRAN and use it to pull the newer Bioconductor build.
echo "[R] Installing BiocManager, upgrading variancePartition, installing crumblr …"
"${MAMBA}" run -n "${ENV_NAME}" Rscript --no-save --no-restore -e "
  lib <- .libPaths()[1]
  install.packages('BiocManager', lib=lib, repos='https://cloud.r-project.org',
                   quiet=TRUE)
  BiocManager::install('variancePartition', ask=FALSE, update=TRUE, force=TRUE)
  remotes::install_github('GabrielHoffman/crumblr', upgrade='never', quiet=FALSE)
  cat('\n=== Verification ===\n')
  for (pkg in c('crumblr','variancePartition','lme4','edgeR','limma',
                'SingleCellExperiment','ggtree','compositions')) {
    ok <- requireNamespace(pkg, quietly=TRUE)
    cat(sprintf('  %-30s %s\n', pkg, if(ok) 'OK' else 'MISSING'))
  }
"

# ── 4. Python venv (system Python + pip bootstrap) ───────────────────────────
if [ -f "${VENV_DIR}/bin/activate" ]; then
  echo "[Python] venv already exists at ${VENV_DIR} — skipping."
else
  echo "[Python] Creating venv …"
  # Locate a python3 that has the venv module
  PYTHON3=""
  for _py in python3 python; do
    if command -v "${_py}" >/dev/null 2>&1 && "${_py}" -c "import venv" 2>/dev/null; then
      PYTHON3=$(command -v "${_py}")
      break
    fi
  done
  if [ -z "${PYTHON3}" ]; then
    echo "ERROR: Could not find python3 with the 'venv' module on PATH."
    exit 1
  fi
  echo "[Python] Using: ${PYTHON3}"
  "${PYTHON3}" -m venv --without-pip "${VENV_DIR}"
  curl -sS https://bootstrap.pypa.io/get-pip.py | "${VENV_DIR}/bin/python3"
fi

source "${VENV_DIR}/bin/activate"
echo "[Python] Installing packages …"
pip install --quiet -r "${SCRIPT_DIR}/requirements.txt"
deactivate

echo ""
echo "════════════════════════════════════════════════════════"
echo " Setup complete."
echo " Run analysis (paths auto-detected — no env vars needed):"
echo "   source ${VENV_DIR}/bin/activate"
echo "   python ${SCRIPT_DIR}/run_workflow.py \\"
echo "     --csv data.csv --proportion-col cell_type --sample-col sample_id \\"
echo "     --contrast case_control --fixed-effects age sex ... \\"
echo "     --random-effects donor_id --output-dir results/"
echo "════════════════════════════════════════════════════════"
