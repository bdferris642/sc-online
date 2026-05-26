#!/usr/bin/env bash
# setup.sh
# One-time setup: installs micromamba and creates the cnmf-analysis-env conda
# environment with R 4.4 + all pipeline dependencies.
# Everything is confined to cc-sandbox.
#
# Usage (from any directory):
#   bash /path/to/cc-sandbox/cnmf-analysis/setup.sh

set -euo pipefail

SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT_DIR="${SANDBOX}/cnmf-analysis"
BIN_DIR="${SANDBOX}/bin"
MAMBA_ROOT="${SANDBOX}/micromamba_root"
ENV_NAME="cnmf-analysis-env"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT}"
MAMBA="${BIN_DIR}/micromamba"

echo "════════════════════════════════════════════════════════"
echo " cnmf-analysis pipeline — environment setup"
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

# ── 2. conda environment with R 4.4 + packages ───────────────────────────────
if "${MAMBA}" env list | grep -q "^${ENV_NAME}"; then
  echo "[conda] Environment '${ENV_NAME}' already exists — skipping creation."
else
  echo "[conda] Creating environment '${ENV_NAME}' with R 4.4 …"
  # CONDA_PKGS_DIRS stored within MAMBA_ROOT (persistent, no /tmp dependency)
  CONDA_PKGS_DIRS="${MAMBA_ROOT}/conda-pkgs" \
  "${MAMBA}" create -n "${ENV_NAME}" -y \
    -c conda-forge -c bioconda \
    "r-base=4.4" \
    r-optparse \
    r-ggplot2 \
    r-ggrepel \
    r-glue \
    r-getopt \
    r-dplyr \
    r-stringr \
    r-scales \
    r-igraph \
    r-betareg \
    r-nnls \
    r-svglite \
    r-tidyr \
    r-seuratobject \
    r-qs \
    r-qs2 \
    bioconductor-fgsea \
    python \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    numpy
fi

# ── 3. Package verification ───────────────────────────────────────────────────
echo "[R] Verifying installed R packages …"
"${MAMBA}" run -n "${ENV_NAME}" \
  Rscript --no-save --no-restore --no-site-file --no-init-file -e "
  cat('\n=== R package verification ===\n')
  pkgs <- c(
    'optparse', 'ggplot2', 'ggrepel', 'glue', 'getopt',
    'dplyr', 'stringr', 'scales', 'igraph',
    'betareg', 'nnls', 'svglite', 'tidyr',
    'SeuratObject', 'qs', 'qs2',
    'fgsea'
  )
  missing <- c()
  for (pkg in pkgs) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf('  %-35s %s\n', pkg, if (ok) 'OK' else 'MISSING'))
    if (!ok) missing <- c(missing, pkg)
  }
  if (length(missing) > 0) {
    stop(paste('Missing packages:', paste(missing, collapse = ', ')))
  }
  cat('\nAll R packages OK.\n')
"

echo "[Python] Verifying installed Python packages …"
"${MAMBA}" run -n "${ENV_NAME}" \
  python -c "
import sys
pkgs = ['pandas', 'matplotlib', 'seaborn', 'scipy', 'numpy']
missing = []
for pkg in pkgs:
    try:
        __import__(pkg)
        print(f'  {pkg:<35} OK')
    except ImportError:
        print(f'  {pkg:<35} MISSING')
        missing.append(pkg)
if missing:
    sys.exit(f'Missing Python packages: {missing}')
print()
print('All Python packages OK.')
"

echo ""
echo "════════════════════════════════════════════════════════"
echo " Setup complete."
echo " Run the full pipeline (no env activation needed):"
echo "   python ${SCRIPT_DIR}/run_gep_analysis.py \\"
echo "     --manifest /path/to/manifest.tsv \\"
echo "     --outdir   /path/to/output"
echo "════════════════════════════════════════════════════════"
