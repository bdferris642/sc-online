#!/usr/bin/env bash
# setup.sh
# One-time setup: installs micromamba and creates the de-depletion-env conda
# environment with R 4.4 + all pipeline dependencies.
# Everything is confined to cc-sandbox.
#
# Usage (from any directory):
#   bash /path/to/cc-sandbox/de-depletion-analysis/setup.sh

set -euo pipefail

SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT_DIR="${SANDBOX}/de-depletion-analysis"
BIN_DIR="${SANDBOX}/bin"
MAMBA_ROOT="${SANDBOX}/micromamba_root"
ENV_NAME="de-depletion-env"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT}"
MAMBA="${BIN_DIR}/micromamba"
CONDA_PKGS_DIRS="${MAMBA_ROOT}/conda-pkgs"    # persistent cache, NOT /tmp

echo "════════════════════════════════════════════════════════"
echo " de-depletion-analysis pipeline — environment setup"
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
    r-dplyr \
    r-tidyr \
    r-tibble \
    r-ggplot2 \
    r-ggrepel \
    r-cowplot \
    r-stringr \
    r-scales \
    r-svglite
fi

# ── 3. R package verification ─────────────────────────────────────────────────
echo "[R] Verifying installed R packages …"
"${MAMBA}" run -n "${ENV_NAME}" \
  Rscript --no-save --no-restore --no-site-file --no-init-file -e "
  cat('\n=== R package verification ===\n')
  pkgs <- c(
    'optparse', 'dplyr', 'tidyr', 'tibble',
    'ggplot2', 'ggrepel', 'cowplot',
    'stringr', 'scales', 'svglite', 'parallel'
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

echo ""
echo "════════════════════════════════════════════════════════"
echo " Setup complete."
echo " Run the analysis (no env activation needed):"
echo "   Rscript ${SCRIPT_DIR}/run_depletion_analysis.R \\"
echo "     --crumblr  /path/to/crumblr_results.csv \\"
echo "     --manifest /path/to/manifest.csv \\"
echo "     --output-dir /path/to/output/"
echo "════════════════════════════════════════════════════════"
