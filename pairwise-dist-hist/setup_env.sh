#!/usr/bin/env bash
# setup_env.sh
# One-time setup: creates a Python venv (adata-venv) with anndata and
# companion packages needed for pairwise_dist_hist.py.
# Everything is confined to cc-sandbox.
#
# Usage (from any directory):
#   bash /path/to/pairwise-dist-hist/setup_env.sh

set -euo pipefail

SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"
VENV_DIR="${SANDBOX}/adata-venv"

echo "════════════════════════════════════════════════════════"
echo " pairwise-dist-hist — environment setup"
echo " Sandbox : ${SANDBOX}"
echo " Venv    : ${VENV_DIR}"
echo "════════════════════════════════════════════════════════"

# ── 1. Locate a python3 with the venv module ──────────────────────────────────
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
echo "[Python] Using: ${PYTHON3} ($(${PYTHON3} --version))"

# ── 2. Create venv ────────────────────────────────────────────────────────────
if [ -f "${VENV_DIR}/bin/activate" ]; then
  echo "[Python] venv already exists at ${VENV_DIR} — skipping creation."
else
  echo "[Python] Creating venv at ${VENV_DIR} …"
  "${PYTHON3}" -m venv --without-pip "${VENV_DIR}"
  curl -sS https://bootstrap.pypa.io/get-pip.py | "${VENV_DIR}/bin/python3"
fi

# ── 3. Install packages (only if any are missing) ────────────────────────────
VENV_PY="${VENV_DIR}/bin/python3"

# Map: import name → install spec
declare -A PKGS=(
  [anndata]="anndata>=0.10"
  [numpy]="numpy>=1.24"
  [pandas]="pandas>=2.0"
  [scipy]="scipy>=1.10"
  [sklearn]="scikit-learn>=1.3"
  [matplotlib]="matplotlib>=3.7"
  [seaborn]="seaborn>=0.13"
  [h5py]="h5py>=3.9"
  [diptest]="diptest>=0.10"
)

MISSING=()
for import_name in "${!PKGS[@]}"; do
  if ! "${VENV_PY}" -c "import ${import_name}" 2>/dev/null; then
    MISSING+=("${PKGS[${import_name}]}")
    echo "[Python] Missing: ${import_name} (${PKGS[${import_name}]})"
  fi
done

if [ ${#MISSING[@]} -eq 0 ]; then
  echo "[Python] All packages already installed — nothing to do."
else
  echo "[Python] Installing ${#MISSING[@]} missing package(s) …"
  source "${VENV_DIR}/bin/activate"
  pip install --quiet --upgrade pip
  pip install --quiet "${MISSING[@]}"
  deactivate
fi

echo ""
echo "════════════════════════════════════════════════════════"
echo " Setup complete."
echo " Run analysis (venv activates automatically):"
echo "   python ${SANDBOX}/pairwise-dist-hist/pairwise_dist_hist.py \\"
echo "     --h5ad  data.h5ad \\"
echo "     --obs-col cell_type \\"
echo "     --obsm-key X_pca \\"
echo "     --output  pairwise_distances.pdf"
echo "════════════════════════════════════════════════════════"
