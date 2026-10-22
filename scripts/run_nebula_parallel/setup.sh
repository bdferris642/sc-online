#!/usr/bin/env bash
# setup.sh
# One-time setup for the run_nebula_parallel workflow: creates the micromamba env
# run-nebula-parallel-env (R 4.5 + Seurat/qs/nebula 1.5.5, Python with h5py/anndata) inside
# cc-sandbox, then verifies every R and Python package the workflow loads.
#
# nebula is not on conda-forge: it is installed explicitly at version 1.5.5 (the version the
# existing sn-vta DE results were fit with) from its GitHub tag, compiled with the env's toolchain.
#
# Usage (from any directory):
#   bash /path/to/cc-sandbox/run_nebula_parallel/setup.sh

set -euo pipefail

SANDBOX="$(cd "$(dirname "$0")/.." && pwd)"
BIN_DIR="${SANDBOX}/bin"
MAMBA_ROOT="${SANDBOX}/micromamba_root"
ENV_NAME="run-nebula-parallel-env"
ENV_PREFIX="${MAMBA_ROOT}/envs/${ENV_NAME}"
NEBULA_VERSION="1.5.5"
# 1.5.5 was never released on CRAN (1.5.3 -> 1.5.6); it exists only as a GitHub tag
NEBULA_URL="https://github.com/lhe17/nebula/archive/refs/tags/v${NEBULA_VERSION}.tar.gz"

export MAMBA_ROOT_PREFIX="${MAMBA_ROOT}"
export CONDA_PKGS_DIRS="${MAMBA_ROOT}/conda-pkgs"   # persistent, same FS (never /tmp)
MAMBA="${BIN_DIR}/micromamba"

echo "════════════════════════════════════════════════════════"
echo " run_nebula_parallel — environment setup"
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

# ── 2. conda environment ──────────────────────────────────────────────────────
if [ -x "${ENV_PREFIX}/bin/Rscript" ]; then
  echo "[conda] Environment '${ENV_NAME}' already exists — ensuring packages."
  "${MAMBA}" install -n "${ENV_NAME}" -y -c conda-forge r-data.table r-jsonlite
else
  echo "[conda] Creating environment '${ENV_NAME}' …"
  "${MAMBA}" create -n "${ENV_NAME}" -y -c conda-forge -c bioconda \
    "r-base=4.5" \
    "r-seurat>=5.5" r-seuratobject "r-qs>=0.27" \
    r-matrix r-dplyr r-getopt r-glue r-future r-data.table r-jsonlite \
    r-rcpp r-rcppeigen r-nloptr r-rfast r-trust r-parallelly r-dofuture r-foreach r-dorng \
    bioconductor-singlecellexperiment \
    "python=3.12" numpy scipy pandas h5py anndata \
    parallel make gxx_linux-64 gcc_linux-64 gfortran_linux-64 curl
fi

# ── 3. nebula (CRAN archive, pinned) ─────────────────────────────────────────
installed_nebula="$("${ENV_PREFIX}/bin/Rscript" -e 'cat(if (requireNamespace("nebula", quietly=TRUE)) as.character(packageVersion("nebula")) else "none")' 2>/dev/null)"
if [ "${installed_nebula}" == "${NEBULA_VERSION}" ]; then
  echo "[nebula] ${NEBULA_VERSION} already installed"
else
  echo "[nebula] Installing ${NEBULA_VERSION} (found: ${installed_nebula})"
  tmp="$(mktemp -d "${MAMBA_ROOT}/nebula_build.XXXXXX")"
  curl -fsSL -o "${tmp}/nebula_${NEBULA_VERSION}.tar.gz" "${NEBULA_URL}" \
    || { echo "[nebula] ERROR: download failed: ${NEBULA_URL}"; exit 1; }
  tar -xzf "${tmp}/nebula_${NEBULA_VERSION}.tar.gz" -C "${tmp}"
  src="$(find "${tmp}" -maxdepth 2 -name DESCRIPTION -printf '%h\n' | head -1)"
  grep -q "^Version: ${NEBULA_VERSION}$" "${src}/DESCRIPTION" \
    || { echo "[nebula] ERROR: ${NEBULA_URL} is not version ${NEBULA_VERSION}"; exit 1; }
  PATH="${ENV_PREFIX}/bin:${PATH}" "${ENV_PREFIX}/bin/R" CMD INSTALL "${src}"
  rm -rf "${tmp}"
fi

# ── 4. verification ──────────────────────────────────────────────────────────
echo "[verify] R packages"
"${ENV_PREFIX}/bin/Rscript" -e '
pkgs <- c("Seurat", "SeuratObject", "qs", "Matrix", "dplyr", "getopt", "glue", "future",
          "data.table", "jsonlite", "nebula")
for (p in pkgs) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat(sprintf("  %-14s %s\n", p, as.character(packageVersion(p))))
}
if (as.character(packageVersion("nebula")) != "'"${NEBULA_VERSION}"'") stop("nebula version mismatch")
cat("  R", R.version.string, "\n")
'
echo "[verify] Python packages"
LD_LIBRARY_PATH="${ENV_PREFIX}/lib" "${ENV_PREFIX}/bin/python" -c '
import importlib.metadata as md
for p in ["numpy", "scipy", "pandas", "h5py", "anndata"]:
    __import__(p)
    print(f"  {p:<14} {md.version(p)}")
'
command -v "${ENV_PREFIX}/bin/parallel" > /dev/null || { echo "GNU parallel missing in env"; exit 1; }
echo "  parallel       $("${ENV_PREFIX}/bin/parallel" --version | head -1)"
echo "[done] ${ENV_NAME} ready at ${ENV_PREFIX}"
