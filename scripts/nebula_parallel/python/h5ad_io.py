"""Shared helpers for streaming CSR layers and obs columns out of very large h5ad files.

Everything here works on raw h5py handles so that no step ever materialises the full
matrix (which can exceed 2^31 non-zeros and therefore any R dgCMatrix).
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

WF_DIR = Path(__file__).resolve().parent.parent          # run_nebula_parallel/
SANDBOX = WF_DIR.parent
ENV_BIN = SANDBOX / "micromamba_root" / "envs" / "run-nebula-parallel-env" / "bin"
INT32_MAX = 2**31 - 1


def ensure_env_python():
    """Re-exec the calling script under the workflow env python (run setup.sh first)."""
    env_python = ENV_BIN / "python"
    if Path(sys.prefix).resolve() != ENV_BIN.parent.resolve():
        import os
        if not env_python.exists():
            raise FileNotFoundError(f"{env_python} not found: run {WF_DIR / 'setup.sh'}")
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = str(ENV_BIN.parent / "lib")
        os.execve(str(env_python), [str(env_python), *sys.argv], env)


def layer_key(layer):
    """Return the h5 group path for a layer name ('X' means adata.X)."""
    return "X" if layer == "X" else f"layers/{layer}"


def csr_info(h5ad_path, layer):
    """Return (n_obs, n_var, nnz, indptr) for a CSR layer; fail if the layer is not CSR."""
    with h5py.File(h5ad_path, "r") as f:
        key = layer_key(layer)
        if key not in f:
            raise KeyError(f"Layer '{layer}' ({key}) not found in {h5ad_path}")
        g = f[key]
        enc = g.attrs.get("encoding-type")
        if enc != "csr_matrix":
            raise ValueError(f"Layer '{layer}' must be csr_matrix (cells x genes); found encoding-type={enc!r}")
        n_obs, n_var = (int(x) for x in g.attrs["shape"])
        indptr = g["indptr"][:].astype(np.int64)
        nnz = int(indptr[-1])
    return n_obs, n_var, nnz, indptr


def plan_row_blocks(indptr, block_nnz):
    """Split rows into contiguous blocks of roughly `block_nnz` stored entries each.

    Returns a list of (row_start, row_end) half-open intervals covering all rows.
    """
    n_obs = len(indptr) - 1
    targets = np.arange(block_nnz, indptr[-1], block_nnz, dtype=np.int64)
    cuts = np.unique(np.searchsorted(indptr, targets, side="left"))
    cuts = cuts[(cuts > 0) & (cuts < n_obs)]
    bounds = np.concatenate([[0], cuts, [n_obs]])
    return [(int(a), int(b)) for a, b in zip(bounds[:-1], bounds[1:]) if b > a]


def read_csr_block(h5ad_path, layer, r0, r1, want_data=True):
    """Read rows [r0, r1) of a CSR layer.

    Returns (local_indptr int64, indices int64, data float32 or None).
    """
    with h5py.File(h5ad_path, "r") as f:
        g = f[layer_key(layer)]
        indptr = g["indptr"][r0:r1 + 1].astype(np.int64)
        s, e = int(indptr[0]), int(indptr[-1])
        indices = g["indices"][s:e].astype(np.int64, copy=False)
        data = g["data"][s:e] if want_data else None
    return indptr - s, indices, data


def read_var_names(h5ad_path):
    """Return the var index (feature names) as a numpy str array."""
    import anndata
    with h5py.File(h5ad_path, "r") as f:
        idx_key = f["var"].attrs["_index"]
        return np.asarray(anndata.io.read_elem(f["var"][idx_key])).astype(str)


def read_obs(h5ad_path, cols):
    """Read the obs index plus the requested obs columns into a DataFrame (fail on missing columns)."""
    import anndata
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        missing = [c for c in cols if c not in obs]
        if missing:
            raise KeyError(f"obs columns not found in {h5ad_path}: {missing}")
        idx_key = obs.attrs["_index"]
        index = np.asarray(anndata.io.read_elem(obs[idx_key])).astype(str)
        df = pd.DataFrame({c: anndata.io.read_elem(obs[c]) for c in cols}, index=index)
    if not df.index.is_unique:
        raise ValueError("obs index (cell names) is not unique")
    return df


REGION_DUMMIES = {"region_SN": "SN", "region_VTA": "VTA"}


def read_obs_derived(h5ad_path, cols, derive_region_dummies=False):
    """read_obs, optionally rebuilding region_SN / region_VTA from `region`.

    Coding: region_SN = region contains "SN", region_VTA = region contains "VTA", so
    SN -> (1, 0), VTA -> (0, 1), SN_VTA -> (1, 1), anything else -> (0, 0). Fails if the rebuilt values disagree with any
    non-missing stored value; logs how many missing values were filled. Needed for the
    oligodendrocyte object, whose spatial SN cells have region_VTA = NaN.
    """
    derived = [c for c in cols if derive_region_dummies and c in REGION_DUMMIES]
    if not derived:
        return read_obs(h5ad_path, cols)
    extra = ["region"] if "region" not in cols else []
    obs = read_obs(h5ad_path, list(cols) + extra)
    region = obs["region"].astype(str)
    for c in derived:
        new = region.str.contains(REGION_DUMMIES[c], regex=False).astype(np.int64)
        old = obs[c]
        known = old.notna()
        bad = int((old[known].astype(float) != new[known]).sum())
        if bad:
            raise ValueError(f"--derive-region-dummies: rebuilt {c} disagrees with {bad} stored non-missing values")
        print(f"[derive-region-dummies] {c} = (region contains '{REGION_DUMMIES[c]}'); "
              f"filled {int((~known).sum())} missing values, all stored values reproduced", flush=True)
        obs[c] = new
    return obs.drop(columns=extra)


def split_covariates(covs_raw):
    """Split a comma-separated covariate string; interaction terms ('a:b', 'a*b') expand to their parts.

    Mirrors chunk-sobj-by-genes.R / run-nebula.R so the same metadata columns are carried.
    Returns (formula_terms, base_columns).
    """
    terms = [c for c in (covs_raw or "").split(",") if c]
    base = []
    for t in terms:
        parts = t.replace("*", ":").split(":")
        base.extend(p for p in parts if p)
    return terms, list(dict.fromkeys(base))


def r_type_of(series):
    """Map a pandas column to the R type it must become: 'numeric', 'logical' or 'character'."""
    if pd.api.types.is_bool_dtype(series):
        return "logical"
    if isinstance(series.dtype, pd.CategoricalDtype) or not pd.api.types.is_numeric_dtype(series):
        return "character"
    return "numeric"
