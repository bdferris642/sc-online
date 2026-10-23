#!/usr/bin/env python
"""Control-cell PCA for an h5ad, projected onto all cells, with donor-level PC selection.

Subcommands
-----------
compute   Fit the PCA on control cells, project every cell, test each PC for case/control at the
          donor level, and (unless --no-write) replace the h5ad's control PCs in place.
select    Print the PCs to add as NEBULA covariates (one line, comma-separated), read from the
          h5ad's uns/control_pca: the union of the requested sets (elbow, significant) with
          index <= --pc-max. Also prints the uns/control_pca 'created' stamp on a second line.

compute, step by step (memory-minimal: CSR row blocks streamed in parallel, nothing dense
larger than (block rows x n_hvgs) is ever built, and the PCA itself is done from a
n_hvgs x n_hvgs Gram matrix):
  1. Control cells = <contrast-col> == --control-val. Each row of the counts layer is
     normalised to 1e4 by its total over ALL genes, then
     pass 1 accumulates per-gene detection count, sum and sum of squares over control cells.
     Genes detected in > 1% of control cells are eligible; HVGs = top --n-hvgs by scanpy's
     flavor="seurat" normalised dispersion (computed exactly from those moments).
  2. Pass 2 over control cells: log1p(normalised) restricted to the HVGs (sparse), accumulating
     the column sums, sums of squares and X^T X. Scaling uses mean / sd (ddof=0, sd + 1e-8); the
     top --n-pcs eigenvectors of the scaled Gram matrix are the loadings (identical to PCA /
     TruncatedSVD on the scaled dense matrix). Signs are fixed so each loading vector's
     largest-|value| entry is positive. The elbow (pick_pcs) is taken on variance explained over
     the fitted PCs. The model is saved to <pca-model-dir>/<slogan>_control_pca.npz (+ .json).
  3. Pass 3 projects ALL cells with the saved model:
     PC = X_hvg @ (L / sd) - (mean / sd) @ L, with X_hvg sparse log1p(normalised) counts.
     Columns are named control_PC_01 .. control_PC_<n>.
  4. Donor-level test per PC: OLS of the donor-mean PC on the contrast (case vs control) plus
     donor-level covariates (--pc-test-covs; numeric = donor mean, categorical = per-donor level
     fractions with the most common level as reference). Donors with < --min-cells-per-donor
     cells are dropped, covariates constant across donors are dropped (logged). BH FDR over the
     PCs; FDR < --pc-fdr = "significant". Results: <pca-model-dir>/<slogan>_control_pc_test.csv.
  5. In-place h5ad update (skipped with --no-write). Only obs control_PC* datasets and
     uns/control_pca change:
       - backup first: the current obs group and uns/control_pca are copied with h5py into
         <h5ad dir>/backup_<timestamp>/<slogan>_obs_uns.h5 (and the metadata obs CSV, if any);
       - new PC columns and the new uns/control_pca are written to temporary root groups and
         read back; then every existing obs/control_PC* dataset (old unpadded or padded) is
         deleted, the new ones moved in, obs 'column-order' updated, uns/control_pca swapped;
       - everything else (X, layers, obsm, other obs columns, other uns entries) is
         fingerprinted before and checked after; leftover temporary groups are an error;
       - <h5ad dir>/metadata/<slogan>_obs.csv, if present, gets the same column change.

The fit uses full-row totals for normalisation in both fitting and projection. (The previous
clean_h5ad/compute_pcs.py normalised control cells by their total over >1%-expressed genes
during fitting but by the full total when projecting.)
"""

import argparse
import datetime
import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from h5ad_io import ensure_env_python  # noqa: E402

PC_PREFIX = "control_PC_"
TMP_PCS = "__control_pcs_tmp__"
TMP_UNS = "__uns_control_pca_tmp__"
MIN_FREE_GB = 10
TARGET_SUM = 1e4
MIN_DETECT_FRAC = 0.01
N_MEAN_BINS = 20
DEFAULT_TEST_COVS = "age,sex,log10_nUMI,frac_mito,frac_intronic,chemistry"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def pc_names(n):
    """Zero-padded PC column names: control_PC_01 .. (width >= 2)."""
    width = max(2, len(str(n)))
    return [f"{PC_PREFIX}{i + 1:0{width}d}" for i in range(n)]


def pick_pcs(singular_values):
    """Elbow index (0-based) of the explained-variance curve (max distance to the chord).

    Copy of clean_h5ad/compute_pcs.py::_pick_pcs, kept here so the workflow is self-contained.
    """
    var = singular_values ** 2
    var_expl = var / var.sum()
    n = len(var_expl)
    if n <= 2:
        return 0
    x = np.linspace(0, 1, n)
    y = (var_expl - var_expl.min()) / (var_expl.max() - var_expl.min() + 1e-12)
    dx, dy = x[-1] - x[0], y[-1] - y[0]
    norm = np.sqrt(dx ** 2 + dy ** 2)
    dist = np.abs(dy * x - dx * y + x[-1] * y[0] - y[-1] * x[0]) / norm
    return int(np.argmax(dist))


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


# ───────────────────────────── block workers ─────────────────────────────

def _normalised_block(h5ad, layer, r0, r1):
    """Rows [r0, r1) of the counts layer -> (indptr, indices, values normalised to TARGET_SUM, row totals)."""
    from h5ad_io import read_csr_block
    indptr, indices, data = read_csr_block(h5ad, layer, r0, r1)
    data = data.astype(np.float64)
    lens = np.diff(indptr)
    totals = np.bincount(np.repeat(np.arange(r1 - r0), lens), weights=data, minlength=r1 - r0)
    safe = np.where(totals > 0, totals, 1.0)
    norm = data * np.repeat(TARGET_SUM / safe, lens)
    return indptr, indices, norm, totals


def _pass1(args):
    """Per-gene detection count, sum, sum of squares of normalised counts over control rows; row totals."""
    h5ad, layer, r0, r1, is_ctrl, n_var = args
    indptr, indices, norm, totals = _normalised_block(h5ad, layer, r0, r1)
    keep = np.repeat(is_ctrl, np.diff(indptr))
    gi = indices[keep]
    v = norm[keep]
    return (np.bincount(gi, minlength=n_var), np.bincount(gi, weights=v, minlength=n_var),
            np.bincount(gi, weights=v * v, minlength=n_var), totals)


def _hvg_csr(indptr, indices, norm, hvg_pos, rows_mask=None):
    """Sparse (rows x n_hvg) log1p(normalised) matrix restricted to HVG columns (optionally to some rows)."""
    import scipy.sparse as sp
    n_rows = len(indptr) - 1
    row_of = np.repeat(np.arange(n_rows), np.diff(indptr))
    col = hvg_pos[indices]
    keep = col >= 0
    if rows_mask is not None:
        keep &= rows_mask[row_of]
    m = sp.csr_matrix((np.log1p(norm[keep]), (row_of[keep], col[keep])),
                      shape=(n_rows, int(hvg_pos.max()) + 1))
    if rows_mask is not None:
        m = m[np.flatnonzero(rows_mask)]
    return m


def _pass2(args):
    """Column sums, sums of squares and Gram matrix of log1p HVG values over control rows."""
    h5ad, layer, r0, r1, is_ctrl, hvg_pos = args
    if not is_ctrl.any():
        n = int(hvg_pos.max()) + 1
        return np.zeros(n), np.zeros(n), np.zeros((n, n)), 0
    indptr, indices, norm, _ = _normalised_block(h5ad, layer, r0, r1)
    m = _hvg_csr(indptr, indices, norm, hvg_pos, is_ctrl)
    s = np.asarray(m.sum(axis=0)).ravel()
    sq = np.asarray(m.multiply(m).sum(axis=0)).ravel()
    gram = (m.T @ m).toarray()
    return s, sq, gram, m.shape[0]


def _pass3(args):
    """Project all rows of the block: X_hvg @ (L / sd) - (mean / sd) @ L."""
    h5ad, layer, r0, r1, hvg_pos, w, shift = args
    indptr, indices, norm, _ = _normalised_block(h5ad, layer, r0, r1)
    m = _hvg_csr(indptr, indices, norm, hvg_pos)
    return r0, (np.asarray(m @ w) - shift).astype(np.float32)


# ───────────────────────────── HVG / PCA ─────────────────────────────

def seurat_hvgs(n_ctrl, det, s, sq, n_top):
    """scanpy flavor='seurat' HVGs from moments of normalised counts over control cells.

    Mirrors scanpy.pp.highly_variable_genes(flavor='seurat') on log1p-normalised data:
    mean/var (ddof=1) of expm1(X) = normalised counts, dispersion = var/mean, log(dispersion),
    log1p(mean), 20 equal-width mean bins, z-scored dispersion within bin (single-gene bins get
    0 mean / own value as sd, as scanpy), top n by normalised dispersion (ties included).
    Only genes detected in > 1% of control cells are considered.
    Returns (hvg indices into var, eligible mask, table dict).
    """
    import pandas as pd
    eligible = det / n_ctrl > MIN_DETECT_FRAC
    idx = np.flatnonzero(eligible)
    mean = s[idx] / n_ctrl
    var = (sq[idx] / n_ctrl - mean ** 2) * n_ctrl / (n_ctrl - 1)
    mean[mean == 0] = 1e-12
    disp = var / mean
    disp[disp == 0] = np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        disp = np.log(disp)
    lmean = np.log1p(mean)
    df = pd.DataFrame({"means": lmean, "dispersions": disp})
    df["mean_bin"] = pd.cut(df["means"], bins=N_MEAN_BINS)
    grp = df.groupby("mean_bin", observed=True)["dispersions"]
    bmean, bsd = grp.mean(), grp.std(ddof=1)
    one = bsd.isnull()
    bsd[one.values] = bmean[one.values].values
    bmean[one.values] = 0
    dnorm = ((df["dispersions"].to_numpy() - bmean.reindex(df["mean_bin"]).to_numpy())
             / bsd.reindex(df["mean_bin"]).to_numpy())
    n_top = min(n_top, len(idx))
    finite = dnorm[~np.isnan(dnorm)]
    cutoff = np.sort(finite)[::-1][min(n_top, finite.size) - 1]
    chosen = np.nan_to_num(dnorm, nan=-np.inf) >= cutoff
    return idx[chosen], eligible


def fit_pca(n_ctrl, s, sq, gram, n_pcs):
    """Scaled-data PCA from the Gram matrix. Returns (mean, sd, loadings (h x k), singular values)."""
    mean = s / n_ctrl
    sd = np.sqrt(np.maximum(sq / n_ctrl - mean ** 2, 0)) + 1e-8
    cov = gram - n_ctrl * np.outer(mean, mean)          # centred X^T X
    z = cov / np.outer(sd, sd)                         # scaled X^T X
    evals, evecs = np.linalg.eigh(z)
    order = np.argsort(evals)[::-1][:n_pcs]
    evals, load = np.clip(evals[order], 0, None), evecs[:, order]
    flip = np.sign(load[np.abs(load).argmax(axis=0), np.arange(load.shape[1])])
    load = load * flip
    return mean, sd, load, np.sqrt(evals)


# ───────────────────────────── donor test ─────────────────────────────

def donor_test(obs, pcs, id_col, contrast_col, control_val, case_val, test_covs, min_cells):
    """Per-PC OLS on donor means. Returns (results DataFrame, model string, dropped covariates, n donors)."""
    import pandas as pd
    from scipy import stats
    from h5ad_io import r_type_of

    need = [id_col, contrast_col, *test_covs]
    for c in need:
        if obs[c].isna().any():
            raise ValueError(f"PC test: {int(obs[c].isna().sum())} missing values in obs '{c}'")
    donor = obs[id_col].astype(str).to_numpy()
    donors, codes = np.unique(donor, return_inverse=True)
    n_cells = np.bincount(codes)
    cc = obs[contrast_col].astype(str).to_numpy()
    bad = sorted(set(cc) - {control_val, case_val})
    if bad:
        raise ValueError(f"PC test: '{contrast_col}' has levels other than {control_val}/{case_val}: {bad}")
    case_frac = np.bincount(codes, weights=(cc == case_val).astype(float)) / n_cells
    if not np.isin(case_frac, [0.0, 1.0]).all():
        raise ValueError(f"PC test: '{contrast_col}' varies within donors {list(donors[(case_frac > 0) & (case_frac < 1)][:10])}")
    keep = n_cells >= min_cells
    if (~keep).any():
        log(f"PC test: dropping {int((~keep).sum())} donors with < {min_cells} cells")

    cols, names, dropped = [case_frac], [f"{contrast_col}{case_val}"], []
    for c in test_covs:
        s = obs[c]
        if r_type_of(s) == "numeric":
            v = np.bincount(codes, weights=s.to_numpy(dtype=float)) / n_cells
            if np.ptp(v[keep]) == 0:
                dropped.append(c)
                continue
            cols.append(v); names.append(c)
        else:
            lv = s.astype(str).to_numpy()
            levels, lc = np.unique(lv, return_inverse=True)
            tab = np.zeros((len(donors), len(levels)))
            np.add.at(tab, (codes, lc), 1.0)
            frac = tab / n_cells[:, None]
            ref = int(tab[keep].sum(axis=0).argmax())
            added = 0
            for j, l in enumerate(levels):
                if j == ref or np.ptp(frac[keep, j]) == 0:
                    continue
                cols.append(frac[:, j]); names.append(f"{c}{l}"); added += 1
            if not added:
                dropped.append(c)
    for c in dropped:
        log(f"PC test: covariate '{c}' is constant across donors -> dropped")

    x = np.column_stack([np.ones(len(donors)), *cols])[keep]
    n, p = x.shape
    if np.linalg.matrix_rank(x) < p:
        raise ValueError(f"PC test: donor design is rank deficient ({p} columns: intercept + {names})")
    if n - p < 3:
        raise ValueError(f"PC test: only {n} donors for {p} parameters")
    if case_frac[keep].min() == case_frac[keep].max():
        raise ValueError("PC test: only one contrast group among kept donors")
    y = np.column_stack([np.bincount(codes, weights=pcs[:, k].astype(float)) / n_cells
                         for k in range(pcs.shape[1])])[keep]
    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    resid = y - x @ beta
    df = n - p
    sigma2 = (resid ** 2).sum(axis=0) / df
    se = np.sqrt(sigma2 * xtx_inv[1, 1])
    t = beta[1] / se
    pval = 2 * stats.t.sf(np.abs(t), df)
    res = pd.DataFrame({"pc": pc_names(pcs.shape[1]), "beta": beta[1], "se": se, "t": t, "df": df,
                        "p": pval, "fdr": bh_fdr(pval), "n_donors": n,
                        "n_case_donors": int(case_frac[keep].sum()),
                        "n_control_donors": int(n - case_frac[keep].sum())})
    model = f"donor_mean(PC) ~ {' + '.join(names)}  [OLS, donors with >= {min_cells} cells]"
    return res, model, dropped, n


# ───────────────────────────── h5ad writing ─────────────────────────────

def _fingerprint(h, pc_cols):
    """Structure of everything the in-place update must not change."""
    import h5py
    fp = {"root": sorted(k for k in h.keys() if k not in (TMP_PCS, TMP_UNS))}
    for g in ("layers", "obsm", "varm", "obsp", "varp", "var"):
        fp[g] = sorted(h[g].keys()) if g in h else None
    for k in ("X", "layers/counts", "layers/norm_log"):
        if k in h and isinstance(h[k], type(h["obs"])) and "indptr" in h[k]:
            fp[k] = (int(h[k]["indptr"].shape[0]), int(h[k]["indptr"][-1]), str(h[k]["data"].dtype))
    obs = h["obs"]
    # on-disk object address: stable across opens, changes only if the element is rewritten
    fp["obs_other"] = {k: h5py.h5o.get_info(obs[k].id).addr
                       for k in obs.keys() if k not in pc_cols and not k.startswith("control_PC")}
    fp["obs_index"] = obs.attrs["_index"]
    fp["uns_other"] = sorted(k for k in h["uns"].keys() if k != "control_pca")
    return fp


def _old_pc_cols(obs_group):
    return [k for k in obs_group.keys() if k.startswith("control_PC")]


def write_in_place(h5ad, pcs, uns_new, backup=True, metadata_csv=None):
    """Replace obs control_PC* and uns/control_pca in place (backup first). See module docstring."""
    import anndata as ad
    import h5py
    import pandas as pd

    h5ad = Path(h5ad)
    names = pc_names(pcs.shape[1])
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    slogan = h5ad.name[:-len(".h5ad")]
    backup_dir = h5ad.parent / f"backup_{stamp}"

    with h5py.File(h5ad, "r") as h:
        left = [t for t in (TMP_PCS, TMP_UNS) if t in h]
        if left:
            raise RuntimeError(f"{h5ad}: leftover temporary groups {left} from an aborted run (inspect, then delete)")
        n_obs = len(ad.io.read_elem(h["obs"][h["obs"].attrs["_index"]]))
        if n_obs != pcs.shape[0]:
            raise ValueError(f"PC rows {pcs.shape[0]} != n_obs {n_obs}")
        old_cols = _old_pc_cols(h["obs"])

    if backup:
        free_gb = shutil.disk_usage(h5ad.parent).free / 1e9
        if free_gb < MIN_FREE_GB + 5:
            raise RuntimeError(f"only {free_gb:.1f} GB free at {h5ad.parent}; need > {MIN_FREE_GB + 5} for the backup")
        backup_dir.mkdir()
        bpath = backup_dir / f"{slogan}_obs_uns.h5"
        with h5py.File(h5ad, "r") as src, h5py.File(bpath, "w") as dst:
            src.copy(src["obs"], dst, name="obs")
            if "control_pca" in src["uns"]:
                src.copy(src["uns/control_pca"], dst, name="uns_control_pca")
            dst.attrs["source"] = str(h5ad)
            dst.attrs["restore"] = ("h5py: del f['obs']; f.copy(b['obs'], f, 'obs'); "
                                    "del f['uns/control_pca']; f.copy(b['uns_control_pca'], f['uns'], 'control_pca')")
        log(f"backup of obs + uns/control_pca -> {bpath} ({bpath.stat().st_size / 1e6:.1f} MB)")
        if metadata_csv is not None and Path(metadata_csv).exists():
            shutil.copy2(metadata_csv, backup_dir / Path(metadata_csv).name)
            log(f"backup of {metadata_csv} -> {backup_dir}")

    with h5py.File(h5ad, "r+") as h:
        fp_before = _fingerprint(h, names)
        tmp = h.create_group(TMP_PCS)
        for k, nm in enumerate(names):
            ad.io.write_elem(tmp, nm, np.ascontiguousarray(pcs[:, k]))
        ad.io.write_elem(h, TMP_UNS, uns_new)
        for k, nm in enumerate(names):
            if not np.array_equal(ad.io.read_elem(tmp[nm]), pcs[:, k]):
                raise RuntimeError(f"read-back of {nm} differs")
        if sorted(ad.io.read_elem(h[TMP_UNS]).keys()) != sorted(uns_new.keys()):
            raise RuntimeError("read-back of uns/control_pca keys differs")
        # swap
        obs = h["obs"]
        order = [str(c) for c in obs.attrs["column-order"]]
        for c in old_cols:
            del obs[c]
        for nm in names:
            h.move(f"{TMP_PCS}/{nm}", f"obs/{nm}")
        new_order = [c for c in order if not c.startswith("control_PC")] + names
        obs.attrs["column-order"] = np.array(new_order, dtype=h5py.string_dtype())
        del h[TMP_PCS]
        if "control_pca" in h["uns"]:
            del h["uns/control_pca"]
        h.move(TMP_UNS, "uns/control_pca")

    with h5py.File(h5ad, "r") as h:
        fp_after = _fingerprint(h, names)
        if fp_after != fp_before:
            diff = {k for k in fp_before if fp_before[k] != fp_after.get(k)}
            raise RuntimeError(f"{h5ad}: elements outside control PCs changed: {diff}")
        obs_after = ad.io.read_elem(h["obs"])
        cp_after = ad.io.read_elem(h["uns/control_pca"])
    got = [c for c in obs_after.columns if c.startswith("control_PC")]
    if got != names or not np.array_equal(obs_after[names].to_numpy(), pcs):
        raise RuntimeError(f"post-write check failed: obs PC columns {got[:5]}...")
    if sorted(cp_after.keys()) != sorted(uns_new.keys()):
        raise RuntimeError("post-write check failed: uns/control_pca keys")
    log(f"h5ad updated in place: removed {len(old_cols)} old control_PC* columns, wrote {len(names)} "
        f"({names[0]}..{names[-1]}); uns/control_pca replaced; other elements unchanged")

    if metadata_csv is not None and Path(metadata_csv).exists():
        csv = Path(metadata_csv)
        # round_trip: the default float parser is not exact and would perturb untouched columns
        df = pd.read_csv(csv, index_col=0, low_memory=False, float_precision="round_trip",
                         dtype={"participant_id": str})
        df.index = df.index.astype(str)
        if not df.index.equals(obs_after.index.astype(str)):
            raise ValueError(f"{csv} index does not match the h5ad obs index; not rewritten")
        df = df[[c for c in df.columns if not c.startswith("control_PC")]]
        df = pd.concat([df, obs_after[names]], axis=1)
        tmp_csv = csv.with_suffix(".csv.new")
        df.to_csv(tmp_csv)
        tmp_csv.replace(csv)
        log(f"rewrote {csv} with the new control PCs")
    return backup_dir if backup else None


# ───────────────────────────── commands ─────────────────────────────

def _uns_payload(names, var_expl, n_elbow, res, model, dropped, n_hvgs, args, model_path):
    sig = res.loc[res["fdr"] < args.pc_fdr, "pc"].tolist()
    return {
        "pcs": np.array(names, dtype=str),
        "variance_explained": np.asarray(var_expl, dtype=float),
        "n_elbow_pcs": int(n_elbow),
        "elbow_pcs": np.array(names[:n_elbow], dtype=str),
        "significant_pcs": np.array(sig, dtype=str),
        "fdr_threshold": float(args.pc_fdr),
        "test_method": "pseudobulk_ols",
        "model": model,
        "covariates_dropped": ";".join(dropped),
        "case_control_test": {k: res[k].to_numpy(dtype=float) for k in ("beta", "se", "t", "p", "fdr")}
                             | {"pc": res["pc"].to_numpy(dtype=str)},
        "n_hvgs": int(n_hvgs),
        "control_val": args.control_val,
        "case_val": args.case_val,
        "contrast_col": args.contrast_col,
        "model_path": str(model_path),
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
    }


def cmd_compute(args):
    import pandas as pd
    from h5ad_io import csr_info, plan_row_blocks, read_obs, read_var_names

    t0 = time.time()
    h5ad = Path(args.h5ad).resolve()
    slogan = h5ad.name[:-len(".h5ad")]
    model_dir = Path(args.pca_model_dir) if args.pca_model_dir else h5ad.parent / "pca_models"
    model_dir.mkdir(parents=True, exist_ok=True)
    test_covs = [c for c in args.pc_test_covs.split(",") if c]

    obs_cols = list(dict.fromkeys([args.id_col, args.contrast_col, *[c for c in test_covs if c != "log10_nUMI"]]))
    obs = read_obs(h5ad, obs_cols)
    is_ctrl = (obs[args.contrast_col].astype(str) == args.control_val).to_numpy()
    n_ctrl = int(is_ctrl.sum())
    if n_ctrl < 10:
        raise ValueError(f"only {n_ctrl} cells with {args.contrast_col} == '{args.control_val}'")
    n_obs, n_var, nnz, indptr = csr_info(h5ad, args.layer)
    if n_obs != len(obs):
        raise ValueError(f"layer rows {n_obs} != obs rows {len(obs)}")
    blocks = plan_row_blocks(indptr, int(args.block_nnz))
    log(f"{slogan}: {n_obs:,} cells ({n_ctrl:,} '{args.control_val}'), {n_var:,} genes, nnz {nnz:,}; "
        f"{len(blocks)} blocks, {args.workers} workers")
    ex = ProcessPoolExecutor(args.workers)

    # pass 1
    det, s, sq = np.zeros(n_var), np.zeros(n_var), np.zeros(n_var)
    totals = np.zeros(n_obs)
    for (r0, r1), (d_, s_, q_, t_) in zip(blocks, ex.map(
            _pass1, [(str(h5ad), args.layer, r0, r1, is_ctrl[r0:r1], n_var) for r0, r1 in blocks])):
        det += d_; s += s_; sq += q_; totals[r0:r1] = t_
    hvg_idx, eligible = seurat_hvgs(n_ctrl, det, s, sq, args.n_hvgs)
    log(f"pass 1: {int(eligible.sum()):,} genes detected in > 1% of control cells; {len(hvg_idx)} HVGs")
    hvg_pos = np.full(n_var, -1, dtype=np.int64)
    hvg_pos[hvg_idx] = np.arange(len(hvg_idx))

    # pass 2
    h = len(hvg_idx)
    s2, sq2, gram, n2 = np.zeros(h), np.zeros(h), np.zeros((h, h)), 0
    for a, b, g, n in ex.map(_pass2, [(str(h5ad), args.layer, r0, r1, is_ctrl[r0:r1], hvg_pos) for r0, r1 in blocks]):
        s2 += a; sq2 += b; gram += g; n2 += n
    if n2 != n_ctrl:
        raise RuntimeError(f"pass 2 saw {n2} control rows, expected {n_ctrl}")
    n_pcs = min(args.n_pcs, h - 1, n_ctrl - 1)
    mean, sd, load, sv = fit_pca(n_ctrl, s2, sq2, gram, n_pcs)
    var_expl = sv ** 2 / (sv ** 2).sum()
    n_elbow = pick_pcs(sv) + 1
    names = pc_names(n_pcs)
    var_names = read_var_names(h5ad)
    model_path = model_dir / f"{slogan}_control_pca.npz"
    np.savez(model_path, hvgs=var_names[hvg_idx], hvg_index=hvg_idx, mean=mean, sd=sd, loadings=load,
             singular_values=sv, variance_explained=var_expl, n_elbow_pcs=n_elbow, pcs=np.array(names))
    meta = {"h5ad": str(h5ad), "layer": args.layer, "contrast_col": args.contrast_col,
            "control_val": args.control_val, "n_control_cells": n_ctrl, "n_hvgs": h,
            "n_pcs": n_pcs, "n_elbow_pcs": n_elbow, "variance_explained": var_expl.tolist(),
            "target_sum": TARGET_SUM, "min_detect_frac": MIN_DETECT_FRAC, "hvg_flavor": "seurat",
            "scaling": "mean/sd (ddof=0), sd + 1e-8", "created": datetime.datetime.now().isoformat(timespec="seconds")}
    model_path.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    log(f"pass 2: PCA fit on {n_ctrl:,} control cells x {h} HVGs; elbow at PC {n_elbow}; "
        f"variance explained PC1-3 {np.round(var_expl[:3], 3).tolist()}; model -> {model_path}")

    # pass 3 (projection with the saved model)
    mdl = np.load(model_path)
    w = mdl["loadings"] / mdl["sd"][:, None]
    shift = (mdl["mean"] / mdl["sd"]) @ mdl["loadings"]
    pcs = np.zeros((n_obs, n_pcs), dtype=np.float32)
    for r0, block in ex.map(_pass3, [(str(h5ad), args.layer, r0, r1, hvg_pos, w, shift) for r0, r1 in blocks]):
        pcs[r0:r0 + block.shape[0]] = block
    ex.shutdown()
    log(f"pass 3: projected {n_obs:,} cells onto {n_pcs} PCs")

    # donor-level test
    if "log10_nUMI" in test_covs:
        obs["log10_nUMI"] = np.log10(totals + 1)     # as run-nebula.R: log10(nCount_RNA + 1), all genes
    res, model, dropped, n_donors = donor_test(obs, pcs, args.id_col, args.contrast_col, args.control_val,
                                               args.case_val, test_covs, args.min_cells_per_donor)
    test_csv = model_dir / f"{slogan}_control_pc_test.csv"
    res.to_csv(test_csv, index=False)
    sig = res.loc[res["fdr"] < args.pc_fdr, "pc"].tolist()
    log(f"PC test ({n_donors} donors): {len(sig)} PCs with FDR < {args.pc_fdr}: {sig}; -> {test_csv}")
    log(f"elbow PCs: {names[:n_elbow]}")

    uns_new = _uns_payload(names, var_expl, n_elbow, res, model, dropped, h, args, model_path)
    if args.no_write:
        np.save(model_dir / f"{slogan}_control_pcs.npy", pcs)
        log(f"--no-write: h5ad untouched; PC scores -> {model_dir / (slogan + '_control_pcs.npy')}")
    else:
        csv = h5ad.parent / "metadata" / f"{slogan}_obs.csv"
        write_in_place(h5ad, pcs, uns_new, backup=not args.no_backup, metadata_csv=csv if csv.exists() else None)
    log(f"done in {(time.time() - t0) / 60:.1f} min")


def cmd_select(args):
    import anndata as ad
    import h5py
    with h5py.File(args.h5ad, "r") as h:
        if "uns" not in h or "control_pca" not in h["uns"]:
            raise KeyError(f"{args.h5ad} has no uns/control_pca: run with --compute-control-pcs 1")
        cp = ad.io.read_elem(h["uns/control_pca"])
        obs_cols = set(h["obs"].keys())
    if "created" not in cp or "pcs" not in cp or not str(cp["pcs"][0]).startswith(PC_PREFIX):
        raise ValueError(f"{args.h5ad}: uns/control_pca is in the old format (unpadded control_PC names, "
                         f"no 'created'); recompute with --compute-control-pcs 1")
    names = [str(x) for x in cp["pcs"]]
    chosen = set()
    for kind in args.pcs_to_add:
        if kind == "all":
            chosen |= set(names)
            continue
        key = {"elbow": "elbow_pcs", "significant": "significant_pcs"}[kind]
        chosen |= {str(x) for x in cp[key]}
    out = [n for i, n in enumerate(names) if n in chosen and i + 1 <= args.pc_max]
    missing = [n for n in out if n not in obs_cols]
    if missing:
        raise KeyError(f"PCs listed in uns/control_pca are missing from obs: {missing}")
    print(",".join(out))
    print(cp["created"])


def main():
    ensure_env_python()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute", help="fit, project, test and write control PCs")
    c.add_argument("--h5ad", required=True)
    c.add_argument("--layer", default="counts", help="CSR raw-count layer ('X' for adata.X)")
    c.add_argument("--id-col", default="participant_id")
    c.add_argument("--contrast-col", default="case_control")
    c.add_argument("--control-val", default="ctr")
    c.add_argument("--case-val", default="pd")
    c.add_argument("--n-hvgs", type=int, default=2000)
    c.add_argument("--n-pcs", type=int, default=30)
    c.add_argument("--pc-test-covs", default=DEFAULT_TEST_COVS,
                   help="donor-level OLS covariates; log10_nUMI is computed from the counts if requested")
    c.add_argument("--pc-fdr", type=float, default=0.05)
    c.add_argument("--min-cells-per-donor", type=int, default=5)
    c.add_argument("--pca-model-dir", default=None, help="default: <h5ad dir>/pca_models")
    c.add_argument("--workers", type=int, default=24)
    c.add_argument("--block-nnz", type=float, default=2e7)
    c.add_argument("--no-write", action="store_true", help="fit/test only; leave the h5ad untouched")
    c.add_argument("--no-backup", action="store_true", help="skip the obs/uns sidecar backup (not recommended)")
    s = sub.add_parser("select", help="print PCs to add from uns/control_pca")
    s.add_argument("--h5ad", required=True)
    s.add_argument("--pcs-to-add", nargs="+", choices=["elbow", "significant", "all"], required=True,
                   help="union of these sets; 'all' = every computed PC (control_PC_01..) up to --pc-max")
    s.add_argument("--pc-max", type=int, default=20)
    a = ap.parse_args()
    {"compute": cmd_compute, "select": cmd_select}[a.cmd](a)


if __name__ == "__main__":
    main()
