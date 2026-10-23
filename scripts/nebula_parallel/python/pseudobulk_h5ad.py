#!/usr/bin/env python
"""Donor pseudobulk plus donor-level covariates, for rank purity.

Python replacement for steps 1-4 of compute_rank_purity() in ~/sc-online/add_rank_purity.R,
streamed from the h5ad so the >2^31-nnz matrix is never loaded whole.

Scale (--scale):
    log-rate (default)  log( r + eps_g ), r = mean over the donor's cells of counts / offset and
                        eps_g = 1% of the gene's mean r over donors (same for every donor; the raw
                        r is also saved as pb_rate.npy).
                        This is NEBULA's log-link scale (log E[y] = log offset + X b), so the fitted
                        covariate contributions B X that rank_purity.py subtracts are in matching
                        units. Subtracting log-scale NEBULA coefficients from a log1p-normalised
                        pseudobulk (scale 'mean' on norm_log) over-corrects lowly expressed genes and
                        shifted the oligo cases vs controls for most genes.
    mean                donor mean of --layer as stored (upstream behaviour on its 'data' slot).

Outputs in --out-dir:
    pb_mean.npy          donors x genes float64, donor pseudobulk on the chosen scale
    donors.txt, genes.txt
    donor_n_cells.tsv
    donor_numeric.tsv    per-donor mean (NaN cells ignored) of each numeric covariate
    donor_levels.tsv     long table (donor, covariate, level, fraction): per-donor fraction of
                         cells at each level of each categorical covariate (NaN cells ignored)
    covariate_types.json {covariate: "numeric" | "categorical"}

Unlike upstream (which takes the first cell's value for categoricals), categoricals that vary
within a donor (e.g. sequencing chemistry) get fractional values; for donor-constant
categoricals the fractions are 0/1 and identical to upstream.
"""

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from h5ad_io import (csr_info, ensure_env_python, plan_row_blocks, r_type_of, read_csr_block,  # noqa: E402
                     read_obs_derived, read_var_names, split_covariates)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def block_donor_sums(args):
    """Sum (weighted) rows of one CSR block into donors -> dense (n_donors x n_genes) float64.

    offset None: plain row sums (scale 'mean'). offset 'rowsum': each row is divided by its total
    over all genes (nCount_RNA); an array: divided by those per-row values (scale 'log-rate').
    """
    h5ad, layer, r0, r1, codes, n_donors, n_var, offset = args
    indptr, indices, data = read_csr_block(h5ad, layer, r0, r1)
    m = sp.csr_matrix((data.astype(np.float64), indices, indptr), shape=(r1 - r0, n_var))
    if offset is None:
        w = np.ones(r1 - r0)
    else:
        off = np.asarray(m.sum(axis=1)).ravel() if isinstance(offset, str) else offset
        if not (np.isfinite(off).all() and (off > 0).all()):
            raise ValueError(f"rows {r0}-{r1}: offset must be finite and > 0 (min {off.min()})")
        w = 1.0 / off
    a = sp.csr_matrix((w, (codes, np.arange(r1 - r0))), shape=(n_donors, r1 - r0))
    return (a @ m).toarray()


def level_strings(s):
    """Categorical/bool/string column -> object array of level strings as R would print them (NaN kept)."""
    if pd.api.types.is_bool_dtype(s):
        return s.map({True: "TRUE", False: "FALSE"}).astype(object)
    return s.astype(object).where(s.notna(), np.nan).map(lambda v: v if pd.isna(v) else str(v))


def main():
    ensure_env_python()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--layer", default="counts",
                    help="CSR layer ('X' for adata.X): raw counts for --scale log-rate, normalised for mean")
    ap.add_argument("--scale", choices=["log-rate", "mean"], default="log-rate",
                    help="log-rate: log of donor mean of counts/offset (NEBULA's link scale); "
                         "mean: donor mean of --layer as is")
    ap.add_argument("--offset-col", default="nCount_RNA",
                    help="log-rate offset: nCount_RNA = per-cell total of --layer over all genes; else obs column")
    ap.add_argument("--id-col", required=True)
    ap.add_argument("--covs", default="")
    ap.add_argument("--contrast-col", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--derive-region-dummies", action="store_true",
                    help="rebuild region_SN/region_VTA from region (SN->1,0; VTA->0,1; else 0,0)")
    ap.add_argument("--block-nnz", type=float, default=1e8)
    a = ap.parse_args()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    _, base_covs = split_covariates(a.covs)
    cov_cols = list(dict.fromkeys([*base_covs, a.contrast_col]))
    cov_cols = [c for c in cov_cols if c != a.id_col]
    computed = [c for c in cov_cols if c == "nCount_RNA"]
    if computed:
        raise ValueError("nCount_RNA as a rank-purity covariate is not supported (it is not an obs column)")
    obs = read_obs_derived(a.h5ad, [a.id_col, *cov_cols], a.derive_region_dummies)
    if obs[a.id_col].isna().any():
        raise ValueError(f"id column '{a.id_col}' has missing values")

    donor_str = obs[a.id_col].astype(str).to_numpy()
    donors, first_idx, codes = np.unique(donor_str, return_index=True, return_inverse=True)
    # donor order = first appearance in obs (as upstream's unique())
    order = np.argsort(first_idx)
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    donors = donors[order]
    codes = remap[codes].astype(np.int64)
    n_donors = len(donors)
    n_cells_per_donor = np.bincount(codes, minlength=n_donors)
    log(f"{len(obs):,} cells, {n_donors} donors")

    # ---- donor-level covariates
    types, numeric_cols, level_rows, varying = {}, {}, [], []
    for c in cov_cols:
        s = obs[c]
        if r_type_of(s) == "numeric":
            types[c] = "numeric"
            v = s.to_numpy(dtype=np.float64)
            ok = np.isfinite(v)
            sums = np.bincount(codes[ok], weights=v[ok], minlength=n_donors)
            cnt = np.bincount(codes[ok], minlength=n_donors)
            with np.errstate(invalid="ignore", divide="ignore"):
                numeric_cols[c] = np.where(cnt > 0, sums / cnt, np.nan)
            continue
        types[c] = "categorical"
        lv = level_strings(s)
        ok = lv.notna().to_numpy()
        levels, lcodes = np.unique(lv[ok].to_numpy().astype(str), return_inverse=True)
        tab = np.zeros((n_donors, len(levels)))
        np.add.at(tab, (codes[ok], lcodes), 1.0)
        tot = tab.sum(axis=1, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where(tot > 0, tab / tot, np.nan)
        mixed = donors[(tab > 0).sum(axis=1) > 1]
        if len(mixed):
            varying.append((c, mixed))
        for di, d in enumerate(donors):
            for li, lvl in enumerate(levels):
                level_rows.append((d, c, lvl, frac[di, li]))

    for c, mixed in varying:
        if c == a.contrast_col:
            raise ValueError(f"contrast column '{c}' varies within {len(mixed)} donors, e.g. {list(mixed[:10])}")
        log(f"categorical '{c}' varies within {len(mixed)}/{n_donors} donors -> fractional dummies")

    pd.DataFrame({"donor": donors, "n_cells": n_cells_per_donor}).to_csv(out / "donor_n_cells.tsv", sep="\t", index=False)
    pd.DataFrame({"donor": donors, **numeric_cols}).to_csv(out / "donor_numeric.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(level_rows, columns=["donor", "covariate", "level", "fraction"]).to_csv(
        out / "donor_levels.tsv", sep="\t", index=False, na_rep="NA")
    with open(out / "covariate_types.json", "w") as fh:
        json.dump(types, fh, indent=2)

    # ---- pseudobulk
    n_obs, n_var, nnz, indptr = csr_info(a.h5ad, a.layer)
    if n_obs != len(obs):
        raise ValueError(f"layer has {n_obs} rows but obs has {len(obs)}")
    blocks = plan_row_blocks(indptr, int(a.block_nnz))
    if a.scale == "mean":
        offset = None
    elif a.offset_col == "nCount_RNA":
        offset = "rowsum"
    else:
        offset = read_obs_derived(a.h5ad, [a.offset_col], False)[a.offset_col].to_numpy(dtype=np.float64)
    log(f"pseudobulking layer '{a.layer}' (nnz={nnz:,}, scale {a.scale}) in {len(blocks)} blocks, {a.workers} workers")
    sums = np.zeros((n_donors, n_var))
    jobs = [(a.h5ad, a.layer, r0, r1, codes[r0:r1], n_donors, n_var,
             offset if offset is None or isinstance(offset, str) else offset[r0:r1]) for r0, r1 in blocks]
    with ProcessPoolExecutor(a.workers) as ex:
        for k, s in enumerate(ex.map(block_donor_sums, jobs), 1):
            sums += s
            if k % max(1, len(jobs) // 10) == 0 or k == len(jobs):
                log(f"  {k}/{len(jobs)} blocks")
    pb = sums / n_cells_per_donor[:, None]
    if a.scale == "log-rate":
        # per-gene pseudocount, identical for all donors so it cannot reorder donors within a gene
        # (a donor-dependent one, e.g. 0.5 / donor mean offset, tracks donor library size)
        np.save(out / "pb_rate.npy", pb)
        pb = np.log(pb + np.maximum(1e-2 * pb.mean(axis=0, keepdims=True), 1e-12))
    np.save(out / "pb_mean.npy", pb)
    np.savetxt(out / "donors.txt", donors, fmt="%s")
    np.savetxt(out / "genes.txt", read_var_names(a.h5ad), fmt="%s")
    with open(out / "args.json", "w") as fh:
        json.dump(vars(a), fh, indent=2)
    log(f"wrote {out}/pb_mean.npy ({n_donors} x {n_var})")


if __name__ == "__main__":
    main()
