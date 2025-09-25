#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simulate subject-level case/control data with gene-set perturbations and cell-type downsampling.

Hierarchy:
  GeneSet g = 1..G
    └── Subject split s (= g): 50/50 (by some subject_id) with fixed RNG per g
        ├── control.h5ad  (unperturbed controls for that split)
        └── For each (c in C) × (keep_frac in D) × (rep in 1..N)
            ├── perturb n-perturb genes with fixed Δ (evenly spaced in [-2,2], excluding 0), NB resample
            ├── library balancing of non-perturbed (eligible) genes per cell (∝ counts)
            ├── downsample target cell type c by keep_frac (others unchanged)
            ├── save case h5ad
            └── save CSV of {gene_id, log2fc, factor}

Filenames & layout (dirname and slogan extracted from --path):
  {dirname}/simulations/{slogan}/geneset_{g:03d}/
    genes.csv
    split_{g:03d}/
        {slogan}_split_{g:03d}__control.h5ad
        {slogan}_split_{g:03d}__case__downsample_{c}_{keep_frac}__rep_{rep}.h5ad

Notes:
- If `obs[subject_id]` absent, throw an error
- Eligible genes are computed on the full dataset BEFORE splitting:
    total_counts > min_total_counts_gene_thresh AND pct_cells_expressed > min_pct_expressed_gene_thresh.
  PERTURBED genes must be drawn from this eligible set, BUT NO GENES ARE REMOVED FROM ANY AnnData.
- Perturbations apply to ALL CELLS (not just a specific cell type). Downsampling happens after perturbation.
- If eligible genes < n_perturb → error and stop.
- If a selected cell type `c` has 0 cells in the CASE cohort BEFORE downsampling → error.
- keep_frac can be 0..1; sampling is without replacement; d=0 allowed (may drop type c entirely).
"""

from __future__ import annotations
import os
# Limit BLAS threads so multiprocessing can fan out
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import math
from pathlib import Path
from typing import List, Tuple, Dict, Any
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from multiprocessing import Pool, cpu_count
import zlib


# --------------------------
# Utilities
# --------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate subject-level case/control H5ADs with gene perturbations and downsampling."
    )
    p.add_argument("--path", required=True, help="Path to input anndata .h5ad (raw counts in X).")
    p.add_argument("--n-perturbed-genes", type=int, required=True, default=400, 
                    help="Number of genes to perturb per gene set.")
    p.add_argument("--cell-type-col", required=True, default="cell_type",
                help="obs column name containing cell types.")
    p.add_argument("--subject-id-col", required=True, default="donor",
                   help="obs column name containing subject IDs (e.g. 'subject.subjectGuid', 'donor', etc.).")
    p.add_argument("--num-gene-sets", type=int, default=1, help="Number of perturbed gene sets (one split per gene set).")
    p.add_argument("--cell-types", nargs="+", required=True, help="List of target cell types C.")
    p.add_argument("--keep-fracs", type=float, nargs="+", required=True,
                   help="List of keep fractions D (0..1). Example: 1.0 0.5 0.25 0.0")
    p.add_argument("--num-reps", type=int, default=1, help="Number of replicates per (c, keep_frac).")
    p.add_argument("--min-total-counts-gene-thresh", type=int, default=200,
                   help="Minimum total counts across all cells to consider a gene ELIGIBLE (no matrix filtering).")
    p.add_argument("--min-pct-expressed-gene-thresh", type=float, default=0.01,
                   help="Minimum fraction of cells with count>0 to consider a gene ELIGIBLE (no matrix filtering).")
    p.add_argument("--seed", type=int, default=12345, help="Base RNG seed.")
    p.add_argument("--processes", type=int, default=None,
                   help="Number of parallel worker processes for case generation (default: min(available, 40)).")
    p.add_argument("--clobber", action="store_true",
                   help="Overwrite existing input h5ad if subject ID must be created or outputs already exist.")
    p.add_argument("--subject-theta", action="store_true",
                   help="Estimate per-subject dispersion (θ) with shrinkage toward global; "
                        "otherwise use global θ for all subjects.")
    p.add_argument("--theta-shrinkage-k", type=float, default=50.0,
                   help="Shrinkage strength toward global dispersion when --subject-theta is set "
                        "(larger K => more shrinkage).")
    
    # p.add_argument("--subject-lfc-sd", type=float, default=0,
    #                help="Std dev of subject-level log2FC around the global Δ for each perturbed gene. Default 0 (no variation).")

    return p.parse_args()


def basename_noext(path: str) -> str:
    b = os.path.basename(path)
    if b.endswith(".h5ad"):
        b = b[:-5]
    return b


def to_csr(X):
    if sparse.issparse(X):
        return X.tocsr()
    else:
        return sparse.csr_matrix(X)


def to_csc(X):
    if sparse.issparse(X):
        return X.tocsc()
    else:
        return sparse.csc_matrix(X)


def compute_eligible_gene_mask(adata: ad.AnnData, min_total: int, min_pct: float) -> np.ndarray:
    """
    Return a boolean mask (length n_vars) of genes that pass the thresholds,
    computed on raw counts in X across all cells. Does NOT modify adata.
    """
    X = to_csc(adata.raw.X)
    n_cells = adata.n_obs

    totals = np.array(X.sum(axis=0)).ravel()
    expressed = np.array((X > 0).sum(axis=0)).ravel()
    pct = expressed / float(n_cells)

    keep = (totals > min_total) & (pct > min_pct)
    return keep


def split_subjects(subject_ids: pd.Series, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    50/50 split at the subject level, deterministic via `seed`.
    Control gets the floor if odd.
    Returns: (control_subjects, case_subjects)
    """
    uniq = np.asarray(pd.Index(subject_ids.unique()).sort_values())
    rng = np.random.default_rng(seed)
    rng.shuffle(uniq)
    n = len(uniq)
    n_ctrl = n // 2
    ctrl = uniq[:n_ctrl]
    case = uniq[n_ctrl:]
    return ctrl, case


def estimate_nb_dispersion_gene(X_csc: sparse.csc_matrix) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimate per-gene mean (μ), variance (σ²), and dispersion θ for a negative binomial model.

    Parameterization:
      Var(X) = μ + μ² / θ

    - If θ → ∞, the variance → μ (i.e., Poisson).
    - Smaller θ → more overdispersion.

    Arguments
    ---------
    X_csc : sparse.csc_matrix
        Gene expression count matrix (cells × genes) in CSC format, raw counts.

    Returns
    -------
    mu : np.ndarray
        Per-gene mean across cells.
    var : np.ndarray
        Per-gene variance across cells.
    theta : np.ndarray
        Per-gene dispersion parameter θ (>= 1e-8).
    """

    # --- Mean per gene
    mu = np.array(X_csc.mean(axis=0)).ravel()

    # --- Compute second moment E[X²]
    X_sq = X_csc.copy()
    X_sq.data = X_sq.data ** 2
    ex2 = np.array(X_sq.mean(axis=0)).ravel()

    # --- Variance = E[X²] - (E[X])²
    var = ex2 - mu**2

    # --- Dispersion calculation
    # From NB variance formula: Var = μ + μ²/θ
    # Rearranging: θ = μ² / (Var - μ)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = (var - mu)
        theta_calc = np.where(denom > 0, (mu**2) / denom, np.inf)

    # --- Replace invalid/infinite values
    # - If variance <= mean, denom <= 0 → implies ~Poisson behavior
    #   → set θ very large (1e12) so NB ~ Poisson.
    # - Ensure no -inf, nan, or zero.
    theta = np.where(np.isfinite(theta_calc), theta_calc, 1e12)
    theta = np.maximum(theta, 1e-8)  # guard against zero/negative

    return mu, var, theta



def nb_sample_mean_disp(mean: float, theta: float, rng: np.random.Generator) -> int:
    """
    Sample a single NB draw using mean/dispersion parameterization via Gamma–Poisson.
    mean: target mean
    theta: dispersion parameter (>0)
    rng: np.random.Generator instance
    """
    if mean <= 0:
        return 0
    lam = rng.gamma(shape=theta, scale=(mean / theta))
    return rng.poisson(lam)


def make_output_paths(base_path: str, geneset_idx: int) -> Tuple[Path, str, str]:
    """
    Build output directory for a given gene set and return (split_dir, dirname, slogan).
    """
    dirname = os.path.dirname(base_path)
    slogan = basename_noext(base_path)
    root = Path(dirname) / "simulations" / f"{slogan}" / f"geneset_{geneset_idx:03d}" / f"split_{geneset_idx:03d}"
    return root, dirname, slogan


def add_case_control_flag(adata: ad.AnnData, label: str):
    adata.obs["case_control"] = pd.Categorical([label] * adata.n_obs, categories=["control", "case"])


def build_delta_vector(num_genes: int = 400) -> np.ndarray:
    """
    Build evenly spaced Δ in [-2, 2] excluding 0, then shuffle.
    """
    grid_neg = np.linspace(-2, -1e-12, num_genes // 2)
    grid_pos = np.linspace(1e-12, 2, num_genes - len(grid_neg))
    vec = np.concatenate([grid_neg, grid_pos])
    rng = np.random.default_rng(0)
    rng.shuffle(vec)
    return vec


def derive_rng(seed_base: int, geneset: int, c: str, keep_frac: float, rep: int) -> np.random.Generator:
    key = f"{seed_base}|{geneset}|{c}|{keep_frac:.6f}|{rep}".encode()
    h = zlib.adler32(key) & 0xffffffff  # stable across sessions
    return np.random.default_rng(h)

def balance_library_counts_per_cell(
    counts_row: np.ndarray,
    perturbed_mask: np.ndarray,
    eligible_mask: np.ndarray,   # === CHANGED: include eligibility
    delta_total: int,
    rng: np.random.Generator,
):
    """
    Adjust non-perturbed (and eligible) genes so total change equals -delta_total.
    Weights selection proportional to current counts; non-negative constraint enforced.
    """
    if delta_total == 0:
        return

    nonpert_mask = (~perturbed_mask) & eligible_mask  # === CHANGED
    counts_nonpert = counts_row[nonpert_mask]
    if counts_nonpert.size == 0:
        return

    idx_nonpert = np.flatnonzero(nonpert_mask)

    if delta_total > 0:
        available = counts_nonpert.sum()
        to_remove = min(delta_total, int(available))
        if to_remove <= 0:
            return
        cvec = counts_row[idx_nonpert].astype(np.int64, copy=True)
        while to_remove > 0 and cvec.sum() > 0:
            probs = cvec / cvec.sum()
            j = rng.choice(len(cvec), p=probs)
            if cvec[j] > 0:
                cvec[j] -= 1
                to_remove -= 1
        counts_row[idx_nonpert] = cvec
    else:
        to_add = -delta_total
        cvec = counts_row[idx_nonpert].astype(np.int64, copy=True)
        for _ in range(to_add):
            total = cvec.sum()
            if total == 0:
                j = rng.integers(0, len(cvec))
            else:
                probs = cvec / total
                j = rng.choice(len(cvec), p=probs)
            cvec[j] += 1
        counts_row[idx_nonpert] = cvec


def generate_control_h5ad(
    adata_full: ad.AnnData,
    subject_ctrl: np.ndarray,
    split_dir: Path,
    subject_id_col: str,
    slogan: str,
    clobber: bool = False,
):
    ctrl_mask = adata_full.obs[subject_id_col].isin(subject_ctrl).values
    ad_ctrl = adata_full[ctrl_mask].copy()
    add_case_control_flag(ad_ctrl, "control")
    ad_ctrl.uns["simulation"] = {
        "split": int(split_dir.name.split("_")[-1]),
        "case_control": "control",
        "cell_type": None,
        "keep_frac": None,
        "replicate": None,
        "perturbed_genes": {},
    }
    # if there is an adata.raw, replace adata.X with adata.raw.X
    if ad_ctrl.raw is not None:
        X = ad_ctrl.raw.X
        if sparse.issparse(X):
            X = X.astype(np.float64)
        else:
            X = X.astype(np.float64, copy=False)
        ad_ctrl.X = X          # replace main matrix
        ad_ctrl.var = ad_ctrl.raw.var.copy()   # make sure genes match raw
    else: 
        print("Warning: adata.raw is None; using adata.X as-is for control.")

    out = split_dir / f"{slogan}_split_{int(split_dir.name.split('_')[-1]):03d}__control.h5ad"
    if out.exists():
        # check that the control subjects are identical
        ad_existing = ad.read_h5ad(out, backed="r")
        existing_subjects = set(ad_existing.obs[subject_id_col].unique())
        new_subjects = set(ad_ctrl.obs[subject_id_col].unique())
        if existing_subjects != new_subjects:
                raise RuntimeError(
                    f"Control file already exists but subjects differ: {out}\n"
                    f"Existing subjects ({len(existing_subjects)}): {sorted(existing_subjects)}\n"
                    f"New subjects ({len(new_subjects)}): {sorted(new_subjects)}\n"
                    "Use --clobber to overwrite."
                )
    if out.exists() and not clobber:
        print(f"Control file already exists: {out}")
    else:
        print(f"Writing control file: {out}")
        # make sure to create directories if needed
        out.parent.mkdir(parents=True, exist_ok=True)
        write_sanitized_h5ad(ad_ctrl, path=out, clobber=clobber, keep_raw=False)
    return out

def _sanitize_df(df: pd.DataFrame, where: str) -> None:
    # AnnData reserves the column name "_index" in obs/var/raw.var
    if "_index" in df.columns:
        df.rename(columns={"_index": f"{where}_index_orig"}, inplace=True)
        df.set_index(f"{where}_index_orig", inplace=True, drop=False)

def write_sanitized_h5ad(
    adata: ad.AnnData,
    path: Path,
    clobber: bool = False,
    keep_raw: bool = False,
) -> None:
    """
    Write H5AD after removing reserved column names.
    If keep_raw=False (default), drops .raw entirely (simplest & safest).
    If keep_raw=True, rebuilds .raw from sanitized copies of X/var/obs.
    """
    path = Path(path)
    if path.exists() and not clobber:
        raise FileExistsError(f"File already exists: {path}. Use --clobber to overwrite.")

    path.parent.mkdir(parents=True, exist_ok=True)

    # Sanitize main obs/var (this also affects any future raw we rebuild)
    _sanitize_df(adata.obs, "obs")
    _sanitize_df(adata.var, "var")

    if not keep_raw:
        # simplest: drop raw entirely
        adata.raw = None
    else:
        # rebuild a clean raw from sanitized copies (no private attribute poking)
        raw_obs = adata.obs.copy()
        raw_var = adata.var.copy()
        _sanitize_df(raw_obs, "raw_obs")  # defensive; shouldn't have "_index"
        _sanitize_df(raw_var, "raw_var")
        adata.raw = ad.AnnData(X=adata.X.copy(), var=raw_var, obs=raw_obs)

    adata.write(path, compression=None)

def _case_worker(
    args: Tuple[
        int,            # geneset_idx
        Path,           # split_dir
        str,            # slogan
        ad.AnnData,     # ad_case_base
        np.ndarray,     # gene_indices
        np.ndarray,     # delta_vec
        np.ndarray,     # eligible_mask_global
        np.ndarray,     # theta_global
        np.ndarray,     # subj_codes
        np.ndarray,     # mu_subject
        np.ndarray,     # subject_deltas
        Optional[np.ndarray], # theta_subject
        str,            # cell_type_col
        str,            # c_value
        float,          # keep_frac
        int,            # rep
        int             # seed_base
    ]
) -> Tuple[str, str]:
    """
    Simulate one CASE dataset for a given gene set, subject split, and (cell type, keep_frac, replicate).

    Parameters
    ----------
    args : tuple
        A packed tuple containing:

        geneset_idx : int
            Index of the gene set (1..G).
        split_dir : Path
            Output directory for this split.
        slogan : str
            Basename of the input file, used in output filenames.
        ad_case_base : AnnData
            AnnData slice with cells from CASE subjects only, before perturbation.
        gene_indices : np.ndarray, shape (n_perturb,)
            Column indices of perturbed genes (e.g. length 400).
        delta_vec : np.ndarray, shape (n_perturb,)
            Global log2 fold changes assigned to perturbed genes.
        eligible_mask_global : np.ndarray, shape (n_genes,)
            Boolean mask over all genes, True for genes that passed eligibility thresholds.
        theta_global : np.ndarray, shape (n_genes,)
            Global NB dispersion estimates for each gene (used if `theta_subject` is None).
        subj_codes : np.ndarray, shape (n_cells,)
            Integer subject code (0..n_subj-1) for each cell in `ad_case_base`.
        mu_subject : np.ndarray, shape (n_subj, n_genes)
            Subject-level baseline means. Entry [s, g] = average count of gene g across
            all cells belonging to subject s.
        subject_deltas : np.ndarray, shape (n_subj, n_perturb)
            Subject-specific log2 fold changes. Entry [s, j] = Δ for subject s and perturbed gene j.
        theta_subject : np.ndarray, shape (n_subj, n_genes), or None
            Subject-specific dispersion estimates. If None, fall back to `theta_global`.
        cell_type_col : str
            obs column name containing cell types.
        c_value : str
            Target cell type to downsample.
        keep_frac : float
            Fraction of cells of type `c_value` to keep (0..1).
        rep : int
            Replicate index (1..N).
        seed_base : int
            Base RNG seed, combined with job identifiers to derive deterministic RNG.

    Returns
    -------
    case_h5ad_path : str
        Path to the saved perturbed case AnnData file.
    genes_csv_path : str
        Path to the `genes.csv` file containing the global perturbed gene set.
    
    Notes
    -----
    - Perturbations are applied to all cells, with subject-specific means and Δ.
    - Library size is balanced per cell across non-perturbed eligible genes.
    - Downsampling is applied after perturbation.
    - The returned `ad_case.uns['simulation']` metadata includes split index,
      case/control label, cell type, keep_frac, replicate, and perturbed genes.
    """
    (
        geneset_idx,                # int
        split_dir,                  # Path
        slogan,                     # str
        ad_case_base,               # ad.AnnData               
        gene_indices,               # np.ndarray (shape (n_perturb,))   
        delta_vec,                  # np.ndarray (shape (n_perturb,))                    
        eligible_mask_global,       # np.ndarray (shape (n_genes,))  
        theta_global,               # np.ndarray (shape (n_genes,))
        subj_codes,                 # np.ndarray (shape (n_cells,))
        mu_subject,                 # np.ndarray (shape (n_subj, n_genes))  
        theta_subject,              # np.ndarray (shape (n_subj, n_genes)) or None  
        cell_type_col,              # str
        c_value,                    # str
        keep_frac,                  # float
        rep,                        # int
        seed_base                   # int
    ) = args

    rng = derive_rng(seed_base, geneset_idx, c_value, keep_frac, rep)

    # Views & shapes
    Xcsr = ad_case_base.X.tocsr()
    Xcsc = ad_case_base.X.tocsc()   # fast column access for original columns
    n_cells, n_genes = Xcsr.shape
    Gp = len(gene_indices)

    # 1) Build NEW perturbed columns (dense) in a vectorized way
    #    new_cols shape: (n_cells, Gp)
    new_cols = np.empty((n_cells, Gp), dtype=np.int32)

    # Cells grouped by subject
    subj_codes = subj_codes.astype(np.int32, copy=False)
    n_subj = mu_subject.shape[0]
    # Precompute which cells belong to each subject once
    cells_by_subj = [np.where(subj_codes == s)[0] for s in range(n_subj)]

    for j_local, gidx in enumerate(gene_indices):
        # Global Δ for this gene
        d = float(delta_vec[j_local])
        # Subject means for this gene (length n_subj)
        mu_s = mu_subject[:, gidx].astype(float, copy=False)
        # Subject-level perturbed mean for this gene
        mu_s_prime = (2.0 ** d) * mu_s

        # Dispersion: subject-specific if provided, else global
        if theta_subject is not None:
            theta_s = theta_subject[:, gidx].astype(float, copy=False)
        else:
            theta_s = np.full(n_subj, float(theta_global[gidx]), dtype=float)

        # Fill this gene’s column for all cells, grouped by subject (vectorized per subject)
        col_vec = np.empty(n_cells, dtype=np.int32)
        for s in range(n_subj):
            idx_cells = cells_by_subj[s]
            if idx_cells.size == 0:
                continue
            mu_prime = float(mu_s_prime[s])
            th = float(theta_s[s])
            if mu_prime <= 0:
                col_vec[idx_cells] = 0
            else:
                # Vectorized Gamma–Poisson for this subject block
                lam = rng.gamma(shape=th, scale=(mu_prime / th), size=idx_cells.size)
                col_vec[idx_cells] = rng.poisson(lam).astype(np.int32, copy=False)
        new_cols[:, j_local] = col_vec

    # Old perturbed columns (dense) to compute per-row delta quickly
    old_cols = Xcsc[:, gene_indices].toarray().astype(np.int32, copy=False)

    # 2) Per-row total delta across perturbed genes (k' - k)
    delta_row = (new_cols - old_cols).sum(axis=1, dtype=np.int64)  # shape (n_cells,)

    # 3) Switch to LIL for cheap row writes
    X_lil = Xcsr.tolil(copy=True)
    del Xcsr  # free CSR

    eligible_mask = np.asarray(eligible_mask_global, dtype=bool)
    pert_mask = np.zeros(n_genes, dtype=bool)
    pert_mask[gene_indices] = True
    nonpert_eligible_mask = (~pert_mask) & eligible_mask

    # Pre-extract non-perturbed-eligible column indices once
    idx_nonpert_eligible = np.flatnonzero(nonpert_eligible_mask)

    # Also, create a fast lookup for each row's current non-pert-eligible counts as dense
    # We'll get row nonzeros from LIL internals and slice.
    for i in range(n_cells):
        # === Build dense row only for columns we need to touch ===
        # Start with all zeros
        # 3a) Set perturbed gene counts from new_cols
        # 3b) Read current non-pert-eligible counts to dense for balancing
        # 3c) Write back both segments

        # (a) start from previous row nnz to get existing values
        cols_i = X_lil.rows[i]
        data_i = X_lil.data[i]
        # Build a tiny dict for the small set we will read (faster than full dense)
        row_map = {c: int(v) for c, v in zip(cols_i, data_i)}

        # Current counts for non-pert eligible columns as dense vector
        c = np.array([row_map.get(col, 0) for col in idx_nonpert_eligible], dtype=np.int64)

        # Set the perturbed columns new counts and compute delta already done
        # We'll write perturbed columns after balancing.

        # Balance per-row total change over non-pert eligible genes via multinomial
        T = int(delta_row[i])
        if T != 0 and idx_nonpert_eligible.size > 0:
            tot = c.sum()
            if T > 0:
                # Need to REMOVE T from non-pert eligible
                # If tot == 0, nothing to remove; else allocate removals with capping
                if tot > 0:
                    remaining = T
                    # Cap-aware removal: do a few passes to respect c_i bounds
                    # (Typically 1–2 passes suffice unless T is huge relative to tot)
                    for _ in range(4):
                        if remaining <= 0 or c.sum() == 0:
                            break
                        p = c / c.sum()
                        draw = rng.multinomial(remaining, p)
                        # Clip by available counts
                        remove = np.minimum(draw, c)
                        c -= remove
                        remaining -= int(remove.sum())
            else:
                # Need to ADD -T over non-pert eligible; allocate by proportions (include +1 to avoid dead-ends)
                add_T = -T
                w = c
                if w.sum() == 0:
                    # uniform if all zero
                    p = np.ones_like(w, dtype=float) / float(w.size)
                else:
                    p = w / w.sum()
                add = rng.multinomial(add_T, p)
                c += add

        # Now write row back into LIL:
        # 1) start with perturbed columns
        nz_cols = []
        nz_vals = []

        # perturbed
        for j_local, gidx in enumerate(gene_indices):
            v = int(new_cols[i, j_local])
            if v != 0:
                nz_cols.append(int(gidx))
                nz_vals.append(v)

        # eligible non-pert (balanced)
        for off, col in enumerate(idx_nonpert_eligible):
            v = int(c[off])
            if v != 0:
                nz_cols.append(int(col))
                nz_vals.append(v)

        # ineligible genes: keep previous counts (we never touched them)
        for c0, v0 in zip(cols_i, data_i):
            if (c0 in pert_mask.nonzero()[0]) or (nonpert_eligible_mask[c0]):
                # skip, already set
                continue
            if v0 != 0:
                nz_cols.append(int(c0))
                nz_vals.append(int(v0))

        # Assign the new sparse row structure
        X_lil.rows[i] = nz_cols
        X_lil.data[i] = nz_vals

    # Convert once to CSR
    X = X_lil.tocsr()
    del X_lil

    # Build AnnData
    ad_case = ad.AnnData(X=X, obs=ad_case_base.obs.copy(), var=ad_case_base.var.copy(), uns=ad_case_base.uns.copy())

    # Downsample cells of type c_value by keep_frac
    mask_c = (ad_case.obs[cell_type_col].astype(str).values == str(c_value))
    n_c = int(mask_c.sum())
    if n_c == 0:
        raise RuntimeError(
            f"No cells with {cell_type_col} == '{c_value}' in CASE cohort BEFORE downsampling (geneset {geneset_idx})."
        )
    if keep_frac < 1.0:
        n_keep = int(round(keep_frac * n_c))
        idx_c = np.flatnonzero(mask_c)
        keep_idx_c = set(rng.choice(idx_c, size=n_keep, replace=False).tolist()) if n_keep > 0 else set()
        final_mask = np.ones(ad_case.n_obs, dtype=bool)
        for idx in idx_c:
            if idx not in keep_idx_c:
                final_mask[idx] = False
        ad_case = ad_case[final_mask].copy()

    # Flags + metadata (global Δ only, as requested)
    add_case_control_flag(ad_case, "case")
    pert_genes = ad_case.var_names[gene_indices]
    gene_map = {g: {"log2fc": float(delta_vec[i]), "factor": float(2.0 ** delta_vec[i])}
                for i, g in enumerate(pert_genes)}
    ad_case.uns["simulation"] = {
        "split": geneset_idx,
        "case_control": "case",
        "cell_type": c_value,
        "keep_frac": float(keep_frac),
        "replicate": int(rep),
        "perturbed_genes": gene_map,
    }

    # Save
    case_h5ad = split_dir / f"{slogan}_split_{geneset_idx:03d}__case__downsample_{c_value}_{keep_frac}__rep_{rep}.h5ad"
    case_h5ad = case_h5ad.with_name(case_h5ad.name.replace(" ", "_"))
    # sanitize before write if needed
    write_sanitized_h5ad(ad_case, case_h5ad, clobber=True, keep_raw=False)
    ad_case.write(case_h5ad, compression=None)

    genes_csv = split_dir.parent / "genes.csv"
    return str(case_h5ad), str(genes_csv)


def main():
    args = parse_args()

    base_path = args.path
    dirname = os.path.dirname(base_path)
    slogan = basename_noext(base_path)

    # Load and ensure subject id exists (your ensure_subject_id can still auto-create if missing)
    print("Loading data...")
    adata = ad.read_h5ad(base_path, backed=None)  # in-memory

    # Validate required columns
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"--cell-type-col '{args.cell_type_col}' not found in obs.")
    if args.subject_id_col not in adata.obs.columns:
        raise ValueError(f"--subject-id-col '{args.subject_id_col}' not found in obs.")
    if args.keep_fracs is None or len(args.keep_fracs) == 0:
        raise ValueError("At least one --keep-fracs value required.")
    if any ((d < 0.0) or (d > 1.0) for d in args.keep_fracs):
        raise ValueError("--keep-fracs values must be in [0, 1].")

    # Fail fast if any requested cell type is missing
    missing = set(args.cell_types) - set(adata.obs[args.cell_type_col].unique())
    if missing:
        raise ValueError(f"Cell types not found in obs['{args.cell_type_col}']: {missing}")

    # Eligibility mask over ALL genes (no filtering of matrix)
    eligible_mask_global = compute_eligible_gene_mask(
        adata,
        min_total=args.min_total_counts_gene_thresh,
        min_pct=args.min_pct_expressed_gene_thresh,
    )
    n_eligible = int(eligible_mask_global.sum())
    if n_eligible < args.n_perturbed_genes:
        raise ValueError(f"Eligible genes ({n_eligible}) < {args.n_perturbed_genes}. Increase data or relax thresholds.")

    # Coerce key obs cols to str
    adata.obs[args.cell_type_col] = adata.obs[args.cell_type_col].astype(str)
    adata.obs[args.subject_id_col] = adata.obs[args.subject_id_col].astype(str)

    # Processes
    procs = min(cpu_count(), 40) if args.processes is None else max(1, args.processes)

    # Δ template for perturbed genes
    delta_vec_template = build_delta_vector(args.n_perturbed_genes)

    # For each gene set (one split per gene set)
    for g in range(1, args.num_gene_sets + 1):
        split_dir, _, _ = make_output_paths(base_path, g)
        split_dir.mkdir(parents=True, exist_ok=True)

        # Subject split for this geneset (seeded)
        ctrl_subjects, case_subjects = split_subjects(adata.obs[args.subject_id_col], seed=args.seed + g)

        # Save control once per split (your helper can verify same subjects unless --clobber)
        ctrl_h5 = generate_control_h5ad(
            adata_full=adata,
            subject_ctrl=ctrl_subjects,
            split_dir=split_dir,
            subject_id_col = args.subject_id_col,
            slogan = slogan,
            clobber = args.clobber)
        print(f"[geneset {g:03d}] control saved: {ctrl_h5}")

        # CASE base slice (before perturbation)
        case_mask = adata.obs[args.subject_id_col].isin(case_subjects).values
        ad_case_base = adata[case_mask].copy()

        # Pick args.n_perturbed_genes perturbed genes ONLY from eligible set (global indexing)
        rng_gene = np.random.default_rng(args.seed + 1000 * g)
        eligible_indices = np.flatnonzero(eligible_mask_global)
        gene_indices = rng_gene.choice(eligible_indices, size=args.n_perturbed_genes, replace=False)

        # Global Δ for this geneset (shuffled copy of template)
        delta_vec = delta_vec_template.copy()
        rng_gene.shuffle(delta_vec)

        # Validate that each requested cell type exists in CASE base before downsampling
        cell_type_col = args.cell_type_col
        case_counts_by_c = ad_case_base.obs[cell_type_col].value_counts().to_dict()
        for c in args.cell_types:
            if case_counts_by_c.get(str(c), 0) == 0:
                raise RuntimeError(
                    f"[geneset {g:03d}] No cells with {cell_type_col} == '{c}' in CASE cohort BEFORE downsampling."
                )

        # === Subject structure and parameters (CASE cohort only) ===
        subjects_case = ad_case_base.obs[args.subject_id_col].astype(str).values
        subj_codes, subj_uniqs = pd.factorize(subjects_case, sort=True)
        n_subj = len(subj_uniqs)

        # Compute global NB params ONCE for case cohort
        print(f"[geneset {g:03d}] Estimating NB dispersion parameters...")
        X_case_csc = ad_case_base.X.tocsc()
        mu_global, var_global, theta_global = estimate_nb_dispersion_gene(X_case_csc)

        # Per-subject per-gene baseline means μ_{g,s}
        print(f"[geneset {g:03d}] Computing subject-level baseline means...")
        mu_subject = np.zeros((n_subj, ad_case_base.n_vars), dtype=float)
        for s_idx in range(n_subj):
            rows = np.where(subj_codes == s_idx)[0]
            if rows.size:
                mu_subject[s_idx, :] = np.asarray(X_case_csc[rows, :].mean(axis=0)).ravel()

        # Optional: per-subject dispersion with shrinkage toward global
        theta_subject = None
        if getattr(args, "subject_theta", False):
            print(f"[geneset {g:03d}] Estimating subject-level dispersion with shrinkage...")
            
            # Initialize an array of shape (n_subj, n_genes) to hold subject-specific θ
            theta_subject = np.empty_like(mu_subject)

            # Global precision (φ = 1/θ) used as shrinkage target
            phi_global = 1.0 / np.maximum(theta_global, 1e-8)

            K = float(getattr(args, "theta_shrinkage_k", 50.0)) 
            for s_idx in range(n_subj):
                rows = np.where(subj_codes == s_idx)[0]  # cells belonging to subject s
                if rows.size:
                    # Compute per-gene mean and variance for this subject
                    Xs = X_case_csc[rows, :]
                    mu_s = np.asarray(Xs.mean(axis=0)).ravel()
                    Xs_sq = Xs.copy()
                    Xs_sq.data = Xs_sq.data ** 2
                    ex2_s = np.asarray(Xs_sq.mean(axis=0)).ravel()
                    var_s = ex2_s - mu_s**2

                    # Estimate NB dispersion θ_s = μ^2 / (var - μ), with guards
                    with np.errstate(divide="ignore", invalid="ignore"):
                        denom = (var_s - mu_s)
                        theta_s = np.where(denom > 0, (mu_s**2) / denom, np.inf)
                    theta_s = np.where(np.isfinite(theta_s), theta_s, 1e12)
                    theta_s = np.maximum(theta_s, 1e-8)
                else:
                    # If no cells for this subject, assign very large θ (≈Poisson)
                    theta_s = np.full(ad_case_base.n_vars, 1e12, dtype=float)

                # Compute shrinkage weight w based on sample size
                n_s = max(1, rows.size)
                w = n_s / (n_s + K)  # fraction of weight given to subject estimate

                # Work on precision scale (φ = 1/θ), since precisions combine linearly
                phi_s = 1.0 / np.maximum(theta_s, 1e-8)

                # Shrink subject-specific precision toward the global precision
                phi_shrunk = w * phi_s + (1.0 - w) * phi_global

                # Convert back to dispersion θ = 1/φ
                theta_subject[s_idx, :] = 1.0 / np.maximum(phi_shrunk, 1e-8)

        # CODE FOR SUBJECT-SPECIFIC LFCs (currently disabled)        
        # Subject-specific log2FCs for the selected (perturbed) genes:
        # Δ_{s,j} = Δ_j + Normal(0, subject_lfc_sd)
        # print(f"[geneset {g:03d}] Sampling subject-specific log2FCs (sd={args.subject_lfc_sd})...")
        # rng_subj = np.random.default_rng(args.seed + 1000 * g + 7)

        # If subject_lfc_sd=0, all subjects have the same Δ_j; otherwise add Gaussian noise
        # if float(getattr(args, "subject_lfc_sd", 0)) > 0.0:
        #     subject_deltas = (delta_vec[None, :] + rng_subj.normal(
        #         loc=0.0, scale=float(getattr(args, "subject_lfc_sd", 0.25)), size=(n_subj, len(delta_vec))
        #     ))
        # else:
        #     subject_deltas = np.tile(delta_vec[None, :], (n_subj, 1))

        # Free the CSC to lower peak memory before multiprocessing
        del X_case_csc

        # Write genes.csv once per geneset (global Δ/factors)
        geneset_dir = split_dir.parent  # .../geneset_{g:03d}
        genes_csv = geneset_dir / "genes.csv"
        if (not genes_csv.exists()) or args.clobber:
            df_csv = pd.DataFrame({
                "gene_id": ad_case_base.var_names[gene_indices].values,
                "log2fc": [float(x) for x in delta_vec],
                "factor": [float(2.0 ** x) for x in delta_vec],
            })
            df_csv.to_csv(genes_csv, index=False)

        # Prepare jobs for all (c, keep_frac, rep)
        jobs = []
        for c in args.cell_types:
            for d in args.keep_fracs:
                for rep in range(1, args.num_reps + 1):  
                    # Skip redundant keep_fracs
                    if (d == 1.0) and (rep > 1):
                        continue

                    jobs.append((
                        g, 
                        split_dir, 
                        slogan,
                        ad_case_base, 
                        gene_indices, 
                        delta_vec, 
                        eligible_mask_global,
                        theta_global,
                        subj_codes, 
                        mu_subject, 
                        theta_subject,  
                        cell_type_col, 
                        c, 
                        float(d), 
                        rep, 
                        args.seed
                    ))

        # Run in parallel (note: passing big arrays; tune 'procs' per RAM)
        if procs == 1:
            results = [_case_worker(job) for job in jobs]
        else:
            with Pool(processes=procs) as pool:
                results = pool.map(_case_worker, jobs)

        for case_h5ad_path, genes_csv_path in results:
            print(f"[geneset {g:03d}] case saved: {case_h5ad_path}")

    print("Done.")


if __name__ == "__main__":
    main()
