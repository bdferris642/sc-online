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
    p.add_argument("--n-perturb", type=int, required=True, default=400, help="Number of genes to perturb per gene set.")
    p.add_argument("--subject-id-col", required=False, help="obs column name containing subject IDs")
    p.add_argument("--cell-type-col", required=True, help="obs column name containing cell types.")
    p.add_argument("--num-gene-sets", type=int, default=1, help="Number of perturbed gene sets (one split per gene set).")
    p.add_argument("--cell-types", nargs="+", required=True, help="List of target cell types C.")
    p.add_argument("--keep-fracs", type=float, nargs="+", required=True,
                   help="List of keep fractions D (0..1). Example: 1.0 0.5 0.25 0.0")
    p.add_argument("--num-reps", type=int, default=1, help="Number of replicates per (c, keep_frac).")
    p.add_argument("--min-total-counts-gene-thresh", type=int, default=200,
                   help="Minimum total counts across all cells to consider a gene ELIGIBLE (no filtering of matrix).")
    p.add_argument("--min-pct-expressed-gene-thresh", type=float, default=0.01,
                   help="Minimum fraction of cells with count>0 to consider a gene ELIGIBLE (no filtering of matrix).")
    p.add_argument("--seed", type=int, default=12345, help="Base RNG seed.")
    p.add_argument("--processes", type=int, default=None,
                   help="Number of parallel worker processes for case generation (default: min(available, 40)).")
    p.add_argument("--clobber", action="store_true", 
                   help="Overwrite existing input h5ad.")
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


# === CHANGED: compute eligibility mask only; do NOT filter adata
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
        ad_ctrl.write(out, compression=None)
    return out


def _case_worker(args):
    (
        geneset_idx, split_dir, slogan,
        ad_case_base,
        gene_indices, delta_vec, eligible_mask_global,
        mu_global, theta_global,
        subj_codes, mu_subject, subject_deltas, theta_subject,
        cell_type_col, c_value, keep_frac, rep, seed_base
    ) = args

    rng = derive_rng(seed_base, geneset_idx, c_value, keep_frac, rep)

    X = to_csr(ad_case_base.X)
    n_cells, n_genes = X.shape

    pert_mask = np.zeros(n_genes, dtype=bool)
    pert_mask[gene_indices] = True
    eligible_mask = np.asarray(eligible_mask_global, dtype=bool)

    for i in range(n_cells):
        row = X.getrow(i).toarray().ravel().astype(np.int64, copy=True)
        s_idx = int(subj_codes[i])
        delta_sum = 0

        for idx_local, gidx in enumerate(gene_indices):
            k = int(row[gidx])

            # subject-level baseline mean and subject-level delta for this gene
            mu_gs     = float(mu_subject[s_idx, gidx])
            delta_sj  = float(subject_deltas[s_idx, idx_local])

            # choose dispersion: per-subject if available, otherwise global
            th_j = float(theta_subject[s_idx, gidx]) if (theta_subject is not None) else float(theta_global[gidx])

            # subject-specific perturbed mean (no dependence on k to avoid purity)
            mu_prime_cell = (2.0 ** delta_sj) * mu_gs

            k_prime = nb_sample_mean_disp(mu_prime_cell, th_j, rng)
            k_prime = max(0, int(k_prime))
            row[gidx] = k_prime
            delta_sum += (k_prime - k)

        # Balance on non-perturbed ∩ eligible genes (proportional to counts)
        balance_library_counts_per_cell(row, pert_mask, eligible_mask, delta_sum, rng)

        # Clamp and write back
        np.maximum(row, 0, out=row)
        X.data[X.indptr[i]:X.indptr[i+1]] = 0
        nz = np.flatnonzero(row)
        X[i, nz] = row[nz]

    # Replace matrix
    ad_case = ad.AnnData(X=X, obs=ad_case_base.obs.copy(), var=ad_case_base.var.copy(), uns=ad_case_base.uns.copy())

    # Downsample cells of type c_value by keep_frac
    mask_c = (ad_case.obs[cell_type_col].astype(str).values == str(c_value))
    n_c = int(mask_c.sum())
    if n_c == 0:
        raise RuntimeError(
            f"No cells with cell_type_col == '{c_value}' in CASE cohort BEFORE downsampling (geneset {geneset_idx})."
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

    add_case_control_flag(ad_case, "case")

    # Store metadata
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

    # Save h5ad
    case_h5ad = split_dir / f"{slogan}_split_{geneset_idx:03d}__case__downsample_{c_value}_{keep_frac}__rep_{rep}.h5ad"
    case_h5ad = case_h5ad.with_name(case_h5ad.name.replace(" ", "_"))
    ad_case.write(case_h5ad, compression=None)

    return str(case_h5ad)



def _case_worker_1(args):
    """
    Worker to generate one (c, keep_frac, rep) case dataset for a given gene set/split.
    Designed for parallel execution.
    """
    (
        geneset_idx, split_dir, slogan,
        ad_case_base,            # AnnData slice: case subjects only, BEFORE perturbation
        gene_indices,            # indices of selected genes (global var index, same in this slice)
        delta_vec,               # length n-perturb log2FC vector aligned to gene_indices
        eligible_mask_global,    
        mu, theta,                 # precomputed NB params for the case cohort (on ad_case_base)
        cell_type_col, c_value, keep_frac, rep, seed_base
    ) = args

    rng = derive_rng(seed_base, geneset_idx, c_value, keep_frac, rep)

    # Convert X to CSR for row-wise updates (cells)
    X = to_csr(ad_case_base.raw.X)
    n_cells, n_genes = X.shape

    # Build perturbed mask (length n_genes)
    pert_mask = np.zeros(n_genes, dtype=bool)
    pert_mask[gene_indices] = True

    # Eligible mask (same indexing as var)
    eligible_mask = eligible_mask_global.astype(bool)

    eps = 1e-8

    # Apply perturbations row-wise
    for i in range(n_cells):
        row = X.getrow(i).toarray().ravel().astype(np.int64, copy=True)
        delta_sum = 0

        for idx_local, gidx in enumerate(gene_indices):
            k = int(row[gidx]) # current count
            mu_j = float(mu[gidx]) # current mean
            th_j = float(theta[gidx]) # current dispersion
            if mu_j <= 0: # no expression; leave as is
                k_prime = k
            else: # adding some epsilons, but basically mu_prime_cell ≈ k * (2.0 ** delta)
                p_cell = (k + eps) / (mu_j + eps) # cell's expression relative to gene mean
                mu_prime_gene = (2.0 ** float(delta_vec[idx_local])) * mu_j # perturbed mean
                mu_prime_cell = p_cell * mu_prime_gene # perturbed mean for this cell
                k_prime = nb_sample_mean_disp(
                    mean=mu_prime_cell, 
                    theta=th_j, 
                    rng=rng) # new count
            k_prime = max(0, int(k_prime)) # ensure non-negative integer
            row[gidx] = k_prime
            delta_sum += (k_prime - k)

        # Library balancing over non-perturbed ELIGIBLE genes
        balance_library_counts_per_cell(row, pert_mask, eligible_mask, delta_sum, rng)

        np.maximum(row, 0, out=row)

        # Write back (re-sparsify row)
        X.data[X.indptr[i]:X.indptr[i+1]] = 0
        nz_cols = np.flatnonzero(row)
        X[i, nz_cols] = row[nz_cols]

    # Replace matrix
    ad_case = ad.AnnData(X=X, obs=ad_case_base.obs.copy(), var=ad_case_base.var.copy(), uns=ad_case_base.uns.copy())

    # Downsample cells of type c_value by keep_frac
    mask_c = (ad_case.obs[cell_type_col].astype(str).values == str(c_value))
    n_c = int(mask_c.sum())
    if n_c == 0:
        raise RuntimeError(
            f"No cells with cell_type_col == '{c_value}' in CASE cohort BEFORE downsampling (geneset {geneset_idx})."
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

    add_case_control_flag(ad_case, "case")

    # Store metadata
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

    # Save h5ad
    case_h5ad = split_dir / f"{slogan}_split_{geneset_idx:03d}__case__downsample_{c_value}_{keep_frac}__rep_{rep}.h5ad"
    case_h5ad = case_h5ad.with_name(case_h5ad.name.replace(" ", "_"))
    ad_case.write(case_h5ad, compression=None)

    return str(case_h5ad)


def main():
    args = parse_args()

    base_path = args.path
    dirname = os.path.dirname(base_path)
    slogan = basename_noext(base_path)

    # Load and ensure subject_id
    print("Loading data...")
    adata = ad.read_h5ad(base_path, backed=None)  # in-memory

    print("Checking for subject_id...")
    # if subject_id_col is not present, throw an error
    if args.subject_id_col not in adata.obs.columns:
        raise ValueError(f"--subject-id-col '{args.subject_id_col}' not found in obs.")

    # Validate cell_type_col exists
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"--cell-type-col '{args.cell_type_col}' not found in obs.")
    
    # Fail fast if any requested cell type is missing
    missing = set(args.cell_types) - set(adata.obs[args.cell_type_col].unique())
    if missing:
        raise ValueError(f"Cell types not found in obs['{args.cell_type_col}']: {missing}")

    eligible_mask_global = compute_eligible_gene_mask(
        adata,
        min_total=args.min_total_counts_gene_thresh,
        min_pct=args.min_pct_expressed_gene_thresh,
    )
    n_eligible = int(eligible_mask_global.sum())
    if n_eligible < args.n_perturb:
        raise ValueError(
            f"Eligible genes ({n_eligible}) < {args.n_perturb}. Increase data or relax thresholds."
        )

    # Basic types
    adata.obs[args.cell_type_col] = adata.obs[args.cell_type_col].astype(str)
    adata.obs[args.subject_id_col] = adata.obs[args.subject_id_col].astype(str)

    # Processes
    procs = min(cpu_count(), 40) if args.processes is None else max(1, args.processes)

    # Δ vector for args.n_perturb genes
    delta_vec_template = build_delta_vector(args.n_perturb)

    # For each gene set (one split per gene set)
    for g in range(1, args.num_gene_sets + 1):
        split_dir, _, _ = make_output_paths(base_path, g)

        # Subject split for this geneset (seeded)
        ctrl_subjects, case_subjects = split_subjects(adata.obs[args.subject_id_col], seed=args.seed + g)

        # Save control once per split
        ctrl_h5 = generate_control_h5ad(
            adata_full = adata
            subject_ctrl = ctrl_subjects,
            split_dir = split_dir,
            subject_id_col = args.subject_id_col,
            slogan = slogan,
            clobber = args.clobber)
        print(f"[geneset {g:03d}] control saved: {ctrl_h5}")


        # CASE base slice (before perturbation)
        case_mask = adata.obs[args.subject_id_col].isin(case_subjects).values
        ad_case_base = adata[case_mask].copy()

        # CONTROL base slice
        ctrl_mask = adata.obs[args.subject_id_col].isin(ctrl_subjects).values
        ad_ctrl = adata[ctrl_mask].copy()

        # === CHANGED: choose n_perturb genes ONLY from eligible set (do not drop any genes)
        rng_gene = np.random.default_rng(args.seed + 1000 * g)
        eligible_indices = np.flatnonzero(eligible_mask_global)
        gene_indices = rng_gene.choice(eligible_indices, size=args.n_perturb, replace=False)

        # Δ vector for this geneset (shuffled copy)
        delta_vec = delta_vec_template.copy()
        rng_gene.shuffle(delta_vec)

        # Save CSV of genes & factors
        csv_path = split_dir / "genes.csv"
        csv_path = csv_path.with_name(csv_path.name.replace(" ", "_"))
        if not csv_path.exists():
            df_csv = pd.DataFrame({
                "gene_id": ad_case_base.var_names[gene_indices].values,
                "log2fc": [float(x) for x in delta_vec],
                "factor": [float(2.0 ** x) for x in delta_vec],
                "mean_expression_control_subjects": \
                    [float(x) for x in np.array(ad_case_base.raw.X[:, gene_indices].mean(axis=0)).ravel()],
                "mean_expression_case_subjects": \
                    [float(x) for x in np.array(adata[case_mask].raw.X[:, gene_indices].mean(axis=0)).ravel()]
            })
            df_csv.to_csv(csv_path, index=False)

        # Validate that each requested cell type exists in CASE base before downsampling
        cell_type_col = args.cell_type_col
        case_counts_by_c = ad_case_base.obs[cell_type_col].value_counts().to_dict()
        for c in args.cell_types:
            if case_counts_by_c.get(str(c), 0) == 0:
                raise RuntimeError(
                    f"[geneset {g:03d}] No cells with {cell_type_col} == '{c}' in CASE cohort BEFORE downsampling."
                )

        # Precompute NB params ONCE for the case cohort for this geneset
        # This avoids building a CSC per worker.
        print(f"[geneset {g:03d}] Estimating NB dispersion parameters...")        
        X_case_csc = ad_case_base.raw.X.tocsc()
        mu, var, theta = estimate_nb_dispersion_gene(X_case_csc)
        del X_case_csc

        # Prepare parallel jobs for all (c, keep_frac, rep)
        jobs = []
        for c in args.cell_types:
            for d in args.keep_fracs:
                if not (0.0 <= d <= 1.0):
                    raise ValueError(f"keep_frac must be in [0,1]; got {d}")
                for rep in range(1, args.num_reps + 1):
                    if d == 1.0 and rep > 1:
                        print(f"[geneset {g:03d}] Warning: skipping redundant rep {rep} with keep_frac=1.0")
                        continue

                    case_h5ad_path = split_dir / f"{slogan}_split_{g:03d}__case__downsample_{c}_{d}__rep_{rep}.h5ad"
                    case_h5ad_path = case_h5ad_path.with_name(case_h5ad_path.name.replace(" ", "_"))
                    if case_h5ad_path.exists() and not args.clobber:
                        print(f"[geneset {g:03d}] Skipping existing case file: {case_h5ad_path}")
                        continue

                    jobs.append((
                        g, split_dir, slogan,
                        ad_case_base, gene_indices, delta_vec, eligible_mask_global,
                        mu, theta, 
                        cell_type_col, c, float(d), rep, args.seed
                    ))

        # Run in parallel
        if procs == 1:
            results = [_case_worker(job) for job in jobs]
        else:
            with Pool(processes=procs) as pool:
                results = pool.map(_case_worker, jobs)

        for case_h5ad_path in results:
            print(f"[geneset {g:03d}] case saved: {case_h5ad_path}")

    print("Done.")


if __name__ == "__main__":
    main()
