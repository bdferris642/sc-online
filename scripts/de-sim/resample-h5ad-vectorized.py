#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simulate subject-level case/control data with gene-set perturbations and cell-type downsampling.

Hierarchy:
  Subject split s (= g): 50/50 (by some subject_id) with fixed RNG 
    ├── control.h5ad  (unperturbed controls for that split)
    └── GeneSet g = 1..G
        ├── save CSV of {gene_id, log2fc, factor}
        └── For (rep in 1..N) × each (c in C) × (keep_frac in D)  
            ├── perturb n-perturb genes with fixed Δ (evenly spaced in [-2, -.1], [.1, 2]), NB resample
            ├── TODO: library balancing of non-perturbed (eligible) genes per cell (∝ counts)
            ├── downsample target cell type c by keep_frac (others unchanged)
            └── save case h5ad
            
Filenames & layout (dirname and slogan extracted from --path):
    {dirname}/simulations/{slogan}/split_{split_idx:03d}/
        {slogan}__control.h5ad
        geneset_{geneset_idx:03d}/
            genes.csv
            {slogan}__case__rep_{rep:03d}__{c}_keep_{k}.h5ad

Notes:
- Eligible genes are computed on the full dataset BEFORE splitting:
    total_counts > min_total_counts_gene_thresh AND frac_cells_expressed > min_frac_expressed_gene_thresh.
  PERTURBED genes must be drawn from this eligible set, BUT NO GENES ARE REMOVED FROM ANY AnnData.
- Perturbations apply to ALL CELLS (not just a specific cell type). Downsampling happens after perturbation.
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
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any, Dict, List, Optional, NamedTuple, Tuple

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# --- Global handle for the CASE CSR used by workers (Linux/fork-friendly) ---
_GLOBAL_X_CASE_CSR: Optional[sparse.csr_matrix] = None


class SubjectJob(NamedTuple):
    subject_index: int
    rows_for_subject: np.ndarray
    gene_indices_perturbed: np.ndarray
    log2fc_perturbed: np.ndarray
    theta_global: np.ndarray
    use_subject_theta: bool
    theta_shrinkage_k: float
    seed: int


# --------------------------
# Utilities
# --------------------------


def add_case_control_flag(adata: ad.AnnData, label: str):
    adata.obs["case_control"] = pd.Categorical(
        [label] * adata.n_obs, categories=["control", "case"]
    )


def assign_gene_effects_to_subjects(
    n_genes: int,
    case_subject_ids: np.ndarray,  # shape (S,) strings
    logfc_levels: np.ndarray,  # shape (L,), includes neg & pos sides together
    frac_levels: np.ndarray,  # shape (F,)
    seed: int,
) -> Tuple[np.ndarray, List[Dict[str, Any]]]:
    """
    Create a per-gene-by-subject Δ matrix and per-gene metadata.

    Returns
    -------
    delta_mat : np.ndarray, shape (n_genes, S)
        Entry [g, s] = Δ (0 if not applied to subject s).
    meta_list : list of dicts, length n_genes
        Each dict has: {'log2fc': float, 'frac_applied': float, 'n_subjects': int, 'subjects_applied': List[str]}
    """
    rng = np.random.default_rng(seed)
    S = int(case_subject_ids.size)
    combos = [
        (float(d), float(f)) for d in logfc_levels for f in frac_levels
    ]  # length K = L*F
    rng.shuffle(combos)  # deterministic shuffle
    K = len(combos)

    # Fill exactly n_genes: cycle combos and take first n_genes of the cycled list
    repeats = (n_genes + K - 1) // K
    pool = (combos * repeats)[:n_genes]

    delta_mat = np.zeros((n_genes, S), dtype=float)
    meta_list: List[Dict[str, Any]] = []

    for g_idx, (delta_val, frac_val) in enumerate(pool):
        n_apply = int(np.clip(round(frac_val * S), 1, S))
        apply_idx = rng.choice(S, size=n_apply, replace=False)
        delta_mat[g_idx, apply_idx] = delta_val
        subj_applied = [str(case_subject_ids[i]) for i in sorted(apply_idx.tolist())]
        meta_list.append(
            {
                "log2fc": float(delta_val),
                "frac_applied": float(frac_val),
                "n_subjects": int(n_apply),
                "subjects_applied": subj_applied,
            }
        )
    return delta_mat, meta_list


def basename_noext(path: str) -> str:
    b = os.path.basename(path)
    if b.endswith(".h5ad"):
        b = b[:-5]
    return b


def build_delta_vector(num_genes: int = 400) -> np.ndarray:
    """
    Build evenly spaced Δ in [-2, -0.1], [0.1, 2]
    """
    grid_neg = np.linspace(-2, -0.1, num_genes // 2)
    grid_pos = np.linspace(0.1, 2, num_genes - len(grid_neg))
    vec = np.concatenate([grid_neg, grid_pos])
    rng = np.random.default_rng(0)
    rng.shuffle(vec)
    return vec


def build_frac_grid(num_levels: int, fmin: float, fmax: float) -> np.ndarray:
    """
    Evenly spaced fractions in [fmin, fmax], shape (num_levels,).
    """
    return np.linspace(float(fmin), float(fmax), int(num_levels))


def build_logfc_grid(
    num_per_sign: int, min_abs: float, max_abs: float, *, include_endpoints: bool = True
) -> np.ndarray:
    """
    Evenly spaced log2FC values, symmetric around 0:
    [-max_abs .. -min_abs] and [min_abs .. max_abs], with `num_per_sign` per side.
    Returns shape (2*num_per_sign,).
    """
    if include_endpoints:
        neg = np.linspace(-max_abs, -min_abs, num_per_sign)
        pos = np.linspace(min_abs, max_abs, num_per_sign)
    else:
        neg = np.linspace(-max_abs, -min_abs, num_per_sign, endpoint=False)
        pos = np.linspace(min_abs, max_abs, num_per_sign, endpoint=False)
    return np.concatenate([neg, pos])


def compute_eligible_gene_mask(
    adata: ad.AnnData, min_total: int, min_frac: float
) -> np.ndarray:
    """
    Return a boolean mask (length n_vars) of genes that pass the thresholds,
    for fraction of cells expressing the gene and total counts.
    computed on raw counts in X across all cells. Leaves adata unchanged.
    """

    X = adata.raw.X if adata.raw is not None else adata.X
    n = adata.n_obs
    if sparse.issparse(X):
        Xc = X.tocsc()
        totals = np.asarray(Xc.sum(axis=0)).ravel()
        n_nonzero = np.asarray((Xc > 0).sum(axis=0)).ravel()
    else:
        totals = X.sum(axis=0).astype(np.float64, copy=False)
        n_nonzero = (X > 0).sum(axis=0)
    frac = n_nonzero / float(n)
    return (totals > min_total) & (frac > min_frac)


def downsample_case_and_save(
    X_resampled_csr: sparse.csr_matrix,
    ad_case_base: ad.AnnData,
    cell_type_col: str,
    target_cell_type: str,
    keep_fraction: float,
    out_path: Path,
    simulation_meta: Dict[str, Any],
    clobber: bool = False,
    seed: int = 12345,
) -> None:
    """
    From a resampled CASE matrix, downsample cells of the given type
    to the specified keep fraction and write an H5AD.

    - Preserves `obs`/`var` from `ad_case_base`.
    - Drops `.raw` to avoid I/O issues with reserved columns (can be rebuilt if needed).
    """
    assert 0.0 <= keep_fraction <= 1.0
    ad_case = ad.AnnData(
        X=X_resampled_csr.copy(),
        obs=ad_case_base.obs.copy(),
        var=ad_case_base.var.copy(),
        uns=ad_case_base.uns.copy(),
    )
    mask_type = ad_case.obs[cell_type_col].astype(str).values == str(target_cell_type)
    n_type = int(mask_type.sum())
    if n_type == 0:
        raise RuntimeError(
            f"No cells with {cell_type_col} == '{target_cell_type}' in CASE cohort before downsampling."
        )

    if keep_fraction < 1.0:
        rng = np.random.default_rng(seed)
        n_keep = int(round(keep_fraction * n_type))
        idx_type = np.flatnonzero(mask_type)
        keep_idx = (
            set(rng.choice(idx_type, size=n_keep, replace=False).tolist())
            if n_keep > 0
            else set()
        )
        final_mask = np.ones(ad_case.n_obs, dtype=bool)
        for idx in idx_type:
            if idx not in keep_idx:
                final_mask[idx] = False
        ad_case = ad_case[final_mask].copy()

    add_case_control_flag(ad_case, "case")
    ad_case.uns["simulation"] = dict(simulation_meta)  # copy
    out_path = out_path.with_name(out_path.name.replace(" ", "_"))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_sanitized_h5ad(adata=ad_case, path=out_path, clobber=clobber, keep_raw=False)


def estimate_nb_dispersion_gene(
    X_csc: sparse.csc_matrix,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    X_sq.data = X_sq.data**2
    ex2 = np.array(X_sq.mean(axis=0)).ravel()

    # --- Variance = E[X²] - (E[X])²
    var = ex2 - mu**2

    # --- Dispersion calculation
    # From NB variance formula: Var = μ + μ²/θ
    # Rearranging: θ = μ² / (Var - μ)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = var - mu
        theta_calc = np.where(denom > 0, (mu**2) / denom, np.inf)

    # --- Replace invalid/infinite values
    # - If variance <= mean, denom <= 0 → implies ~Poisson behavior
    #   → set θ very large (1e12) so NB ~ Poisson.
    # - Ensure no -inf, nan, or zero.
    theta = np.where(np.isfinite(theta_calc), theta_calc, 1e12)
    theta = np.maximum(theta, 1e-8)  # guard against zero/negative

    return mu, var, theta


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
    library_sizes = np.array(ad_ctrl.X.sum(axis=1)).ravel()
    ad_ctrl.uns["simulation"] = {
        "split": int(split_dir.name.split("_")[-1]),
        "case_control": "control",
        "cell_type": None,
        "keep_frac": None,
        "replicate": None,
        "perturbed_genes": {},
        "library_sizes": library_sizes,
        "perturbation": None,
    }
    # if there is an adata.raw, replace adata.X with adata.raw.X
    if ad_ctrl.raw is not None:
        X = ad_ctrl.raw.X
        if sparse.issparse(X):
            X = X.astype(np.float64)
        else:
            X = X.astype(np.float64, copy=False)
        ad_ctrl.X = X  # replace main matrix
        ad_ctrl.var = ad_ctrl.raw.var.copy()  # make sure genes match raw
    else:
        print("Warning: adata.raw is None; using adata.X as-is for control.")

    out = split_dir / f"{slogan}__control.h5ad"
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


def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate subject-level case/control H5ADs with gene perturbations and downsampling."
    )
    p.add_argument(
        "--path", required=True, help="Path to input anndata .h5ad (raw counts in X)."
    )
    p.add_argument(
        "--n-perturbed-genes",
        type=int,
        default=200,
        help="Number of genes to perturb per gene set.",
    )
    p.add_argument(
        "--cell-type-col",
        required=True,
        default="cell_type",
        help="obs column name containing cell types.",
    )
    p.add_argument(
        "--subject-id-col",
        required=True,
        default="donor",
        help="obs column name containing subject IDs (e.g. 'subject.subjectGuid', 'donor', etc.).",
    )
    p.add_argument(
        "--num-splits", type=int, default=1, help="Number of case/control splits."
    )
    p.add_argument(
        "--num-gene-sets",
        type=int,
        default=1,
        help="Number of perturbed gene sets (one split per gene set).",
    )
    p.add_argument(
        "--num-reps", type=int, default=1, help="Number of replicates per gene set."
    )
    p.add_argument(
        "--cell-types", nargs="+", required=True, help="List of target cell types C."
    )
    p.add_argument(
        "--keep-fracs",
        type=float,
        nargs="+",
        required=True,
        help="List of keep fractions D (0..1). Example: 1.0 0.5 0.25 0.0",
    )
    p.add_argument(
        "--min-total-counts-gene-thresh",
        type=int,
        default=200,
        help="Minimum total counts across all cells to consider a gene ELIGIBLE for perturbation (no matrix filtering).",
    )
    p.add_argument(
        "--min-frac-expressed-gene-thresh",
        type=float,
        default=0.01,
        help="Minimum fraction of cells with count>0 to consider a gene ELIGIBLE for perturbation (no matrix filtering).",
    )
    p.add_argument("--seed", type=int, default=12345, help="Base RNG seed.")
    p.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of parallel worker processes for case generation (default: min(available, 40)).",
    )
    p.add_argument(
        "--clobber",
        action="store_true",
        help="Overwrite existing input h5ad if subject ID must be created or outputs already exist.",
    )
    p.add_argument(
        "--subject-theta",
        action="store_true",
        help="Estimate per-subject dispersion (θ) with shrinkage toward global; "
        "otherwise use global θ for all subjects.",
    )
    p.add_argument(
        "--theta-shrinkage-k",
        type=float,
        default=0,
        help="Shrinkage strength toward global dispersion when --subject-theta is set "
        "(larger K => more shrinkage).",
    )
    p.add_argument(
        "--rebalance",
        action="store_true",
        help="CURRENTLY A NO-OP. TODO: Use iterative multinomial rebalancing to adjust library sizes per cell. ",
    )
    p.add_argument(
        "--num-logfc-per-sign",
        type=int,
        default=5,
        help="Number of evenly spaced log2FC levels on each side (negative and positive).",
    )
    p.add_argument(
        "--logfc-min-abs",
        type=float,
        default=0.1,
        help="Minimum absolute log2FC (exclusive of 0).",
    )
    p.add_argument(
        "--logfc-max-abs", type=float, default=0.2, help="Maximum absolute log2FC."
    )
    p.add_argument(
        "--num-frac-levels",
        type=int,
        default=10,
        help="Number of evenly spaced fractions in [0.1, 1.0] for fraction of case subjects perturbed.",
    )
    p.add_argument(
        "--frac-min",
        type=float,
        default=0.1,
        help="Minimum fraction of case subjects to perturb (inclusive).",
    )
    p.add_argument(
        "--frac-max",
        type=float,
        default=1.0,
        help="Maximum fraction of case subjects to perturb (inclusive).",
    )

    # p.add_argument("--subject-lfc-sd", type=float, default=0,
    #                help="Std dev of subject-level log2FC around the global Δ for each perturbed gene. Default 0 (no variation).")

    return p.parse_args()


def rebalance_counts() -> None:
    pass
    """
    To-Do: In-place rebalancing of a cell's library size.
    """


def resample_case_once_subject_parallel(
    ad_case_base: ad.AnnData,
    gene_indices_perturbed: np.ndarray,
    delta_mat: np.ndarray,
    eligible_mask_global: np.ndarray,
    subject_id_col: str,
    theta_global: np.ndarray,
    use_subject_theta: bool,
    theta_shrinkage_k: float,
    seed_base: int,
    procs: int,
) -> sparse.csr_matrix:
    """
    Resample the CASE cohort **once** for a given (gene set, replicate),
    by parallelizing over subjects. Each worker:
      - computes that subject's μ/var/θ (optionally shrunk to global),
      - draws NB counts for the perturbed genes.

    The parent assembles the full (cells x genes) matrix by replacing
    the perturbed columns with the sampled values. (Rebalancing is a
    no-op placeholder for now.)

    Returns
    -------
    X_new_csr : csr_matrix
        Resampled CASE counts matrix (perturbed genes replaced).
    """
    # Original matrices in convenient formats
    X_case_csr = ad_case_base.X.tocsr(copy=True)
    n_cells, n_genes = X_case_csr.shape
    n_perturb = int(gene_indices_perturbed.size)

    # Subject coding
    subjects_case = ad_case_base.obs[subject_id_col].astype(str).values
    subj_codes, subj_uniqs = pd.factorize(subjects_case, sort=True)
    n_subjects = len(subj_uniqs)

    # Build jobs: one per subject
    rows_by_subject: List[np.ndarray] = [
        np.where(subj_codes == s)[0] for s in range(n_subjects)
    ]
    jobs = [
        SubjectJob(
            subject_index=s,
            rows_for_subject=rows_by_subject[s],
            gene_indices_perturbed=gene_indices_perturbed,
            log2fc_perturbed=delta_mat[:, s],
            theta_global=theta_global,
            use_subject_theta=bool(use_subject_theta),
            theta_shrinkage_k=float(theta_shrinkage_k),
            seed=int(seed_base + s),
        )
        for s in range(n_subjects)
    ]

    # Run subject workers in parallel; attach the CSR once via initializer
    worker_procs = min(max(1, procs), n_subjects)
    if worker_procs == 1:
        _init_subject_worker(X_case_csr)
        results = [_subject_resample_worker(job) for job in jobs]
    else:
        with Pool(
            processes=worker_procs,
            initializer=_init_subject_worker,
            initargs=(X_case_csr,),
        ) as pool:
            results = pool.map(_subject_resample_worker, jobs)

    # Assemble dense block for perturbed genes: (n_cells x n_perturb)
    new_counts_perturbed = np.empty((n_cells, n_perturb), dtype=np.int32)
    for s_idx, rows_s, block_s in results:
        new_counts_perturbed[rows_s, :] = block_s

    # (Rebalancing placeholder goes here; for now, we only replace perturbed columns.)

    # Assemble final sparse matrix by rebuilding rows (LIL → CSR)
    X_lil = X_case_csr.tolil(copy=True)

    # For speed, prepare a set for perturbed indices
    perturbed_set = set(map(int, gene_indices_perturbed))

    # Replace perturbed columns row-wise (fast enough; avoids column-wise costly ops)
    for i in range(n_cells):
        cols_i = X_lil.rows[i]
        vals_i = X_lil.data[i]

        # Keep only non-perturbed entries from the original row
        keep_cols: List[int] = []
        keep_vals: List[int] = []
        for c0, v0 in zip(cols_i, vals_i):
            if c0 not in perturbed_set and v0 != 0:
                keep_cols.append(int(c0))
                keep_vals.append(int(v0))

        # Add perturbed entries from the new dense block
        for j_local, gidx in enumerate(gene_indices_perturbed):
            v = int(new_counts_perturbed[i, j_local])
            if v != 0:
                keep_cols.append(int(gidx))
                keep_vals.append(v)

        X_lil.rows[i] = keep_cols
        X_lil.data[i] = keep_vals

    return X_lil.tocsr()


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
        raise FileExistsError(
            f"File already exists: {path}. Use --clobber to overwrite."
        )

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


def _init_subject_worker(X_case_csr: sparse.csr_matrix) -> None:
    """
    Pool initializer. Attaches the CASE cohort CSR matrix to a module-global
    so each worker can slice its subject's rows without pickling the whole matrix.
    (On Linux 'fork', memory is shared copy-on-write.)
    """
    global _GLOBAL_X_CASE_CSR
    _GLOBAL_X_CASE_CSR = X_case_csr


def _sanitize_df(df: pd.DataFrame, where: str) -> None:
    # AnnData reserves the column name "_index" in obs/var/raw.var
    if "_index" in df.columns:
        df.rename(columns={"_index": f"{where}_index_orig"}, inplace=True)
        df.set_index(f"{where}_index_orig", inplace=True, drop=False)


def _subject_resample_worker(job: SubjectJob) -> Tuple[int, np.ndarray, np.ndarray]:
    subject_index = job.subject_index
    rows_for_subject = job.rows_for_subject
    gene_indices_perturbed = job.gene_indices_perturbed
    log2fc_perturbed = job.log2fc_perturbed
    theta_global = job.theta_global
    use_subject_theta = job.use_subject_theta
    theta_shrinkage_k = job.theta_shrinkage_k
    seed = job.seed
    """
    Compute per-subject μ/var/θ (optionally shrink to global) and sample NB counts
    for the perturbed genes ONLY, returning a dense block for this subject.

    Parameters
    ----------
    args :
      subject_index : int
          0..(n_subjects-1).
      rows_for_subject : np.ndarray, shape (n_cells_subject,)
          Absolute row indices of the CASE matrix belonging to this subject.
      gene_indices_perturbed : np.ndarray, shape (n_perturb,)
          Column indices of perturbed genes (global var index).
      log2fc_perturbed : np.ndarray, shape (n_perturb,)
          Δ per perturbed gene (same order).
      theta_global : np.ndarray, shape (n_genes,)
          Global θ per gene for shrinkage or fallback.
      use_subject_theta : bool
          If True, estimate subject θ and shrink toward global.
      theta_shrinkage_k : float
          Shrinkage constant K (larger ⇒ more shrinkage).
      seed : int
          RNG seed (subject-specific for determinism).

    Returns
    -------
    subject_index : int
        Echo of the subject index.
    rows_for_subject : np.ndarray
        Echo of the row indices (used by the parent to place results).
    sampled_block : np.ndarray, shape (n_cells_subject, n_perturb), dtype=int32
        NB-sampled counts for this subject’s cells and perturbed genes.
    """
    assert _GLOBAL_X_CASE_CSR is not None, "Worker CSR not initialized."

    rng = np.random.default_rng(seed)

    # Slice this subject’s rows once
    Xs_csr = _GLOBAL_X_CASE_CSR[rows_for_subject, :]  # (cells_s x genes)
    n_cells_s, n_genes = Xs_csr.shape
    n_perturb = int(gene_indices_perturbed.size)

    # ---- per-subject means μ_s(g) ----
    # μ = mean over cells (use CSR → CSC squared trick for moments if needed)
    # Mean:
    mu_s = np.asarray(Xs_csr.mean(axis=0)).ravel()  # (n_genes,)

    # ---- per-subject dispersion θ_s(g) (optional) ----
    if use_subject_theta:
        # second moment E[X^2]
        Xs2 = Xs_csr.copy()
        Xs2.data **= 2
        ex2_s = np.asarray(Xs2.mean(axis=0)).ravel()
        var_s = ex2_s - mu_s**2
        denom = var_s - mu_s
        with np.errstate(divide="ignore", invalid="ignore"):
            theta_s = np.where(denom > 0, (mu_s**2) / denom, np.inf)
        theta_s = np.where(np.isfinite(theta_s), theta_s, 1e12)
        theta_s = np.maximum(theta_s, 1e-8)

        # shrink to global on precision scale
        phi_global = 1.0 / np.maximum(theta_global, 1e-8)
        n_s = max(1, n_cells_s)
        w = n_s / (n_s + float(theta_shrinkage_k))
        phi_s = 1.0 / np.maximum(theta_s, 1e-8)
        phi_shrunk = w * phi_s + (1.0 - w) * phi_global
        theta_effective = 1.0 / np.maximum(phi_shrunk, 1e-8)
    else:
        theta_effective = theta_global

    # ---- sample perturbed genes for this subject ----
    # Prepare output block: (cells_s x n_perturb)
    old_block = Xs_csr[:, gene_indices_perturbed].toarray().astype(np.int32, copy=False)
    sampled_block = old_block.copy()

    deltas_subject = np.asarray(log2fc_perturbed, dtype=float)
    nonzero_cols = np.flatnonzero(deltas_subject != 0.0)
    if nonzero_cols.size == 0:
        return (subject_index, rows_for_subject, sampled_block)

    # Sample only where Δ ≠ 0
    for j_local in nonzero_cols:
        gidx = int(gene_indices_perturbed[j_local])
        mu_prime = (2.0 ** float(deltas_subject[j_local])) * float(mu_s[gidx])
        theta_j = float(theta_effective[gidx])
        if mu_prime <= 0.0 or n_cells_s == 0:
            sampled_block[:, j_local] = 0
        else:
            lam = rng.gamma(shape=theta_j, scale=(mu_prime / theta_j), size=n_cells_s)
            sampled_block[:, j_local] = rng.poisson(lam).astype(np.int32, copy=False)

    return (subject_index, rows_for_subject, sampled_block)


# --------------------------


def main():
    args = parse_args()

    base_path = args.path
    slogan = basename_noext(base_path)

    print("Loading data...")
    adata = ad.read_h5ad(base_path, backed=None)  # in-memory

    # Validate required columns
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"--cell-type-col '{args.cell_type_col}' not found in obs.")
    if args.subject_id_col not in adata.obs.columns:
        raise ValueError(f"--subject-id-col '{args.subject_id_col}' not found in obs.")
    if args.keep_fracs is None or len(args.keep_fracs) == 0:
        raise ValueError("At least one --keep-fracs value required.")
    if any((d < 0.0) or (d > 1.0) for d in args.keep_fracs):
        raise ValueError("all --keep-fracs values must be in [0, 1].")

    # Fail fast if any requested cell type is missing
    missing = set(args.cell_types) - set(adata.obs[args.cell_type_col].unique())
    if missing:
        raise ValueError(
            f"Cell types not found in obs['{args.cell_type_col}']: {missing}"
        )

    # Eligibility mask over all genes (boolean mask of length n-genes)
    eligible_mask_global = compute_eligible_gene_mask(
        adata,
        min_total=args.min_total_counts_gene_thresh,
        min_frac=args.min_frac_expressed_gene_thresh,
    )
    n_eligible = int(eligible_mask_global.sum())
    if n_eligible < args.n_perturbed_genes:
        raise ValueError(
            f"Eligible genes ({n_eligible}) < {args.n_perturbed_genes}. Increase data or relax thresholds."
        )

    # Coerce key obs cols to str
    adata.obs[args.cell_type_col] = adata.obs[args.cell_type_col].astype(str)
    adata.obs[args.subject_id_col] = adata.obs[args.subject_id_col].astype(str)

    # Processes
    procs = min(cpu_count(), 40) if args.processes is None else max(1, args.processes)

    for split_idx in range(1, args.num_splits + 1):
        dirname = os.path.dirname(base_path)
        slogan = basename_noext(base_path)
        split_dir = (
            Path(dirname) / "simulations" / f"{slogan}" / f"split_{split_idx:03d}"
        )
        split_dir.mkdir(parents=True, exist_ok=True)

        # Subject split for this geneset (seeded)
        ctrl_subjects, case_subjects = split_subjects(
            adata.obs[args.subject_id_col], seed=args.seed + split_idx
        )

        # Save control once per split (your helper can verify same subjects unless --clobber)
        ctrl_h5 = generate_control_h5ad(
            adata_full=adata,
            subject_ctrl=ctrl_subjects,
            split_dir=split_dir,
            subject_id_col=args.subject_id_col,
            slogan=slogan,
            clobber=args.clobber,
        )
        print(f"[split {split_idx:03d}] control saved: {ctrl_h5}")

        # CASE base slice (before perturbation)
        print(f"[split {split_idx:03d}] Estimating NB dispersion parameters...")
        case_mask = adata.obs[args.subject_id_col].isin(case_subjects).values
        ad_case_base = adata[case_mask].copy()
        X_case_csc = ad_case_base.X.tocsc()
        subjects_case = ad_case_base.obs[args.subject_id_col].astype(str).values
        subj_codes, subj_uniqs = pd.factorize(subjects_case, sort=True)

        # Compute global NB params ONCE for case cohort
        mu_global, var_global, theta_global = estimate_nb_dispersion_gene(X_case_csc)

        for geneset_idx in range(1, args.num_gene_sets + 1):
            # Pick args.n_perturbed_genes perturbed genes ONLY from eligible set (global indexing)
            rng_gene = np.random.default_rng(args.seed + geneset_idx)
            eligible_indices = np.flatnonzero(eligible_mask_global)
            gene_indices = rng_gene.choice(
                eligible_indices, size=args.n_perturbed_genes, replace=False
            )

            # Case-subject IDs (stable order matches factorization used later)
            case_subject_ids = (
                ad_case_base.obs[args.subject_id_col].astype(str).unique()
            )
            case_subject_ids.sort()  # to match factorize(sort=True)

            # Build grids and per-gene, per-subject Δ matrix (deterministic)
            logfc_levels = build_logfc_grid(
                args.num_logfc_per_sign, args.logfc_min_abs, args.logfc_max_abs
            )
            frac_levels = build_frac_grid(
                args.num_frac_levels, args.frac_min, args.frac_max
            )

            delta_mat, gene_meta_list = assign_gene_effects_to_subjects(
                n_genes=args.n_perturbed_genes,
                case_subject_ids=case_subject_ids,
                logfc_levels=logfc_levels,
                frac_levels=frac_levels,
                seed=args.seed + (split_idx + 1) + 2 * (geneset_idx + 1),
            )

            # Write genes.csv (including per-gene metadata + subject IDs)
            geneset_dir = split_dir / f"geneset_{geneset_idx:03d}"
            geneset_dir.mkdir(parents=True, exist_ok=True)
            genes_csv_path = geneset_dir / "genes.csv"
            if (not genes_csv_path.exists()) or args.clobber:
                df_genes = pd.DataFrame(
                    {
                        "gene_id": ad_case_base.var_names[gene_indices].values,
                        "log2fc": [m["log2fc"] for m in gene_meta_list],
                        "frac_applied": [m["frac_applied"] for m in gene_meta_list],
                        "n_subjects_applied": [m["n_subjects"] for m in gene_meta_list],
                        "subjects_applied": [
                            ",".join(m["subjects_applied"]) for m in gene_meta_list
                        ],
                    }
                )
                df_genes.to_csv(genes_csv_path, index=False)

            # Validate that each requested cell type exists in CASE base before downsampling
            cell_type_col = args.cell_type_col
            case_counts_by_c = ad_case_base.obs[cell_type_col].value_counts().to_dict()
            for c in args.cell_types:
                if case_counts_by_c.get(str(c), 0) == 0:
                    raise RuntimeError(
                        f"[geneset {geneset_idx:03d}], split {split_idx:03d} \
                        No cells with {cell_type_col} == '{c}' in CASE cohort BEFORE downsampling."
                    )

            # === Replicates mean: RESAMPLE ONCE, then fan out over (c, d) ===
            for rep in range(1, args.num_reps + 1):
                print(
                    f"[geneset {geneset_idx:03d}] Subject-parallel resampling (rep {rep:03d})..."
                )

                # get the resampled case matrix
                case_resampled_csr = resample_case_once_subject_parallel(
                    ad_case_base=ad_case_base,
                    gene_indices_perturbed=gene_indices,
                    delta_mat=delta_mat,
                    eligible_mask_global=eligible_mask_global,  # currently unused (rebalance is a no-op)
                    subject_id_col=args.subject_id_col,
                    theta_global=theta_global,
                    use_subject_theta=bool(args.subject_theta),
                    theta_shrinkage_k=float(args.theta_shrinkage_k),
                    seed_base=args.seed + (2 * (geneset_idx + 1)) + (3 * (rep + 1)),
                    procs=procs,
                )

                # Fan out over cell types / keep fracs (cheap; can parallelize if desired)
                for ic, c in enumerate(args.cell_types):
                    for ik, k in enumerate(args.keep_fracs):
                        out_path = (
                            geneset_dir
                            / f"{slogan}__case__rep_{rep:03d}__{c}_keep_{k}.h5ad"
                        )
                        if out_path.exists() and not args.clobber:
                            print(f"Skipping existing case file: {out_path}")
                            continue

                        sim_meta = {
                            "split": split_idx,
                            "geneset": geneset_idx,
                            "case_control": "case",
                            "cell_type": c,
                            "keep_frac": float(k),
                            "replicate": int(rep),
                            "perturbed_genes": {
                                gname: {
                                    "log2fc": float(
                                        df_genes.log2fc[df_genes.gene_id == gname]
                                    ),
                                    "factor": float(
                                        2.0
                                        ** df_genes.log2fc[df_genes.gene_id == gname]
                                    ),
                                    "frac_applied": float(
                                        df_genes.frac_applied[df_genes.gene_id == gname]
                                    ),
                                    "n_subjects_applied": int(
                                        df_genes.n_subjects_applied[
                                            df_genes.gene_id == gname
                                        ]
                                    ),
                                    "subjects_applied": df_genes.subjects_applied[
                                        df_genes.gene_id == gname
                                    ].split(",")
                                    if pd.notna(
                                        df_genes.subjects_applied[
                                            df_genes.gene_id == gname
                                        ]
                                    )
                                    else [],
                                }
                                for gname in ad_case_base.var_names[gene_indices].values
                            },
                        }

                        downsample_case_and_save(
                            X_resampled_csr=case_resampled_csr,
                            ad_case_base=ad_case_base,
                            cell_type_col=args.cell_type_col,
                            target_cell_type=c,
                            keep_fraction=float(k),
                            out_path=out_path,
                            simulation_meta=sim_meta,
                            clobber=args.clobber,
                            seed=args.seed
                            + (split_idx + 1)
                            + (2 * (geneset_idx + 1))
                            + (3 * (rep + 1))
                            + (4 * (ic + 1))
                            + (5 * (ik + 1)),
                        )
    print("Done.")


if __name__ == "__main__":
    main()
