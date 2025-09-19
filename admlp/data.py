from __future__ import annotations
import os
from typing import Optional, Tuple, List, Dict
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


def read_h5ad(path: str) -> ad.AnnData:
    """Load an AnnData .h5ad from disk without modifying it."""
    return ad.read_h5ad(path)


def ensure_csr(mat) -> sp.csr_matrix:
    """Ensure cells×genes CSR."""
    if sp.issparse(mat):
        return mat.tocsr()
    return sp.csr_matrix(mat)


def normalize_within_cells(
    X_cxg: sp.csr_matrix,
    target_sum: float = 1e4,
) -> sp.csr_matrix:
    """Library-size normalize each cell (row) to `target_sum` counts.
    Returns a CSR matrix (cells×genes) of float32.
    """
    X = X_cxg.tocsr() if sp.issparse(X_cxg) else sp.csr_matrix(X_cxg)
    cell_sums = np.asarray(X.sum(axis=1)).ravel()
    cell_sums[cell_sums == 0] = 1.0
    scale = (target_sum / cell_sums).astype(np.float32)
    X = X.multiply(scale[:, None]).astype(np.float32)
    return X.tocsr() 


def apply_gene_mask(
    adata: ad.AnnData,
    gene_mask_path: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Resolve the gene mask against adata.var_names.
    Returns (mask_bool_per_var, ordered_gene_names_used).
    If gene_mask_path is None, uses all genes in adata.var_names.
    """
    var_names = np.asarray(adata.var_names.astype(str))
    if gene_mask_path is None:
        mask = np.ones(var_names.shape[0], dtype=bool)
        return mask, var_names

    genes = pd.read_csv(gene_mask_path, header=None, names=["gene"]).gene.astype(str).values
    genes = np.unique(genes)  # dedupe

    # Intersection, keeping adata.var_names order
    mask = np.isin(var_names, set(genes))
    used = var_names[mask]
    return mask, used


def build_X_y_groups(
    adata: ad.AnnData,
    label_col: str,
    subject_col: str,
    target_sum: float = 1e4,
    gene_mask_path: Optional[str] = None,
    log1p: bool = False,
) -> Tuple[sp.csr_matrix, np.ndarray, np.ndarray, List[str]]:
    """Create normalized feature matrix and aligned targets/groups.

    Returns:
      X_cxg_norm_masked: cells×genes sparse CSR
      y_labels: array[str] length n_cells
      groups: array[str] (subject IDs) length n_cells
      genes_used: list[str] in the exact order used for columns
    """
    if label_col not in adata.obs.columns:
        raise KeyError(f"label_col '{label_col}' not found in adata.obs")
    if subject_col not in adata.obs.columns:
        raise KeyError(f"subject_col '{subject_col}' not found in adata.obs")

    X_cxg = ensure_csr(adata.X)
    X_cxg = normalize_within_cells(X_cxg, target_sum=target_sum)
    if log1p:
        X_cxg = X_cxg.copy()
        X_cxg.data = np.log1p(X_cxg.data)

    # ensure sliceable before masking columns
    X_cxg = X_cxg.tocsc()
    X_cxg = X_cxg[:, mask].tocsr()

    y_labels = adata.obs[label_col].astype(str).values
    groups = adata.obs[subject_col].astype(str).values

    return X_cxg, y_labels, groups, list(map(str, genes_used))


def split_subjects(
    groups: np.ndarray,
    test_frac: float = 0.25,
    seed: int = 42,
) -> Tuple[List[str], List[str]]:
    """Split unique subjects into train vs test by fraction.
    Returns (train_subjects, test_subjects) as lists of strings.
    """
    uniq = np.array(sorted(pd.unique(groups)))
    rng = np.random.default_rng(seed)
    rng.shuffle(uniq)
    n_test = max(1, int(np.round(test_frac * len(uniq))))
    test_subjects = uniq[:n_test].tolist()
    train_subjects = uniq[n_test:].tolist()
    if len(train_subjects) == 0:
        # fall back to at least one subject in train
        train_subjects = test_subjects[:1]
        test_subjects = test_subjects[1:]
    return train_subjects, test_subjects
