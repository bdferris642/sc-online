#!/usr/bin/env python3
"""
simulate_data.py  –  generate a small synthetic AnnData .h5ad for testing
                     pairwise_dist_hist.py

Creates 5 cell types with distinct cluster centres in a 50-dim latent space.
TypeA–TypeD have well-separated centres; TypeE is drawn from the exact same
distribution as TypeA (same centre, same noise scale), making A↔E distances
indistinguishable from within-A — a guaranteed unimodal cross-type pair for testing.
PCA (20 dims) is pre-computed and stored in obsm['X_pca'].

Usage:
    python simulate_data.py [--out test_data.h5ad] [--n-cells 3000] [--seed 0]
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Re-exec with adata-venv Python if needed
_VENV_DIR    = Path(__file__).parent.parent / "adata-venv"
_VENV_PYTHON = _VENV_DIR / "bin" / "python3"
if Path(sys.prefix).resolve() != _VENV_DIR.resolve():
    if not _VENV_PYTHON.is_file():
        sys.exit(
            f"ERROR: adata-venv not found at {_VENV_DIR}.\n"
            f"Run setup first:  bash {Path(__file__).parent / 'setup_env.sh'}"
        )
    os.execv(str(_VENV_PYTHON), [str(_VENV_PYTHON)] + sys.argv)

import argparse

import anndata
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def simulate(n_cells: int = 3000, seed: int = 0) -> anndata.AnnData:
    rng = np.random.default_rng(seed)

    # 5 cell types with intentionally unequal sizes (tests balancing logic).
    # TypeE shares TypeA's centre with smaller noise — guaranteed unimodal A↔E.
    cell_types  = ["TypeA", "TypeB", "TypeC", "TypeD", "TypeE"]
    proportions = [0.35,    0.25,    0.15,    0.10,    0.15]
    n_per_type  = [max(1, int(p * n_cells)) for p in proportions]
    n_per_type[-1] = n_cells - sum(n_per_type[:-1])   # fix rounding

    n_dims = 50   # raw feature space
    n_pca  = 20   # PCA dims stored in obsm

    # 4 distinct cluster centres (TypeA and TypeE share centres[0])
    centres     = rng.standard_normal((4, n_dims)) * 4
    noise_scales = [1.0, 1.0, 1.0, 1.0, 1.0]
    type_centres = [centres[0], centres[1], centres[2], centres[3], centres[0]]

    X_parts = []
    labels  = []
    for ct, n, centre, scale in zip(cell_types, n_per_type, type_centres, noise_scales):
        noise = rng.standard_normal((n, n_dims)) * scale
        X_parts.append(centre + noise)
        labels.extend([ct] * n)

    X = np.vstack(X_parts).astype(np.float32)

    obs = pd.DataFrame(
        {"cell_type": labels},
        index=[f"cell_{i}" for i in range(n_cells)],
    )

    adata = anndata.AnnData(X=X, obs=obs)

    # Pre-compute PCA and store as obsm['X_pca']
    pca = PCA(n_components=n_pca, random_state=seed)
    adata.obsm["X_pca"] = pca.fit_transform(X).astype(np.float32)

    return adata


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out",     default="test_data.h5ad",
                   help="Output .h5ad path (default: test_data.h5ad)")
    p.add_argument("--n-cells", type=int, default=3000,
                   help="Total number of cells (default: 3000)")
    p.add_argument("--seed",    type=int, default=0,
                   help="Random seed (default: 0)")
    args = p.parse_args()

    print(f"Simulating {args.n_cells} cells …")
    adata = simulate(n_cells=args.n_cells, seed=args.seed)
    print(f"  obs shape : {adata.obs.shape}")
    print(f"  X shape   : {adata.X.shape}")
    print(f"  X_pca     : {adata.obsm['X_pca'].shape}")
    print(f"  Cell types: {dict(adata.obs['cell_type'].value_counts())}")

    out = Path(args.out)
    adata.write_h5ad(out)
    size_mb = out.stat().st_size / 1e6
    print(f"Saved → {out.resolve()}  ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()
