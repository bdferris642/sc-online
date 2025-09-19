#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import torch
from admlp import data as ZD
from admlp import model as ZM
from admlp import io_utils as ZIO


def parse_args():
    ap = argparse.ArgumentParser(description="Score an AnnData with a trained classifier, output softmax and predictions")
    ap.add_argument("--h5ad", required=True, type=str)
    ap.add_argument("--model", required=True, type=str, help="Path to classifier.pt saved by training script")
    ap.add_argument("--out-csv", required=True, type=str, help="Output CSV path")
    ap.add_argument("--target-sum", type=float, default=None, help="Override normalization target (defaults to training metadata)")
    ap.add_argument("--log1p", action="store_true", help="Apply log1p if used at train time (you should match training)")
    return ap.parse_args()


def main():
    args = parse_args()

    pkg = ZIO.load_classifier(args.model)
    class_names = pkg["class_names"]
    genes_used = pkg["genes_used"]
    train_args = pkg.get("args", {})

    # Use training target_sum unless explicitly overridden
    target_sum = args.target_sum if args.target_sum is not None else float(train_args.get("target_sum", 10000.0))
    log1p = bool(train_args.get("log1p", False)) if not args.log1p else True

    adata = ad.read_h5ad(args.h5ad)
    # Build normalized matrix using only the stored genes (order-critical). Fill missing genes with 0s.
    var_index = pd.Index(adata.var_names.astype(str))
    present_mask = var_index.isin(genes_used)

    X_cxg = ZD.ensure_csr(adata.X)
    X_cxg = ZD.normalize_within_cells(X_cxg, target_sum=target_sum)
    if log1p:
        X_cxg = X_cxg.copy()
        X_cxg.data = np.log1p(X_cxg.data)

    # Build columns in the order of genes_used
    # Start with zeros, then fill columns that exist
    n_cells = X_cxg.shape[0]
    n_genes = len(genes_used)
    X_ordered = sp.csr_matrix((n_cells, n_genes), dtype=np.float32)

    # Map existing genes to their positions
    adata_pos = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    used_pos = {g: j for j, g in enumerate(genes_used)}
    # Stack column slices that exist
    slices = []
    cols = []
    for g, j in used_pos.items():
        i = adata_pos.get(g, None)
        if i is not None:
            slices.append(X_cxg[:, i])
            cols.append(j)
    if slices:
        X_subset = sp.hstack(slices, format="csr")
        # reorder columns to match their intended positions
        perm = np.argsort(cols)
        X_subset = X_subset[:, perm]
        cols_sorted = np.array(cols)[perm]
        # place into X_ordered
        X_ordered = X_ordered.tolil()
        X_ordered[:, cols_sorted] = X_subset
        X_ordered = X_ordered.tocsr()

    # Build model and load weights
    in_dim = X_ordered.shape[1]
    n_classes = len(class_names)
    hidden = train_args.get("hidden", [512, 256])
    dropout = float(train_args.get("dropout", 0.5))
    device = "cuda" if torch.cuda.is_available() else "cpu"

    model = ZM.DropoutMLP(in_dim=in_dim, n_classes=n_classes, hidden=list(map(int, hidden)), dropout=dropout)
    model.load_state_dict(pkg["state_dict"])
    model.eval()
    model.to(device)

    # Inference
    loader = torch.utils.data.DataLoader(ZM.TensorDataset(ZM._to_tensor(X_ordered), torch.zeros(X_ordered.shape[0], dtype=torch.int64)), batch_size=4096, shuffle=False)
    # The above uses a convenience import of TensorDataset through ZM; if IDE complains, it's fine at runtime.

    logits = []
    with torch.no_grad():
        for xb, _ in loader:
            xb = xb.to(device)
            lg = model(xb).cpu().numpy()
            logits.append(lg)
    y_log = np.vstack(logits) if logits else np.zeros((0, n_classes), dtype=np.float32)
    y_prob = ZM.softmax(y_log) if y_log.size else np.zeros((0, n_classes), dtype=np.float32)
    y_pred_idx = y_prob.argmax(axis=1) if y_prob.size else np.array([], dtype=int)
    y_pred = [class_names[i] for i in y_pred_idx]

    # Build output DataFrame
    df = pd.DataFrame({"cell_barcode": adata.obs_names.astype(str)})
    for k, cname in enumerate(class_names):
        df[f"prob_{cname}"] = y_prob[:, k]
    df["prediction"] = y_pred

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print("Wrote:", out_path)


if __name__ == "__main__":
    main()
