
#!/usr/bin/env python3
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from admlp import data as ZD
from admlp import model as ZM
from admlp import io_utils as ZIO


def parse_args():
    ap = argparse.ArgumentParser(description="Train a dropout MLP on an AnnData .h5ad")
    ap.add_argument("--h5ad", required=True, type=str)
    ap.add_argument("--out-dir", required=True, type=str)
    ap.add_argument("--label-col", required=True, type=str, help="Column in adata.obs with labels (e.g., Final_Zone_Assignments)")
    ap.add_argument("--subject-col", required=True, type=str, help="Column in adata.obs with subject IDs (e.g., donor_id)")
    ap.add_argument("--test-frac", type=float, default=0.25, help="Fraction of unique subjects for held-out test")
    ap.add_argument("--cval-k", type=int, default=5, help="K for GroupKFold within training subjects")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--target-sum", type=float, default=10000.0, help="Per-cell normalization target sum")
    ap.add_argument("--gene-mask", type=str, default=None, help="Path to line-separated gene list to use")
    ap.add_argument("--log1p", action="store_true", help="Apply log1p after normalization")

    # MLP hyperparams
    ap.add_argument("--hidden", type=int, nargs="+", default=[512, 256], help="Hidden layer sizes")
    ap.add_argument("--dropout", type=float, default=0.5)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=2048)
    ap.add_argument("--lr", type=float, default=1e-3)
    return ap.parse_args()


def main():
    args = parse_args()

    adata = ZD.read_h5ad(args.h5ad)
    X, y_labels, groups, genes_used = ZD.build_X_y_groups(
        adata,
        label_col=args.label_col,
        subject_col=args.subject_col,
        target_sum=args.target_sum,
        gene_mask_path=args.gene_mask,
        log1p=args.log1p,
    )

    # Subject split
    train_subjects, test_subjects = ZD.split_subjects(groups, test_frac=args.test_frac, seed=args.seed)

    # Train + CV + Test
    cv_summary, cv_per_fold, test_metrics, le, state = ZM.train_eval_single_split(
        X_cxg=X,
        y_labels=y_labels,
        groups=groups,
        train_subjects=train_subjects,
        hidden=list(args.hidden),
        dropout=args.dropout,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        seed=args.seed,
    )

    # Persist artifacts
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ZIO.save_artifacts(
        out_dir=str(out_dir),
        model_state=state,
        class_names=le.classes_.tolist(),
        genes_used=genes_used,
        train_subjects=train_subjects,
        test_subjects=test_subjects,
        args_dict={
            "h5ad": args.h5ad,
            "label_col": args.label_col,
            "subject_col": args.subject_col,
            "test_frac": args.test_frac,
            "cval_k": args.cval_k,
            "seed": args.seed,
            "target_sum": args.target_sum,
            "gene_mask": args.gene_mask,
            "log1p": args.log1p,
            "hidden": list(map(int, args.hidden)),
            "dropout": float(args.dropout),
            "epochs": int(args.epochs),
            "batch_size": int(args.batch_size),
            "lr": float(args.lr),
        },
        cv_summary=cv_summary,
        cv_per_fold=cv_per_fold,
        test_metrics=test_metrics,
    )

    # Also dump a small human-readable summary
    with open(out_dir / "README.txt", "w") as f:
        f.write("Artifacts saved to this directory.\n")
        f.write("- classifier.pt : model+preproc metadata\n")
        f.write("- metadata.json : subjects/args snapshot\n")
        f.write("- metrics_cv.json : cross-validation summary and per-fold\n")
        f.write("- metrics_test.json : held-out test metrics on remaining subjects\n")

    print("Training complete. Artifacts written to:", out_dir)


if __name__ == "__main__":
    main()