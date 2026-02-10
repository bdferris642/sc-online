import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

from data_utils import apply_cnmf_usage_df


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad-path", "-p", required=True, help="Input .h5ad path")
    ap.add_argument("--usages-path", "-u", required=True, help="Input CNMF usages .csv path")
    ap.add_argument("--features-to-plot", "-f", nargs="+", help="List of features to plot",
                    default=["frac_mito", "frac_intronic", "log10_nUMI", "dropsift_frac_contamination"])
    
    args = ap.parse_args()
    # print arguments for reference
    for k, v in vars(args).items():
        print(f"{k}: {v}")
    return args


def main():
    args = parse_args()
    
    # check that usages exists 
    if not os.path.exists(args.usages_path):
        raise FileNotFoundError(f"CNMF usages file not found: {args.usages_path}")
    
    metadata_dir = os.path.join(os.path.dirname(args.usages_path), "qc")
    if not os.path.exists(metadata_dir):
        os.makedirs(metadata_dir)

    print(f"\nqc-cnmf-factors.py:\tLoading input AnnData from: {args.h5ad_path}")
    adata = sc.read_h5ad(args.h5ad_path)

    missing_cols = [col for col in args.features_to_plot if col not in adata.obs.columns]
    if len(missing_cols) > 0:
        raise ValueError(f"Features {missing_cols} not found in adata.obs")
    
    print(f"\nqc-cnmf-factors.py:\tApplying cNMF usages from: {args.usages_path}")
    adata = apply_cnmf_usage_df(adata, cnmf_path=args.usages_path)

    cnmf_cols = sorted(
        [col for col in adata.obs.columns if re.match(r"cnmf.*_k_[0-9]+_gep.*_norm$", col)])
    print(f"\nFound {len(cnmf_cols)} CNMF usage columns in adata.obs")
    if len(cnmf_cols) == 0:
        raise ValueError("No CNMF usage columns found in adata.obs")

    corrs_df = pd.DataFrame(index=cnmf_cols)
    for qc_metric in args.features_to_plot:
        corrs_df[qc_metric] = [
            adata.obs[qc_metric].corr(adata.obs[col], method="pearson")
            for col in cnmf_cols
        ]

    corrs_df["abs_mean"] = corrs_df[args.features_to_plot].abs().mean(axis=1)
    corrs_df.reset_index(names="cnmf_col").to_csv(
        f"{metadata_dir}/cnmf_usage_qc_metric_correlations.csv",
        index=False
    )
    print(f"\nSaved CNMF usage - QC metric correlations to: {metadata_dir}/cnmf_usage_qc_metric_correlations.csv")

    n_qc = len(args.features_to_plot)
    batch_size = 4

    for start in range(0, len(cnmf_cols), batch_size):
        end = min(start + batch_size, len(cnmf_cols))
        batch_cols = cnmf_cols[start:end]

        fig, axs = plt.subplots(
            nrows=len(batch_cols),
            ncols=n_qc,
            figsize=(3 * n_qc, 2 * len(batch_cols)),
            squeeze=False
        )

        for i, cnmf_col in enumerate(batch_cols):
            for j, qc_metric in enumerate(args.features_to_plot):
                ax = axs[i, j]
                ax.scatter(
                    adata.obs[qc_metric],
                    adata.obs[cnmf_col],
                    s=1,
                    alpha=0.02
                )
                ax.set_xlabel(qc_metric if i == len(batch_cols) - 1 else "")
                ax.set_ylabel(cnmf_col if j == 0 else "")

        plt.tight_layout()

        # 1-indexed, zero-padded figure numbering
        fig_path = (
            f"{metadata_dir}/cnmf_usage_qc_metric_scatterplots_"
            f"{start + 1:02d}_to_{end:02d}.png"
        )
        plt.savefig(fig_path, dpi=200)
        plt.close(fig)

        print(f"\nSaved CNMF usage - QC metric scatterplots to: {fig_path}")

    print("Done.")

if __name__ == "__main__":
    main()