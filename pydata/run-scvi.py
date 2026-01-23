"""
Script to run scVI on an input AnnData object (optionally with batch correction), 
performing clustering and UMAP on that embedding,
and save the output AnnData without overwriting any default layers / obsm / cluster names

Assumes the presence of a 'counts' layer in the input AnnData for scVI training.
And a norm_log layer which is returned to adata.X in the final output.
"""


import argparse
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import scvi
from data_utils import add_back_og_counts, find_markers, is_integer_valued


import torch
torch.set_num_threads(68)
torch.set_num_interop_threads(4)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-path", "-i", required=True, help="Input .h5ad path")
    ap.add_argument("--output-path", "-o", required=False, help="Output .h5ad path. If not set, overwrites input.", default=None)
    ap.add_argument("--n-latent", type=int, default=20)
    ap.add_argument("--n-neighbors", type=int, default=20)
    ap.add_argument("--cluster-name-prefix", type=str, default="scvi")
    ap.add_argument(
        "--og-counts-path",
        type=str,
        default=None,
        help="Path to original counts .h5ad file to restore counts from. \
            If not set, uses counts from input .h5ad.",
    )
    ap.add_argument(
        "--resolutions", "-r", nargs="+", type=float, default=[0.2, 0.5, 0.8]
    )
    ap.add_argument("--categorical-batch-keys", nargs="+", default=None, 
        help="List of obs columns to use as categorical covariate keys for scVI.")
    ap.add_argument("--continuous-covariate-keys", nargs="+", default=None, 
        help="List of obs columns to use as continuous covariate keys for scVI.")
    ap.add_argument(
        "--find-markers",
        action="store_true",
        default=False,
        help="Whether to find marker genes for each cluster after clustering."
    )
    ap.add_argument("--convert-counts-to-dense", action="store_true", default=False,
        help="Whether to convert counts layer to dense before passing to scVI. True is faster but takes more memory." \
        "See https://discourse.scverse.org/t/increase-scvi-integration-speed/1772/3"
    )
    ap.add_argument('--batch-size', type=int, default=1024,
        help='Batch size for scVI training. Larger batch sizes are faster but take more memory.'
    )
    ap.add_argument('--max-epochs', type=int, default=100,
        help='Maximum number of epochs for scVI training.'
    )
    ap.add_argument('--early-stopping-patience', type=int, default=8,
        help='Number of epochs with no improvement on validation set before early stopping.'
    )

    args = ap.parse_args()

    if args.output_path is None:
        args.output_path = args.input_path

    # print arguments for reference
    for k, v in vars(args).items():
        print(f"{k}: {v}")

    return args


def main():
    args = parse_args()

    print(f"\nrun-scvi.py:\tLoading input AnnData from: {args.input_path}")
    adata = sc.read_h5ad(args.input_path)
    print(f"\nrun-scvi.py:\tNumber of cells: {adata.n_obs}, Number of genes: {adata.n_vars}")

    if not "norm_log" in adata.layers:
        raise ValueError("Input AnnData must have a 'norm_log' layer for final output!")
    
    if args.og_counts_path is not None:

        print(f"\nrun-scvi.py:\tLoading ORIGINAL AnnData from {args.og_counts_path}...")
        og_adata = sc.read_h5ad(args.og_counts_path)
        print(f"\nrun-scvi.py:\tRestoring original counts to adata object...")
        adata = add_back_og_counts(adata, og_adata)
        print(
            f"\nrun-scvi.py:\tNumber of cells after restoring counts: {adata.n_obs}, Number of genes: {adata.n_vars}"
        )
        del og_adata  # free memory
    
    if not "counts" in adata.layers:
        raise ValueError("Input AnnData must have a 'counts' layer for scVI!")
    if not is_integer_valued(adata.layers["counts"]):
        raise ValueError("Layer 'counts' must be integer valued for scVI!")
    
    if args.convert_counts_to_dense:
        # it ought not be dense already, but check just in case
        if sp.issparse(adata.layers["counts"]):
            adata.layers["counts"] = (
                adata.layers["counts"]
                .toarray()
                .astype("float32", copy=False)
            )


    missing_cat_keys = []
    if args.categorical_batch_keys is not None:
        for key in args.categorical_batch_keys:
            if key not in adata.obs:
                missing_cat_keys.append(key)
        if len(missing_cat_keys) > 0:
            raise ValueError(f"The following categorical batch keys are missing from adata.obs: {missing_cat_keys}")
    missing_cont_keys = []
    if args.continuous_covariate_keys is not None:
        for key in args.continuous_covariate_keys:
            if key not in adata.obs:
                missing_cont_keys.append(key)
        if len(missing_cont_keys) > 0:
            raise ValueError(f"The following continuous covariate keys are missing from adata.obs: {missing_cont_keys}")
        

    print("\nrun-scvi.py:\tSetting up scVI model...")
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts", # Ojo hard coding layer name for now
        categorical_covariate_keys=args.categorical_batch_keys,  # Ojo: can use batch_key here. Read into this.
        continuous_covariate_keys=args.continuous_covariate_keys)
    vae = scvi.model.SCVI(
        adata, 
        n_latent=args.n_latent)
    print(f"\nrun-scvi.py:\tTraining scVI with n_latent={args.n_latent}...")
    vae.train(
        batch_size=args.batch_size,
        max_epochs=args.max_epochs, 
        early_stopping=True, 
        early_stopping_patience=args.early_stopping_patience)
    adata.obsm[f"X_{args.cluster_name_prefix}"] = vae.get_latent_representation()

    # neighbors / leiden / umap on that embedding
    # 1) Build a *new* neighbor graph
    print("\nrun-scvi.py:\tComputing neighbors on scVI embedding...")
    sc.pp.neighbors(
        adata,
        use_rep=f"X_{args.cluster_name_prefix}",
        n_neighbors=args.n_neighbors,
        key_added=f"{args.cluster_name_prefix}_neighbors",
    )

    print("\nrun-scvi.py:\tRunning Leiden clustering...")
    # 2) Leiden clustering on *that* graph
    for res in args.resolutions:
        cluster_key = f"{args.cluster_name_prefix}_leiden_{str(res).replace('.', 'p')}"
        print(f"\nrun-scvi.py:\tClustering at resolution {res}, storing in {cluster_key}...")
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=cluster_key,
            neighbors_key=f"{args.cluster_name_prefix}_neighbors",
        )

        if args.find_markers:
            print(f"\nrun-scvi.py:\tFinding markers for clusters in {cluster_key}...")
            adata = find_markers(adata, cluster_key=cluster_key, layer = "norm_log")

    print("\nrun-scvi.py:\tRunning UMAP...")
    # 3) UMAP using *that* same graph
    sc.tl.umap(
        adata,
        neighbors_key=f"{args.cluster_name_prefix}_neighbors",
        key_added=f"{args.cluster_name_prefix}_umap",
    )

    adata.X = adata.layers["norm_log"].copy()  # set back to normalized data

    # save layer counts as sparse
    if args.convert_counts_to_dense:
        counts = adata.layers["counts"]

        # sanity check: values should still be integer-valued
        if not np.allclose(counts, np.round(counts)):
            raise ValueError(
                "Counts are no longer integer-valued; cannot safely cast back to int."
            )

        counts_int = np.rint(counts).astype(np.int32, copy=False)
        adata.layers["counts"] = sp.csr_matrix(counts_int)
    
    print(f"\nrun-scvi.py:\tSaving output AnnData to: {args.output_path}")
    adata.write_h5ad(args.output_path)

    print("\nrun-scvi.py:\tDone!")

if __name__ == "__main__":
    main()