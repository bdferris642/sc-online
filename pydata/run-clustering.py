import argparse
import scanpy as sc
from data_utils import add_back_og_counts, normalize_scale_pca_cluster_umap

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-path", "-i", required=True, help="Input .h5ad path")
    ap.add_argument(
        "--og-counts-path",
        type=str,
        default=None,
        help="Path to original counts .h5ad file to restore counts from. \
            If not set, uses counts from input .h5ad.",
    )
    # if output-patch is not provided, set to input-path
    ap.add_argument("--output-path", "-o", 
                    required=False, 
                    default=None,
                    help="Output .h5ad path. If not set, overwrites input path.")
    ap.add_argument(
        "--resolutions", "-r", nargs="+", type=float, default=[0.2, 0.5, 0.8]
    )
    ap.add_argument("--n-hvg", type=int, default=2500)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-neighbors", type=int, default=20)
    ap.add_argument("--cluster-name-prefix", type=str, default="leiden_")
    ap.add_argument("--vars-to-regress", nargs="+", default=None,
                    help="List of obs columns to regress out after normalization.")
    ap.add_argument(
        "--run-harmony-on",
        nargs="+",
        default=None,
        help="List of obs columns to run Harmony integration on (e.g. batch effects). \
            If not set, no Harmony is run.",
    )
    ap.add_argument(
        "--find-markers",
        action="store_true",
        default=False,
        help="Whether to find marker genes for each cluster after clustering."
    )

    args = ap.parse_args()
    if args.output_path is None:
        args.output_path = args.input_path

    # print arguments for reference
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")
    
    return args


def main():
    args = parse_args()

    print(f"\nrun-clustering.py:\tLoading input AnnData from: {args.input_path}")
    adata = sc.read_h5ad(args.input_path)
    print(f"\nrun-clustering.py:\tNumber of cells: {adata.n_obs}, Number of genes: {adata.n_vars}")

    if args.og_counts_path is not None:

        print(f"\nrun-clustering.py:\tLoading ORIGINAL AnnData from {args.og_counts_path}...")
        og_adata = sc.read_h5ad(args.og_counts_path)
        print(f"\nrun-clustering.py:\tRestoring original counts to adata object...")
        adata = add_back_og_counts(adata, og_adata)
        print(
            f"\nrun-clustering.py:\tNumber of cells after restoring counts: {adata.n_obs}, Number of genes: {adata.n_vars}"
        )
        del og_adata  # free memory

    print("\nrun-clustering.py:\tCalling Clustering Routine...")
    adata = normalize_scale_pca_cluster_umap(
        adata,
        n_hvg=args.n_hvg,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolutions=args.resolutions,
        target_sum_norm=1e6,
        max_value_scale=10,
        cluster_name_prefix=args.cluster_name_prefix,
        vars_to_regress=args.vars_to_regress,
        run_harmony_on=args.run_harmony_on,
        run_find_markers=args.find_markers
    )

    print(f"\nrun-clustering.py:\tSaving output to: {args.output_path}")
    adata.write_h5ad(args.output_path)

    print("\nrun-clustering.py:\tDone.")

if __name__ == "__main__":
    main()
