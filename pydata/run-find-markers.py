import argparse
import scanpy as sc
from data_utils import add_back_og_counts, find_markers

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-path", "-i", required=True, help="Input .h5ad path")
    ap.add_argument(
        "--og-counts-path",
        type=str,
        default=None,
        required=False,
        help="Path to original counts .h5ad file to restore counts from. \
            If not set, uses counts from input .h5ad `.raw.X`.",
    )
    ap.add_argument(
        "--cluster-keys",
        type=str,
        required=True,
        nargs="+",
        default=["leiden_0p5"],
        help="set of cluster keys in adata.obs to find markers for. Can specify multiple keys.",
    )
    
    # print arguments for reference
    for arg in vars(ap.parse_args()):
        print(f"{arg}: {getattr(ap.parse_args(), arg)}")

    return ap.parse_args()

def main():

    args = parse_args()

    print(f"\nrun-find-markers.py:\tLoading input AnnData from: {args.input_path}")
    adata = sc.read_h5ad(args.input_path)
    print(f"\nrun-find-markers.py:\tNumber of cells: {adata.n_obs}, Number of genes: {adata.n_vars}")

    # if args.og_counts_path is not None:
    #     print(f"\nrun-find-markers.py:\tLoading ORIGINAL AnnData from {args.og_counts_path}...")
    #     og_adata = sc.read_h5ad(args.og_counts_path)
    #     print(f"\nrun-find-markers.py:\tRestoring original counts to adata object...")
    #     adata = add_back_og_counts(adata, og_adata)
    #     print(
    #         f"\nrun-find-markers.py:\tNumber of cells after restoring counts: {adata.n_obs}, Number of genes: {adata.n_vars}"
    #     )
    #     del og_adata  # free memory

    # # need X to be normalized + log1p transformed of raw counts 
    # print("\nrun-find-markers.py:\tResetting adata.X to normalized + log1p of raw counts for marker finding...")
    # adata.X = adata.raw.X.copy()
    # print("Normalizing total counts per cell...")
    # sc.pp.normalize_total(adata, target_sum=1e6)
    # print("Log1p transform...")
    # sc.pp.log1p(adata)

    for cluster_key in args.cluster_keys:
        print(f"\nrun-find-markers.py:\tFinding markers for clusters in {cluster_key}...")
        adata = find_markers(adata, cluster_key=cluster_key)

        # sc.tl.rank_genes_groups(
        #     adata,
        #     groupby=cluster_key,
        #     method="wilcoxon",
        #     corr_method="benjamini-hochberg",
        #     pts=True,
        #     use_raw=False,
        #     key_added=f"rank_genes_{cluster_key}"
        # )


    print(f"\nrun-find-markers.py:\tSaving output AnnData with markers to: {args.input_path}")
    adata.write_h5ad(args.input_path)
    print("Done.")

if __name__ == "__main__":
    main()