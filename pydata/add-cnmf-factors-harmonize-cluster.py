# Adds cNMF factors from a df to ana adata (skipping geps if specified), normalizes to total usage per cell, 
# then runs harmony and clustering. 
# in the above step, cells are removed if they have no usages applied (e.g. due to low expression)
# or if they have all zero usages after skipping specified geps.

import os
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import re
import scanpy.external as sce
from data_utils import apply_cnmf_usage_df, neighbor_cluster_umap, is_integer_valued


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad-path", "-p", required=True, help="Input .h5ad path")
    ap.add_argument("--usages-path", "-u", required=True, help="Input CNMF usages .csv path")
    ap.add_argument("--output-path", "-o", required=False, help="Output .h5ad path", default=None)
    ap.add_argument("--geps-to-skip", "-g", nargs="+", help="List of GEPs to skip when normalizing usages \
                    \n(e.g. if you know some GEPs are technical artifacts and want to exclude them from the normalized usage calculation)",
                    default=None)
    ap.add_argument("--dimension-reduction-name", "-d", default="cnmf_norm_clean", 
                    help="Name for the dimension reduction to store in adata.obsm (default: cnmf_norm_clean)")
    ap.add_argument("--batch-key", "-b", default="donor_id", help="Key in adata.obs for donor/batch IDs (default: donor_id)")
    ap.add_argument("--max-iter-harmony", "-m", type=int, default=20, help="Maximum iterations for Harmony (default: 20)")
    ap.add_argument("--resolutions", "-r", nargs="+", type=float, default=[0.1, 0.2, 0.3, 0.4, 0.5], 
                    help="List of resolutions for Leiden clustering (default: [0.1, 0.2, 0.3, 0.4, 0.5])")
    ap.add_argument("--num-neighbors", "-n", type=int, default=30, help="Number of neighbors for clustering/UMAP (default: 30)")
    ap.add_argument("--prune-cnmf-cols", action="store_true", 
                    help="Prune cNMF columns before applying usages (default: False). \
                        \nThis can be helpful if your adata already has cNMF columns from a previous run \
                        \nand you want to make sure they don't interfere with the new ones being added.")
    ap.add_argument("--force", action="store_true", 
                    help="Force the script to run even if more than 5 percent of cells are removed after applying cNMF usages. \
                    \nUse with caution and make sure to check the output adata carefully if you use this flag.")
    ap.add_argument("--run-find-markers", action="store_true",
                    help="Whether to run sc.tl.rank_genes_groups to find marker genes for each cluster after clustering. \
                    \nThis can add significant runtime, so only use if you need marker genes for the new clusters.")
    ap.add_argument("--distance-metric", default="cosine", 
                    help="Distance metric for neighbor graph construction in scanpy.pp.neighbors (default: cosine)")
    
    args = ap.parse_args()
    # print arguments for reference
    if args.output_path is None:
        args.output_path = args.h5ad_path

    for k, v in vars(args).items():
        print(f"{k}: {v}")
    return args

def main():
    args = parse_args()
    
    if args.geps_to_skip is not None:
        args.geps_to_skip = [int(gep) for gep in args.geps_to_skip]
    
    if args.resolutions is not None:
        args.resolutions = [float(res) for res in args.resolutions]

    # check that usages exists 
    if not os.path.exists(args.usages_path):
        raise FileNotFoundError(f"CNMF usages file not found: {args.usages_path}")
    
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tLoading input AnnData from: {args.h5ad_path}")
    adata = sc.read_h5ad(args.h5ad_path)

    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tadata shape: {adata.shape}, \
          \n\tadata.layers: {list(adata.layers.keys())}, \
          \n\tadata.obsm keys: {list(adata.obsm.keys())}, adata.obs columns: {list(adata.obs.columns)}")
    
    if not ("counts" in adata.layers):
        if is_integer_valued(adata.X):
            print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tadata.X appears to be raw counts but no 'counts' layer found. \
                    \nCopying adata.X to adata.layers['counts']")
            adata.layers["counts"] = adata.X.copy()
        else:
            raise ValueError(f"adata does not have a 'counts' layer and adata.X does not appear to be raw counts. \
                             \nPlease provide raw counts in adata.layers['counts'] or ensure adata.X contains raw counts.")

    if ("counts" in adata.layers) and not ("norm_log" in adata.layers):
        print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tadata has raw counts in adata.layers['counts'] \
                \nbut no normalized log data in adata.layers['norm_log']. \
                \nNormalizing and log-transforming counts to adata.layers['norm_log']")
        adata.X = adata.layers["counts"].copy()  # set adata.X to raw counts for normalization
        sc.pp.normalize_total(adata, target_sum=1e6)
        sc.pp.log1p(adata)
        adata.layers["norm_log"] = adata.X.copy()  # store the normalized log data in a new layer
        
    orig_n_obs = adata.n_obs
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tOriginal adata shape: {adata.shape}")
    
    for k in list(adata.obsm.keys()):
        if "cnmf" in k.lower():
            print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tDeleting {k} from adata.obsm")
            del adata.obsm[k]  

    if args.prune_cnmf_cols:
        cnmf_cols = [col for col in adata.obs.columns if re.match(r"cnmf.*_k_[0-9]+_gep.*", col)]
        print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tPruning existing cNMF columns from adata.obs: {cnmf_cols}")
        adata.obs.drop(columns=cnmf_cols, inplace=True)  
    
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tApplying cNMF usages from: {args.usages_path}")
    adata = apply_cnmf_usage_df(adata, cnmf_path=args.usages_path, geps_to_skip=args.geps_to_skip)
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tadata shape after applying cNMF usages: {adata.shape}")

    if adata.n_obs / orig_n_obs < 0.95 and not args.force:
        raise ValueError(f"{100 * (1 - adata.n_obs / orig_n_obs):.2f}% of cells were removed after applying cNMF usages. \
                         \n\tOriginal n_obs: {orig_n_obs}, new n_obs: {adata.n_obs}. \
                         \n\tPlease check your CNMF usages and geps_to_skip parameters as well as adata.obs_names!")
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tNumber of cells removed after applying cNMF usages: {orig_n_obs - adata.n_obs}")


    # add newly added obs cols as an obsm matrix
    usage_cols = sorted(
        [col for col in adata.obs.columns if re.match(r"cnmf.*_k_[0-9]+_gep.*_norm_clean$", col)])
    X = adata.obs[usage_cols].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float) # make sure numeric
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\
          \tAdding normalized CNMF usage columns to adata.obsm['X_{args.dimension_reduction_name}']")
    adata.obsm[f"X_{args.dimension_reduction_name}"] = X

    # first cluster on the native cNMF factors
    adata = neighbor_cluster_umap(
        adata_sub = adata,
        adata_full=adata,
        resolutions=args.resolutions,
        n_neighbors=args.num_neighbors,
        n_pcs=adata.obsm[f"X_{args.dimension_reduction_name}"].shape[1],
        cluster_name_prefix=f"{args.dimension_reduction_name}_leiden",
        dim_reduction_key=f"{args.dimension_reduction_name}",
        umap_key=f"umap_{args.dimension_reduction_name}",
        metric = args.distance_metric,
        run_find_markers=args.run_find_markers
    )
    
    # then run harmony on the cNMF factors, cluster, and find markers on those as well
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tRunning harmony integration on {args.dimension_reduction_name}")
    sce.pp.harmony_integrate(
        adata, 
        key=args.batch_key, 
        basis=f'X_{args.dimension_reduction_name}', 
        adjusted_basis=f'X_harmony_{args.dimension_reduction_name}', 
        max_iter_harmony=args.max_iter_harmony)

    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tClustering and UMAP on Harmonized {args.dimension_reduction_name} factors")
    adata = neighbor_cluster_umap(
        adata_sub = adata,
        adata_full=adata,
        resolutions=args.resolutions,
        n_neighbors=args.num_neighbors,
        n_pcs=adata.obsm[f"X_harmony_{args.dimension_reduction_name}"].shape[1],
        cluster_name_prefix=f"harmony_{args.dimension_reduction_name}_leiden",
        dim_reduction_key=f"harmony_{args.dimension_reduction_name}",
        umap_key=f"umap_harmony_{args.dimension_reduction_name}",
        metric = args.distance_metric,
        run_find_markers=args.run_find_markers
    )
    
    print(f"\nadd-cnmf-factors-harmonize-cluster.py:\tWriting output AnnData to: {args.output_path}")
    adata.write_h5ad(args.output_path)


if __name__ == "__main__":
    main()