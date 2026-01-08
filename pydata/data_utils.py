import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
import scanpy as sc
import scanpy.external as sce
import mygene
mg = mygene.MyGeneInfo()


def add_back_og_counts(adata, og_adata):
    """
    given an AnnData object 'adata' that may have modified counts (e.g. normalized, log1p, etc),
    and an original AnnData object 'og_adata' with raw counts, restore the counts from og_adata into adata.
    Assumes that the cells in adata are a subset of those in og_adata, matched by obs_names.
    """

    # subset the og_adata to only include cells present in adata
    og_adata_subset = og_adata[adata.obs_names, :]

    new_adata = ad.AnnData(
        X=og_adata_subset.X,
        obs=adata.obs.copy(),
        var=og_adata_subset.var.copy(),
        uns=adata.uns.copy(),
        obsm=adata.obsm.copy(),
    )

    new_adata.raw = new_adata.copy()

    return new_adata

def ensg_to_symbol(ensg_list):
    """
    Convert a list of ENSG IDs -> gene symbols using mygene.info
    Returns a dict {ENSG: SYMBOL or None}
    """
    # mygene expects no version suffix, so strip ".1", ".2", etc
    clean_ids = [g.split(".")[0] for g in ensg_list]

    results = mg.querymany(
        clean_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=False
    )

    out = {}
    for q in results:
        ensg = q["query"]
        if "notfound" in q:
            out[ensg] = None
        else:
            out[ensg] = q.get("symbol", None)

    # return with original ENSG keys
    return {orig: out[orig.split(".")[0]] for orig in ensg_list}



def find_markers(adata, cluster_key, layer=None, **kwargs):
    """
    Find marker genes for each cluster in adata, given the cluster assignment in adata.obs
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        pts=True,
        use_raw=False,
        key_added=f"rank_genes_{cluster_key}",
        layer=layer,
        **kwargs
    )

    return adata

def filter_markers(adata, 
                   cluster_key, 
                   logfc_threshold=0, 
                   adj_pval_threshold=0.05, 
                   pct_nz_group_threshold=0, 
                   pct_nz_reference_threshold=1, 
                   pct_nz_diff_threshold=-np.inf,
                   top_n=None,
                   order_by_col="logfoldchanges",
                   filter_positive=True,
                   filter_out_ensg=True):
    """
    Filter marker genes found by sc.tl.rank_genes_groups based on logFC and adjusted p-value thresholds.
    Returns a dataframe of filtered markers with one row per gene per cluster.

    Parameters:
    - adata: AnnData object with marker results in adata.uns
    - cluster_key: str, the obs column used for clustering
    - logfc_threshold: float, minimum log fold change to consider a gene a marker
    - adj_pval_threshold: float, maximum adjusted p-value to consider a gene a marker
    - pct_nz_group_threshold: float, minimum percent of non-zero expression in the cluster
    - pct_nz_reference_threshold: float, maximum percent of non-zero expression in the non-cluster cells
    - pct_nz_diff_threshold: float, minimum difference in percent non-zero between cluster and reference
    - top_n: int or None, if set, only return the top N markers per cluster based on order_by
    - order_by_col: str, column name to order markers by when selecting top_n
    - filter_positive: bool, if True, only keep markers with positive logFC
    - filter_out_ensg: bool, if True, filter out genes whose gene names start with 'ENSG'

    Returns:
    - filtered_markers: pd.DataFrame of filtered marker genes
    """

    markers = sc.get.rank_genes_groups_df(
            adata,
            group=None,
            key=f"rank_genes_{cluster_key}",
        ).dropna(subset=["logfoldchanges"])

    filtered_markers = markers[
        (markers['logfoldchanges'] >= logfc_threshold) &
        (markers['pvals_adj'] <= adj_pval_threshold) & 
        (markers['pct_nz_group'] >= pct_nz_group_threshold) &
        (markers['pct_nz_reference'] <= pct_nz_reference_threshold) & 
        (markers['pct_nz_group'] - markers['pct_nz_reference'] >= pct_nz_diff_threshold) &
        ((markers['logfoldchanges'] > 0) if filter_positive else True)
    ]

    filtered_markers = (filtered_markers.merge(adata.var, left_on='names', right_index=True, how='left'))

    if filter_out_ensg:
        gn = filtered_markers["gene_name"].astype("string")
        filtered_markers = filtered_markers[~gn.fillna("").str.startswith("ENSG")].reset_index(drop=True)

    if top_n is not None:
        filtered_markers = (filtered_markers
            .sort_values(by=['group', order_by_col], ascending=[True, False])
            .groupby('group')
            .head(top_n)
            .reset_index(drop=True)
        )

    return filtered_markers[[
        "group", 'gene_name', "logfoldchanges", "pvals_adj", "pct_nz_group", "pct_nz_reference"
    ]]


def is_integer_valued(X, tol=1e-8, max_check=200_000):
    import numpy as np
    from scipy import sparse

    if sparse.issparse(X):
        data = X.data
    else:
        data = np.asarray(X).ravel()

    if data.size == 0:
        return True

    if data.size > max_check:
        data = np.random.default_rng(0).choice(data, max_check, replace=False)

    if np.issubdtype(data.dtype, np.integer):
        return True

    return bool(np.all(np.abs(data - np.round(data)) < tol))


def merge_h5ads(list_of_paths):
    """
    Merge multiple AnnData .h5ad files with:
      - union of obs metadata columns (missing values -> NA)
      - union of genes
      - ensure unique cell IDs
    """
    adatas = []
    for i, p in enumerate(list_of_paths):
        if (i % 10 == 0) and i > 0:
            print(f"  loaded {i} / {len(list_of_paths)} files...")
        a = sc.read_h5ad(p)
        a.obs_names = a.obs.library__barcode
        adatas.append(a)

    print (f"  loaded all {len(list_of_paths)} files, now concatenating...")
    merged = ad.concat(
        adatas,
        join="outer",       # union of genes + obs columns
        merge="same",       # for var columns with same name, ensure values are identical 
        axis=0
    )

    return merged

def neighbor_cluster_umap(
        adata_sub,
        adata_full,
        resolutions,
        n_neighbors=30,
        n_pcs=30,
        cluster_name_prefix="leiden",
        dim_reduction_key="pca",
        umap_key="umap"
):  
    
    print(f"\ndata_utils.py:\tComputing Neighbors on {dim_reduction_key} space...")
    sc.pp.neighbors(
        adata_sub,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep = f"X_{dim_reduction_key}")
    adata_full.obsp[f"connectivities_{dim_reduction_key}"] = adata_sub.obsp[f"connectivities"]
    adata_full.obsp[f"distances_{dim_reduction_key}"] = adata_sub.obsp[f"distances"]
    adata_full.uns[f"neighbors_{dim_reduction_key}"] = adata_sub.uns.get(f"neighbors", {})

    print(f"\ndata_utils.py:\tComputing UMAP on HVG {dim_reduction_key} space...")
    sc.tl.umap(adata_sub)
    adata_full.obsm[f"X_{umap_key}"] = adata_sub.obsm[f"X_umap"]
    adata_full.uns[umap_key] = adata_sub.uns.get("umap", {})

    for res in resolutions:
        key = f"{cluster_name_prefix}_{str(res).replace('.', 'p')}"
        print(f"\ndata_utils.py:\tClustering with Leiden, resolution = {res}\n\tkey = {key}")
        sc.tl.leiden(adata_sub, resolution=res, key_added=key)
        adata_full.obs[key] = adata_sub.obs[key]
    
    return adata_full
    

def normalize_scale_pca_cluster_umap(
    adata,
    n_hvg=2500,
    n_pcs=30,
    n_neighbors=20,
    resolutions=(0.2, 0.5, 1.0),
    target_sum_norm = 1e6,
    max_value_scale = 10,
    cluster_name_prefix = "leiden",
    run_harmony_on = None,          # e.g. ['donor_id', ...]
    run_find_markers=False,
    counts_layer: str = "counts",
    lognorm_layer: str = "norm_log",
    **kwargs
):
    """
    Normalize, log-transform, select HVGs, scale, PCA, neighbors, Leiden, UMAP.
    """

    adata.raw = None # clear raw to save memory

     # (1) Ensure counts layer exists and is integer-like
    if counts_layer in adata.layers:
        if is_integer_valued(adata.layers[counts_layer]):
            adata.X = adata.layers[counts_layer].copy()
        else:
            raise ValueError(
                f"\ndata_utils.py:\tLayer '{counts_layer}' exists but contains non-integer values! "
                "Provide raw counts in a layer, then rerun."
            )
    elif counts_layer not in adata.layers:
        if is_integer_valued(adata.X):
            adata.layers[counts_layer] = adata.X.copy()
        else:
            raise ValueError(
                f"\ndata_utils.py:\tLayer '{counts_layer}' is missing and adata.X contains non-integer values. "
                "Provide raw counts in adata.X or in a layer, then rerun."
            )
    
    print(f"\ndata_utils.py:\tUsing adata.X for normalization and downstream analysis. Example values:\n")
    num_genes = min(2000, adata.X.shape[1])
    min_num_genes = max(0, num_genes - 10)
    ex = adata.X[1:10, min_num_genes:num_genes]
    print(ex.toarray() if hasattr(ex, "toarray") else ex)

    # 2) Normalize total
    print(f"\ndata_utils.py:\tNormalizing total counts per cell...")
    sc.pp.normalize_total(adata, target_sum=target_sum_norm)

    # 3) log1p
    print(f"\ndata_utils.py:\tLog1p transform...")
    sc.pp.log1p(adata)

    if lognorm_layer is not None:
        adata.layers[lognorm_layer] = adata.X.copy()

    # 4) HVGs (do NOT subset by default, to keep all genes)
    print(f"\ndata_utils.py:\tSelecting highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg,
        flavor="seurat_v3",
        subset=False,
        layer=counts_layer
    )

    # Scale, PCA, UMAP, Leiden on HVG subset only
    hvgs = adata.var["highly_variable"].to_numpy()
    if hvgs.sum() == 0:
        raise ValueError("\ndata_utils.py:\tNo HVGs selected (adata.var['highly_variable'] has 0 True entries).")

    # Copy only HVG slice (much smaller than copying all genes)
    adata_hvg = adata[:, hvgs].copy()

    # 5) Scale
    print(f"\ndata_utils.py:\tScaling data on HVGs...")
    sc.pp.scale(adata_hvg, max_value=max_value_scale)

    # 6) PCA
    print(f"\ndata_utils.py:\tRunning PCA on HVGs...")
    sc.tl.pca(adata_hvg, svd_solver="arpack")

    # Copy PCA embeddings back to original adata
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns.get("pca", {})

    #7) Neighbors, Leiden, UMAP
    adata = neighbor_cluster_umap(
        adata_sub = adata_hvg,
        adata_full=adata,
        resolutions=resolutions,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        cluster_name_prefix=cluster_name_prefix,
        dim_reduction_key="pca",
        umap_key="umap"
    )

    # 8) optionally run harmony
    if run_harmony_on is not None:
        print(f"\ndata_utils.py:\tRunning Harmony integration...")
        sce.pp.harmony_integrate(adata_hvg, key=run_harmony_on)
        adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]
        adata.uns["pca_harmony"] = adata_hvg.uns.get("pca_harmony", {})
        adata = neighbor_cluster_umap(
            adata_sub = adata_hvg,
            adata_full=adata,
            resolutions=resolutions,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            cluster_name_prefix=f"harmony_{cluster_name_prefix}",
            dim_reduction_key="pca_harmony",
            umap_key="umap_harmony"
        )   
        
    # 9) optionally find markers
    if run_find_markers:
        for res in resolutions:
            cluster_key = f"{cluster_name_prefix}_{str(res).replace('.', 'p')}"
            print(f"\ndata_utils.py:\tFinding markers for clusters in {cluster_key}...")
            adata = find_markers(adata, cluster_key=cluster_key, layer = lognorm_layer, **kwargs)

            if run_harmony_on is not None:
                harmony_cluster_key = f"harmony_{cluster_key}"
                print(f"\ndata_utils.py:\tFinding markers for clusters in {harmony_cluster_key}...")
                adata = find_markers(adata, cluster_key=harmony_cluster_key, layer = lognorm_layer, **kwargs)

    return adata