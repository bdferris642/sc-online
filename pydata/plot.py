import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scanpy as sc
from data_utils import get_gini_impurity


def plot_purity_by_cluster(adata, cluster_key, batch_key="donor_id"):
    # for every cluster, print the gini impurity / expected impurity 

    purities = {}
    gini_total = get_gini_impurity(adata.obs[batch_key].value_counts())
    print(f"Overall Gini Impurity across all cells: {gini_total:.4f}")

    for cluster in sorted(adata.obs[cluster_key].unique()):
        cluster_df = adata.obs[adata.obs[cluster_key] == cluster]
        donor_counts = cluster_df[batch_key].value_counts()
        gini_impurity = get_gini_impurity(donor_counts)
        purities[cluster] = 1 - gini_impurity


    # make a bar plot of the purities with order based on order of clusters in the umap dendrogram
    sc.tl.dendrogram(adata, groupby=cluster_key)
    dendro_order = adata.uns[f'dendrogram_{cluster_key}']['categories_ordered']
    purity_values = [purities[cluster] for cluster in dendro_order]
    plt.figure(figsize=(8,4))
    plt.bar(range(len(dendro_order)), purity_values)
    plt.xticks(range(len(dendro_order)), dendro_order, rotation=45, ha='right')
    plt.ylabel("1 - Gini Impurity")
    plt.title(f"{cluster_key} Purity by {batch_key}")

    # red dotted line for overall purity
    plt.axhline(1 - gini_total, color='red', linestyle='dotted', label='Overall Purity')
    plt.legend()
    plt.ylim((0, 1))
    
    plt.tight_layout()
    plt.show()



def feature_plot(
        adata, features, basis='umap_harmony', ncols=4, layer="norm_log"        
):
    """
    Plot feature expression on embedding.
    Assumes that adata.var contains a 'gene_name' column mapping var_names to gene names.
    
    :param adata: an anndata object
    :param features: list of gene names to plot
    :param basis: embedding basis to use for plotting
    :param ncols: number of columns in the plot grid
    """ 

    n_rows = int(np.ceil(len(features) / ncols))
    fig, axes = plt.subplots(n_rows, ncols, figsize=(ncols * 3.5, n_rows * 3))
    axes = axes.flatten()
    for i, feat in enumerate(features):
        if feat in adata.var['gene_name'].values:
            sc.pl.embedding(adata, color=feat, gene_symbols="gene_name", 
                            basis=basis, ax=axes[i], layer=layer,
                            title=f'{feat}', 
                            show=False, cmap='viridis', vmax='p99', use_raw=False)
        elif feat in adata.var_names:
            sc.pl.embedding(adata, color=feat, 
                            basis=basis, ax=axes[i], layer=layer,
                            title=f'{feat}', 
                            show=False, cmap='viridis', vmax='p99', use_raw=False)
        elif feat in adata.obs.columns:
            sc.pl.embedding(adata, color=feat, 
                            basis=basis, ax=axes[i],
                            title=f'{feat}', 
                            show=False, cmap='viridis', vmax='p99', use_raw=False)
        else:
            print(f"{feat} not found in var.gene_name, var_names, or obs.columns. Skipping.")
            axes[i].axis('off')
            continue
    # hide any unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    plt.tight_layout()
    plt.show()