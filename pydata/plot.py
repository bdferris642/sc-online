import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scanpy as sc

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
    for i, gene in enumerate(features):
        if gene not in adata.var['gene_name'].values:
            print(f"{gene} not found in dataset, skipping.")
            axes[i].axis('off')
            continue
        sc.pl.embedding(adata, color=gene, gene_symbols="gene_name", 
                        basis=basis, ax=axes[i], layer=layer,
                        title=f'{gene} expression', 
                        show=False, cmap='viridis', vmax='p99', use_raw=False)
    # hide any unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    plt.tight_layout()
    plt.show()