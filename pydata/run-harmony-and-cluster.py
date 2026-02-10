import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
from pathlib import Path
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scvi

base="/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets"


# import functions from ~/sc-online/pydata/data_utils.py
import sys
from pathlib import Path
sys.path.append(str(Path.home() / "sc-online" / "pydata"))
from data_utils import neighbor_cluster_umap, find_markers

from plot import feature_plot

from importlib import reload as reload_module
import data_utils
import plot
reload_module(data_utils)
reload_module(plot)
import scanpy.external as sce

print("reading inh")
nurr_inh = sc.read_h5ad(f"{base}/sn-vta-merged-04-neuron-inh-nurr.h5ad")
print(f"nurr_inh shape: {nurr_inh.shape}")
print("reading exc")
nurr_exc = sc.read_h5ad(f"{base}/sn-vta-merged-04-neuron-exc-nurr.h5ad")
print(f"nurr_exc shape: {nurr_exc.shape}")

for k in ["X_cnmf_norm_clean_harmony", "X_harmony_cnmf_norm_clean", "umap_harmony", ]:
    if k in nurr_inh.obsm.keys():   
        print(f"deleting {k} from nurr_inh")
        del nurr_inh.obsm[k]
    if k in nurr_exc.obsm.keys():
        print(f"deleting {k} from nurr_exc")
        del nurr_exc.obsm[k]


print("running harmony on inh")
sce.pp.harmony_integrate(
    nurr_inh, 
    key="donor_id", 
    basis='X_nurr_inh_04_cnmf', 
    adjusted_basis='X_harmony_cnmf_norm_clean', 
    max_iter_harmony=20)

nurr_inh.write_h5ad(f"{base}/sn-vta-merged-04-neuron-inh-nurr.h5ad")
print(f"nurr_inh shape after harmony: {nurr_inh.shape}")
print(f"nurr_inh obsm keys after harmony: {nurr_inh.obsm.keys()}")
print(f"nurr_inh obsm shape after harmony: {nurr_inh.obsm['X_harmony_cnmf_norm_clean'].shape}")
print(f"number of NaNs in harmony result: {np.sum(np.any(np.isnan(nurr_inh.obsm['X_harmony_cnmf_norm_clean']), axis=1 ))}")
print("writing inh")


print("running harmony on exc")
sce.pp.harmony_integrate(
    nurr_exc, 
    key="donor_id", 
    basis='X_cnmf_norm_clean', 
    adjusted_basis='X_harmony_cnmf_norm_clean', 
    max_iter_harmony=20)
print(f"nurr_exc shape after harmony: {nurr_exc.shape}")
print(f"nurr_exc obsm keys after harmony: {nurr_exc.obsm.keys()}")
print(f"nurr_exc obsm shape after harmony: {nurr_exc.obsm['X_harmony_cnmf_norm_clean'].shape}")
print(f"number of NaNs in harmony result: {np.sum(np.any(np.isnan(nurr_exc.obsm['X_harmony_cnmf_norm_clean']), axis=1 ))}")
print("writing exc")    
nurr_exc.write_h5ad(f"{base}/sn-vta-merged-04-neuron-exc-nurr.h5ad")

print("clustering and UMAP for exc")
nurr_exc = neighbor_cluster_umap(
    adata_sub = nurr_exc,
    adata_full=nurr_exc,
    resolutions=[0.2, 0.5, 0.8],
    n_neighbors=30,
    n_pcs=nurr_exc.obsm["X_harmony_cnmf_norm_clean"].shape[1],
    cluster_name_prefix=f"harmony_cnmf_norm_clean_leiden",
    dim_reduction_key="harmony_cnmf_norm_clean",
    umap_key="umap_harmony_cnmf_norm_clean",
    run_find_markers=True
)
print("writing exc")
nurr_exc.write_h5ad(f"{base}/sn-vta-merged-04-neuron-exc-nurr.h5ad")

print("clustering and UMAP for inh")
nurr_inh = neighbor_cluster_umap(
    adata_sub = nurr_inh,
    adata_full=nurr_inh,
    resolutions=[0.2, 0.5, 0.8],
    n_neighbors=30,
    n_pcs=nurr_inh.obsm["X_harmony_cnmf_norm_clean"].shape[1],
    cluster_name_prefix=f"harmony_cnmf_norm_clean_leiden",
    dim_reduction_key="harmony_cnmf_norm_clean",
    umap_key="umap_harmony_cnmf_norm_clean",
    run_find_markers=True
)
print("writing inh")

nurr_inh.write_h5ad(f"{base}/sn-vta-merged-04-neuron-inh-nurr.h5ad")


