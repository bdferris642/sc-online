import argparse 
import decoupler as dc
import numpy as np
import os
import pandas as pd
import scanpy as sc
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input-dir",
    type=str,
    help="Path to the directory containing the input h5ad files")
parser.add_argument("--output-dir",
    type=str,
    default="pseudobulk",
    help="Directory in which to save the output files")
parser.add_argument("--min-num-cells",
    type=int,
    default=10,
    help="Minimum number of cells per sample")
parser.add_argument("--gene-log-expr-threshold",
    type=float,
    default=1.5,
    help="Minimum mean log-transformed expression of genes to keep")
parser.add_argument("--ct-id",
    type=str,
    default="cell_class",
    help="Cell type identifier")
parser.add_argument("--sample-id",
    type=str,
    default="participant_id",
    help="Sample identifier")
parser.add_argument("--strs-to-skip",
    type=str,
    default="endo.h5ad,opc.h5ad",
    help="String to skip in the file name")

args = parser.parse_args()
CT_ID = args.ct_id
SAMPLE_ID = args.sample_id
MIN_NUM_CELLS = args.min_num_cells
GENE_LOG_EXPR_THRESHOLD = args.gene_log_expr_threshold
BASE = args.input_dir
OUTPUT_DIR = args.output_dir
STRS_TO_SKIP = args.strs_to_skip.split(",")

# CONSTNATS
sc.settings.n_jobs = 30
MIN_COUNTS = 100
MAX_SCALE_EXPR_VALUE = 10


os.makedirs(OUTPUT_DIR, exist_ok=True)

input_files = [f"{BASE}/{f}" for f in os.listdir(BASE) if f.endswith(".h5ad")]
print(input_files)

for file in input_files:
    print(f"Processing {file}...")

    for s in STRS_TO_SKIP:
        if s in file:
            print(f"Skipping {s} file")
            continue

    # Load the h5ad file
    adata = sc.read_h5ad(file)
    if adata.raw is not None:
        print("adata.raw exists! Using adata.raw instead.")
        adata.X = adata.raw.X  # Use raw counts if available

    adata.X = adata.X.toarray() # convert from sparse to dense 

    # count cells per SAMPLE_ID
    sample_counts = adata.obs[SAMPLE_ID].value_counts()

    # get SAMPLE_IDs with at least MIN_NUM_CELLS and subset
    sample_ids = sample_counts[sample_counts >= MIN_NUM_CELLS].index

    print(f"removing {(adata.X.sum(axis=0) == 0).sum()} genes with zero counts") 
    print(f" removing {(adata.X.sum(axis=1) == 0).sum()} cells have zero counts")

    adata = adata[:, adata.X.sum(axis=0) > 0]  # Remove genes with 0 expression
    adata = adata[adata.X.sum(axis=1) > 0, :]  # Remove cells with 0 expression

    num_nas_sample = adata.obs[SAMPLE_ID].isna().sum()
    num_nas_ct = adata.obs[CT_ID].isna().sum()  # Count NaNs in CT_ID
    print(f"removing {num_nas_sample} samples with missing sample ID and {num_nas_ct} cells with missing cell type ID")
    adata = adata[~adata.obs[SAMPLE_ID].isna() & ~adata.obs[CT_ID].isna(), :]

    print(f"adata shape: {adata.shape}")
    
    data_scanpy_1 = dat.from_scanpy(
        adata,
        cell_type_identifier=CT_ID,
        sample_identifier=SAMPLE_ID
    )
    print(f"data_scanpy_1 shape: {data_scanpy_1.shape}")

    
    combined_df_corrected = pd.concat(
        [data_scanpy_1.obs, 
        pd.DataFrame(
            data_scanpy_1.X, 
            index=data_scanpy_1.obs.index, 
            columns=data_scanpy_1.var.index)], axis=1)
    print(f"combined_df_corrected shape: {combined_df_corrected.shape}")

    # if adding `_{CT_ID}` to index names
    # combined_df_corrected.index = combined_df_corrected.index + "_" + combined_df_corrected.columns[0]
    
    #new_column_names = ['cluster' + str(col) for col in data_scanpy_1.var.index]
    # Update the relevant columns in combined_df_corrected with the new names
    #combined_df_corrected.columns = list(data_scanpy_1.obs.columns) + new_column_names
    
    # Pseudobulk aggregation
    pdata = dc.get_pseudobulk(
        adata,
        sample_col=SAMPLE_ID, 
        groups_col=CT_ID,
        skip_checks=True,
        mode="sum",
        min_cells=MIN_NUM_CELLS,
        min_counts=MIN_COUNTS,
        remove_empty=True
    )
    print(f"pseudobulk (pdata) shape: {pdata.shape}")
    
    # Normalize total counts and log-transform
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)
    
    # Filter genes with mean log-transformed expression >= GENE_LOG_EXPR_THRESHOLD
    gene_filter = pdata.X.mean(axis=0) >= GENE_LOG_EXPR_THRESHOLD
    pdata = pdata[:, gene_filter]
    
    # Scale the data to a maximum value
    sc.pp.scale(pdata, max_value=MAX_SCALE_EXPR_VALUE)
    
    # Save the processed data as a CSV file
    data=pd.DataFrame(pdata.X, index=pdata.obs_names, columns=pdata.var_names)
    print(f"cleaned pseudobulk (data) shape: {data.shape}")

    output_file1 = os.path.join(OUTPUT_DIR, f"{os.path.splitext(os.path.basename(file))[0]}_expression_matrix_ds.csv")
    output_file2 = os.path.join(OUTPUT_DIR, f"{os.path.splitext(os.path.basename(file))[0]}_composition_matrix_ds.csv")

    # remove cell type suffix from index
    data.index = ["".join(i.split("_")[:-1]) for i in data.index]

    common_samples = data.index.intersection(combined_df_corrected.index)
    print(f"Number of common samples: {len(common_samples)}")

    # Filter both datasets to only include common samples
    data_aligned = data.loc[common_samples]
    combined_df_corrected_aligned = combined_df_corrected.loc[common_samples]
    print(f"Aligned data shape: {data_aligned.shape}")
    print(f"Aligned combined_df_corrected shape: {combined_df_corrected_aligned.shape}")
    
    data_aligned.to_csv(output_file1)
    combined_df_corrected_aligned.to_csv(output_file2)
    
    print(f"Saved {output_file1}")
    print(f"Saved {output_file2}")
print("Processing completed.")