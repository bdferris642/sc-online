# 01-make-eqtl-pseudobulk.py — Aggregate raw counts per sample × cell type, normalize, and filter genes.
#
# For each .h5ad file listed in --input-files:
#   1. Loads counts layer as X; validates required obs columns exist
#   2. Removes genes/cells with zero counts; removes NaN sample/cell-type labels
#   3. Casts --sample-id to str and strips whitespace
#   4. If --id-map provided: restricts to mapped participants and remaps IDs
#   5. Validates that each h5ad contains exactly one cell type (pipeline assumption)
#   6. Pseudobulk sum aggregation with decoupler.get_pseudobulk
#   7. CPM normalization (target_sum=1e6) + log1p
#   8. Gene filter: mean log-expression >= --gene-log-expr-threshold
#   9. Scale to max value MAX_SCALE_EXPR_VALUE
#
# Outputs per h5ad: named by cell_class from adata.obs[CT_ID], not the h5ad filename stem.
#   {output-dir}/{cell_class}_expression_matrix_ds.csv
#   {output-dir}/{cell_class}_composition_matrix_ds.csv

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
parser.add_argument("--input-files",
    type=str,
    help="Path to a newline-delimited text file listing h5ad files to process")
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
    default="",
    help="String to skip in the file name")
parser.add_argument("--id-map",
    type=str,
    default=None,
    help="Path to CSV with columns nucseq_participant_id and vcf_sample_id. "
         "Participant IDs are remapped and the object is restricted to participants "
         "present in this file. If omitted, no remapping is performed.")

args = parser.parse_args()
CT_ID = args.ct_id
SAMPLE_ID = args.sample_id
MIN_NUM_CELLS = args.min_num_cells
GENE_LOG_EXPR_THRESHOLD = args.gene_log_expr_threshold
INPUT_FILES_TXT = args.input_files
OUTPUT_DIR = args.output_dir
STRS_TO_SKIP = args.strs_to_skip.split(",")
ID_MAP = args.id_map

# CONSTANTS
sc.settings.n_jobs = 30
MIN_COUNTS = 100        # Minimum total counts per pseudobulk sample (passed to decoupler)
MAX_SCALE_EXPR_VALUE = 10  # Maximum value after scaling (sc.pp.scale max_value)


os.makedirs(OUTPUT_DIR, exist_ok=True)

with open(INPUT_FILES_TXT) as fh:
    input_files = [line.strip() for line in fh if line.strip()]
if not input_files:
    raise ValueError(f"No file paths found in {INPUT_FILES_TXT}.")
print(input_files)

for file in input_files:
    print(f"Processing {file}...")

    if not (STRS_TO_SKIP == [""] or STRS_TO_SKIP == []):
        if any(s in file for s in STRS_TO_SKIP):
            print(f"Skipping {file} (matched strs-to-skip)")
            continue

    # Load the h5ad file
    adata = sc.read_h5ad(file)
    if "counts" not in adata.layers:
        raise ValueError(f"'counts' layer not found in {file}. "
                         f"Available layers: {list(adata.layers.keys())}")
    adata.X = adata.layers["counts"]

    # Validate required obs columns exist before any downstream use
    for col in [SAMPLE_ID, CT_ID]:
        if col not in adata.obs.columns:
            raise ValueError(f"{col!r} not found in adata.obs. "
                             f"Available columns: {list(adata.obs.columns)}")

    adata.X = adata.X.toarray() # convert from sparse to dense

    print(f"removing {(adata.X.sum(axis=0) == 0).sum()} genes with zero counts")
    print(f" removing {(adata.X.sum(axis=1) == 0).sum()} cells have zero counts")

    adata = adata[:, adata.X.sum(axis=0) > 0]  # Remove genes with 0 expression
    adata = adata[adata.X.sum(axis=1) > 0, :]  # Remove cells with 0 expression

    num_nas_sample = adata.obs[SAMPLE_ID].isna().sum()
    num_nas_ct = adata.obs[CT_ID].isna().sum()  # Count NaNs in CT_ID
    print(f"removing {num_nas_sample} samples with missing sample ID and {num_nas_ct} cells with missing cell type ID")
    adata = adata[~adata.obs[SAMPLE_ID].isna() & ~adata.obs[CT_ID].isna(), :]

    # Cast SAMPLE_ID to str and strip whitespace
    adata.obs[SAMPLE_ID] = adata.obs[SAMPLE_ID].astype(str).str.strip()
    print(f"Cast {SAMPLE_ID} to str and stripped whitespace.")

    # If an ID map was provided, restrict to mapped participants and remap IDs
    if ID_MAP is not None:
        id_map_df = pd.read_csv(ID_MAP)
        id_map_df['nucseq_participant_id'] = id_map_df['nucseq_participant_id'].astype(str).str.strip()
        id_map_df['vcf_sample_id'] = id_map_df['vcf_sample_id'].astype(str).str.strip()
        mapping = dict(zip(id_map_df['nucseq_participant_id'], id_map_df['vcf_sample_id']))

        n_before = adata.n_obs
        adata = adata[adata.obs[SAMPLE_ID].isin(mapping), :]
        n_after = adata.n_obs
        print(f"Restricted to {n_after} / {n_before} cells with participant IDs in id_map.")
        if n_after == 0:
            raise ValueError(f"No cells remain after restricting to participants in {ID_MAP}. "
                             f"Check that {SAMPLE_ID} values match nucseq_participant_id column.")

        adata.obs[SAMPLE_ID] = adata.obs[SAMPLE_ID].map(mapping)
        print(f"Remapped {SAMPLE_ID} to vcf_sample_id using {ID_MAP}.")

    # ONE-OFF: if a participant has multiple ages recorded, use the oldest
    age_per_participant = adata.obs.groupby(SAMPLE_ID)["age"].max()
    adata.obs["age"] = adata.obs[SAMPLE_ID].map(age_per_participant)
    print("ONE-OFF: resolved multiple ages per participant by taking max.")

    # ONE-OFF: cap age at 90
    adata.obs["age"] = adata.obs["age"].clip(upper=90)
    print("ONE-OFF: capped age at 90.")

    # ONE-OFF: normalize study names — any value containing "Estiar" → "Estiar"
    mask = adata.obs["study"].str.contains("Estiar", na=False)
    adata.obs.loc[mask, "study"] = "Estiar"
    print(f"ONE-OFF: normalized {mask.sum()} cells with study containing 'Estiar' to 'Estiar'.")

    # Validate single cell type per h5ad (pipeline assumption: one file → one cell type).
    # The cell_class value is used for output file naming so it must be unambiguous.
    cell_classes = adata.obs[CT_ID].unique().tolist()
    if len(cell_classes) != 1:
        raise ValueError(f"{file} contains {len(cell_classes)} values of {CT_ID!r}: {cell_classes}. "
                         "Each input h5ad must contain exactly one cell type.")
    cell_class = str(cell_classes[0])
    print(f"Cell class: {cell_class}")

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

    # Pseudobulk aggregation — samples with < MIN_NUM_CELLS cells are dropped
    n_samples_before = adata.obs[SAMPLE_ID].nunique()
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
    n_samples_after = pdata.n_obs
    print(f"pseudobulk: {n_samples_after} / {n_samples_before} samples kept "
          f"(dropped {n_samples_before - n_samples_after} with < {MIN_NUM_CELLS} cells or < {MIN_COUNTS} counts)")

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

    # Name outputs by cell_class from adata.obs[CT_ID], not the h5ad filename stem.
    # This avoids project-specific filename prefix stripping in step 2.
    output_file1 = os.path.join(OUTPUT_DIR, f"{cell_class}_expression_matrix_ds.csv")
    output_file2 = os.path.join(OUTPUT_DIR, f"{cell_class}_composition_matrix_ds.csv")

    # decoupler.get_pseudobulk constructs obs_names as "{sample_id}_{groups_col}".
    # Strip the exact known suffix rather than splitting on all underscores, which breaks
    # if participant IDs or cell class names themselves contain underscores.
    suffix = f"_{cell_class}"
    for idx in data.index:
        if not idx.endswith(suffix):
            raise ValueError(f"Pseudobulk index '{idx}' does not end with expected "
                             f"suffix '{suffix}'. Check decoupler obs_names construction.")
    data.index = [i[:-len(suffix)] for i in data.index]

    common_samples = data.index.intersection(combined_df_corrected.index)
    print(f"Number of common samples: {len(common_samples)}")

    # Filter both datasets to only include common samples
    data_aligned = data.loc[common_samples]
    combined_df_corrected_aligned = combined_df_corrected.loc[common_samples]
    print(f"Aligned data shape: {data_aligned.shape}")
    print(f"Aligned combined_df_corrected shape: {combined_df_corrected_aligned.shape}")

    data_aligned.to_csv(output_file1)
    combined_df_corrected_aligned.to_csv(output_file2)

    # Write per-sample obs metadata for use by step 02.
    # One row per surviving sample; uses SAMPLE_ID as the participant_id column.
    surviving_samples = set(common_samples)
    sample_meta = (
        adata.obs[adata.obs[SAMPLE_ID].isin(surviving_samples)]
        .drop_duplicates(subset=[SAMPLE_ID])
        .rename(columns={SAMPLE_ID: "participant_id"})
        .reset_index(drop=True)
    )
    output_metadata = os.path.join(OUTPUT_DIR, f"{cell_class}_obs_metadata.csv")
    sample_meta.to_csv(output_metadata, index=False)

    print(f"Saved {output_file1}")
    print(f"Saved {output_file2}")
    print(f"Saved {output_metadata}")
print("Processing completed.")
