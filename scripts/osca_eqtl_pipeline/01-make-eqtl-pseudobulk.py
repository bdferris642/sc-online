# 01-make-eqtl-pseudobulk.py — Aggregate raw counts per sample × cell type, normalize, and filter genes.
#
# For each .h5ad file listed in --input-files:
#   1. Loads counts layer as X; validates required obs columns exist
#   2. Removes genes/cells with zero counts; removes NaN sample/cell-type labels
#   3. Casts --sample-id to str and strips whitespace
#   4. If --id-map provided: restricts to mapped participants and remaps IDs
#   5. Validates that each h5ad contains exactly one cell type (pipeline assumption)
#   6. Pseudobulk sum aggregation with decoupler.pp.pseudobulk (decoupler >= 2.x)
#   7. Filter samples with < MIN_NUM_CELLS cells or < MIN_COUNTS total counts
#   8. CPM normalization (target_sum=1e6) + log1p
#   9. Gene filter: mean log-expression >= --gene-log-expr-threshold
#  10. Scale to max value MAX_SCALE_EXPR_VALUE
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
    default=None,
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

    if not (args.strs_to_skip is None or args.strs_to_skip == ""):
        if any(s in file for s in STRS_TO_SKIP):
            print(f"Skipping {file} (matched strs-to-skip)")
            continue

    # Load the h5ad file
    print(f"  Loading {file} ...")
    adata = sc.read_h5ad(file)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes | "
          f"layers: {list(adata.layers.keys())} | "
          f"obs cols: {list(adata.obs.columns)}")

    if "counts" not in adata.layers:
        raise ValueError(f"'counts' layer not found in {file}. "
                         f"Available layers: {list(adata.layers.keys())}")
    adata.X = adata.layers["counts"]
    print(f"  X set to 'counts' layer (sparse={hasattr(adata.X, 'toarray')})")

    # Validate required obs columns exist before any downstream use
    for col in [SAMPLE_ID, CT_ID]:
        if col not in adata.obs.columns:
            raise ValueError(f"{col!r} not found in adata.obs. "
                             f"Available columns: {list(adata.obs.columns)}")

    gene_sums = np.asarray(adata.X.sum(axis=0)).ravel()
    cell_sums = np.asarray(adata.X.sum(axis=1)).ravel()
    print(f"  Removing {(gene_sums == 0).sum()} genes with zero counts")
    print(f"  Removing {(cell_sums == 0).sum()} cells with zero counts")

    adata = adata[:, gene_sums > 0]  # Remove genes with 0 expression
    adata = adata[cell_sums > 0, :]  # Remove cells with 0 expression
    print(f"  After zero-count filter: {adata.n_obs} cells x {adata.n_vars} genes")

    num_nas_sample = adata.obs[SAMPLE_ID].isna().sum()
    num_nas_ct = adata.obs[CT_ID].isna().sum()
    print(f"  Removing {num_nas_sample} cells with missing {SAMPLE_ID!r} and "
          f"{num_nas_ct} cells with missing {CT_ID!r}")
    adata = adata[~adata.obs[SAMPLE_ID].isna() & ~adata.obs[CT_ID].isna(), :]
    print(f"  After NaN filter: {adata.n_obs} cells x {adata.n_vars} genes")

    # Cast SAMPLE_ID to str and strip whitespace
    adata.obs[SAMPLE_ID] = adata.obs[SAMPLE_ID].astype(str).str.strip()
    print(f"  Cast {SAMPLE_ID!r} to str and stripped whitespace | "
          f"{adata.obs[SAMPLE_ID].nunique()} unique samples")

    # If an ID map was provided, restrict to mapped participants and remap IDs
    if ID_MAP is not None:
        id_map_df = pd.read_csv(ID_MAP)
        id_map_df['nucseq_participant_id'] = id_map_df['nucseq_participant_id'].astype(str).str.strip()
        id_map_df['vcf_sample_id'] = id_map_df['vcf_sample_id'].astype(str).str.strip()
        mapping = dict(zip(id_map_df['nucseq_participant_id'], id_map_df['vcf_sample_id']))

        n_before = adata.n_obs
        adata = adata[adata.obs[SAMPLE_ID].isin(mapping), :]
        n_after = adata.n_obs
        print(f"  ID map: restricted to {n_after} / {n_before} cells with participant IDs in id_map.")
        if n_after == 0:
            raise ValueError(f"No cells remain after restricting to participants in {ID_MAP}. "
                             f"Check that {SAMPLE_ID} values match nucseq_participant_id column.")

        adata.obs[SAMPLE_ID] = adata.obs[SAMPLE_ID].map(mapping)
        print(f"  Remapped {SAMPLE_ID!r} to vcf_sample_id using {ID_MAP}.")

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

    # ONE-OFF: if brain_bank contains "aryland", normalize to "Maryland"
    mask = adata.obs["brain_bank"].str.contains("aryland", na=False)
    adata.obs.loc[mask, "brain_bank"] = "Maryland"
    print(f"ONE-OFF: normalized {mask.sum()} cells with brain_bank containing 'aryland' to 'Maryland'.")

    # Determine cell type. If multiple are present, keep only the majority class.
    cell_classes = adata.obs[CT_ID].unique().tolist()
    if len(cell_classes) != 1:
        majority_class = adata.obs[CT_ID].value_counts().idxmax()
        n_majority = (adata.obs[CT_ID] == majority_class).sum()
        print(f"WARNING: {file} contains {len(cell_classes)} cell types: {cell_classes}. "
              f"Keeping majority class '{majority_class}' ({n_majority}/{len(adata)} cells).")
        adata = adata[adata.obs[CT_ID] == majority_class].copy()
    cell_class = str(adata.obs[CT_ID].iloc[0])
    # Sanitize for use in filenames — spaces break shell pipelines in step 3.
    # Keep original cell_class for decoupler suffix stripping below (obs_names use the raw value).
    cell_class_safe = cell_class.replace(" ", "_")
    print(f"  Cell class: {cell_class!r} (safe name: {cell_class_safe!r})")

    output_file1    = os.path.join(OUTPUT_DIR, f"{cell_class_safe}_expression_matrix_ds.csv")
    output_file2    = os.path.join(OUTPUT_DIR, f"{cell_class_safe}_composition_matrix_ds.csv")
    output_metadata = os.path.join(OUTPUT_DIR, f"{cell_class_safe}_obs_metadata.csv")
    if all(os.path.exists(f) for f in [output_file1, output_file2, output_metadata]):
        print(f"  Skipping {cell_class!r}: all 3 output CSVs already exist.")
        continue

    print(f"  Pre-pseudobulk adata shape: {adata.shape}")

    # Composition matrix: cell count per sample (one column = this cell type).
    # Indexed by SAMPLE_ID; used later to align with the expression matrix.
    cell_counts = adata.obs.groupby(SAMPLE_ID).size().rename(cell_class)
    sample_obs = (
        adata.obs.drop_duplicates(subset=[SAMPLE_ID])
        .set_index(SAMPLE_ID)
    )
    combined_df_corrected = pd.concat([sample_obs, cell_counts], axis=1)
    print(f"  Composition matrix shape: {combined_df_corrected.shape}")

    # Pseudobulk aggregation (decoupler >= 2.x API)
    n_samples_before = adata.obs[SAMPLE_ID].nunique()
    print(f"  Running pseudobulk aggregation on {n_samples_before} samples ...")
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col=SAMPLE_ID,
        groups_col=CT_ID,
        skip_checks=True,
        mode="sum",
    )
    print(f"  Pseudobulk raw shape: {pdata.n_obs} samples x {pdata.n_vars} genes")

    # Filter samples with too few cells or too few counts (was done inside get_pseudobulk in v1)
    keep = (pdata.obs["psbulk_cells"] >= MIN_NUM_CELLS) & (pdata.obs["psbulk_counts"] >= MIN_COUNTS)
    pdata = pdata[keep].copy()
    n_samples_after = pdata.n_obs
    print(f"  Sample QC: {n_samples_after} / {n_samples_before} kept "
          f"(dropped {n_samples_before - n_samples_after} with "
          f"< {MIN_NUM_CELLS} cells or < {MIN_COUNTS} counts)")

    # Normalize total counts and log-transform
    print(f"  Normalizing (CPM) and log1p-transforming ...")
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)

    # Filter genes with mean log-transformed expression >= GENE_LOG_EXPR_THRESHOLD
    gene_filter = np.asarray(pdata.X.mean(axis=0)).ravel() >= GENE_LOG_EXPR_THRESHOLD
    n_genes_before = pdata.n_vars
    pdata = pdata[:, gene_filter]
    print(f"  Gene expression filter (mean log-expr >= {GENE_LOG_EXPR_THRESHOLD}): "
          f"{pdata.n_vars} / {n_genes_before} genes kept")

    # Scale the data to a maximum value
    print(f"  Scaling (max_value={MAX_SCALE_EXPR_VALUE}) ...")
    sc.pp.scale(pdata, max_value=MAX_SCALE_EXPR_VALUE)

    # Save the processed data as a CSV file
    data = pd.DataFrame(pdata.X, index=pdata.obs_names, columns=pdata.var_names)
    print(f"  Final pseudobulk matrix: {data.shape[0]} samples x {data.shape[1]} genes")

    # decoupler.pp.pseudobulk constructs obs_names as "{sample_id}_{groups_col_value}".
    # Strip the exact known suffix rather than splitting on all underscores, which breaks
    # if participant IDs or cell class names themselves contain underscores.
    suffix = f"_{cell_class}"
    for idx in data.index:
        if not idx.endswith(suffix):
            raise ValueError(f"Pseudobulk index '{idx}' does not end with expected "
                             f"suffix '{suffix}'. Check decoupler obs_names construction.")
    data.index = [i[:-len(suffix)] for i in data.index]

    common_samples = data.index.intersection(combined_df_corrected.index)
    print(f"  Common samples between expression and composition: {len(common_samples)}")

    # Filter both datasets to only include common samples
    data_aligned = data.loc[common_samples]
    combined_df_corrected_aligned = combined_df_corrected.loc[common_samples]
    print(f"  Aligned expression matrix: {data_aligned.shape}")
    print(f"  Aligned composition matrix: {combined_df_corrected_aligned.shape}")

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
    sample_meta.to_csv(output_metadata, index=False)

    print(f"  Saved expression matrix ({data_aligned.shape}): {output_file1}")
    print(f"  Saved composition matrix ({combined_df_corrected_aligned.shape}): {output_file2}")
    print(f"  Saved obs metadata ({sample_meta.shape}): {output_metadata}")
print("Processing completed.")
