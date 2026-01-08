import argparse
import os
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from cellbender.remove_background.downstream import anndata_from_h5

from constants import (
    ADATA_OUT_BASENAME,
    CB_OUT_FILTERED_BASENAME_TEMPLATE,
    DROPSIFT_BASENAME,
    MOLECULE_INFO_BASENAME,
    VIREO_OUT_BASENAME,
)


def df_decode_bytes_inplace(df: pd.DataFrame) -> None:
    """Convert any bytes/bytearray in object columns to real strings, in place."""
    for c in df.select_dtypes(include="object").columns:
        col = df[c]
        if col.map(lambda x: isinstance(x, (bytes, bytearray))).any():
            df[c] = col.map(
                lambda x: x.decode() if isinstance(x, (bytes, bytearray)) else x
            )


def load_single_library(basepath: str, bcl: str, rna_index: str) -> None:
    """
    Loads a single library to anndata, merging in molecule info, dropsift, and vireo data.
    Saves that anndata to disk in path/to/<bcl>/<rna_index>/adata.h5ad

    Assumes directory structure
    ─── gene-expression
        └── A22C7WYLT3 <-- a BCL
            └── pCalico_GTExsHSrSNE9iDAPId231120 <-- an RNA index
                ├── cellbender
                │   └── pCalico_GTExsHSrSNE9iDAPId231120_out_filtered.h5 <-- CellBender output h5 assumed to have RNA index embedded in name
                ├── dropsift_outputs
                │   └── dropsift_output.csv
                ├── outs
                │   └── molecule_info.h5
                └── vireo <-- optional vireo directory
                    └── donor_ids.tsv

    We only gave vireo assignments to cells that passed CellRanger filtering
    We will therefore only consider cells that passed both CellBender and CellRanger filtering

    While downstream we will only consider cells that were unambiguously assigned to a single donor by vireo
    For this step, heterotypic doublets are retained to aid in helping identify homotypic doublets

    Similarly, no hard filtering on nUMI, pct_intronic, pct_mito, or dropsift is_cell etc. are applied at this stage

    Returns
    -------
    None
        This function does not return a value; it saves the anndata object to disk
    """

    h5_filt_path = f"{basepath}/{bcl}/{rna_index}/{CB_OUT_FILTERED_BASENAME_TEMPLATE.substitute(index = rna_index)}"
    molecule_info_path = f"{basepath}/{bcl}/{rna_index}/{MOLECULE_INFO_BASENAME}"
    dropsift_path = f"{basepath}/{bcl}/{rna_index}/{DROPSIFT_BASENAME}"
    vireo_path = f"{basepath}/{bcl}/{rna_index}/{VIREO_OUT_BASENAME}"
    adata_out_path = f"{basepath}/{bcl}/{rna_index}/{ADATA_OUT_BASENAME}"

    adata = anndata_from_h5(h5_filt_path)

    # make var names unique
    adata.var["gene_name"] = adata.var_names.values
    adata.var_names = adata.var.gene_id.values

    # cell barcodes can collide across libraries, so we use library__barcode as globally unique obs_names
    adata.obs["barcode"] = adata.obs_names.values
    adata.obs["bcl"] = bcl
    adata.obs["rna_index"] = rna_index
    adata.obs["library"] = f"{bcl}__{rna_index}"
    adata.obs["library__barcode"] = adata.obs["library"] + "__" + adata.obs["barcode"]
    adata.obs_names = adata.obs["library__barcode"]

    # frac mito is fraction of counts coming from mitochondrial genes (all of whose names start with "MT-")
    # cell x gene, so sum within cells across genes is axis = 1
    adata.obs["frac_mito"] = (
        np.array(
            adata.X[:, [g.startswith("MT-") for g in adata.var.gene_name]].sum(axis=1)
        ).flatten()
        / np.array(adata.X.sum(axis=1)).flatten()
    )

    df_decode_bytes_inplace(adata.obs)
    df_decode_bytes_inplace(adata.var)

    # Read Molecule Info HDF5
    with h5py.File(molecule_info_path, "r") as f:
        barcodes_all = f["barcodes"].asstr()[:]  # (N_barcodes,) str
        features_all = f["features/name"].asstr()[:]  # (N_features,) str
        barcode_idx = f["barcode_idx"][:]  # (N_molecules,) int
        feature_idx = f["feature_idx"][:]  # (N_molecules,) int
        umi = f["umi"][:]  # (N_molecules,) int
        umi_type = f["umi_type"][:]  # (N_molecules,) int
        count = f["count"][:]  # (N_molecules,) int

    # Map indices to strings (R needed to add +1; Python uses 0-based indexing so no +1 here)
    mol_barcode = barcodes_all[barcode_idx]
    mol_feature = features_all[feature_idx]

    # 10x barcodes have a "-1" suffix
    mol_barcode = np.char.add(mol_barcode, "-1")
    mol_df = pd.DataFrame(
        {
            "barcode": mol_barcode,
            "feature": mol_feature,
            "umi": umi,
            "umi_type": umi_type.astype(int),
            "count": count.astype(int),
        }
    )

    # Keep only molecules whose barcodes are present in the filtered matrix
    mol_df = mol_df[mol_df["barcode"].isin(adata.obs.barcode)]

    # Summarize per barcode
    # mol_info_nUMI = number of distinct molecule rows for that barcode
    # mol_info_nRead = sum of 'count'
    # frac_intronic = fraction with umi_type == 0
    mol_grouped = mol_df.groupby("barcode").agg(
        mol_info_nUMI=("umi", "size"),
        mol_info_nRead=("count", "sum"),
        frac_intronic=("umi_type", lambda x: (x == 0).mean()),
    )
    df_decode_bytes_inplace(mol_grouped)

    # Join into adata.obs on barcode column
    adata.obs = pd.merge(
        adata.obs, mol_grouped, left_on="barcode", right_index=True, how="left"
    )

    adata.obs[["mol_info_nUMI", "mol_info_nRead", "frac_intronic"]] = adata.obs[
        ["mol_info_nUMI", "mol_info_nRead", "frac_intronic"]
    ].fillna({"mol_info_nUMI": 0, "mol_info_nRead": 0, "frac_intronic": 0.0})

    adata.obs["mol_info_nUMI"] = adata.obs["mol_info_nUMI"].astype(int)
    adata.obs["mol_info_nRead"] = adata.obs["mol_info_nRead"].astype(int)

    # Read vireo output
    # left join vireo df to obs on barcode / obs_names
    try: 
        vireo_df = pd.read_csv(vireo_path, sep="\t")
        df_decode_bytes_inplace(vireo_df)

        # Force prefix on everything except 'cell' and 'donor_id'
        vireo_keep = {"cell", "donor_id"}
        vireo_old_cols = list(vireo_df.columns)
        vireo_df = vireo_df.set_axis(
            [c if c in vireo_keep else f"vireo_{c}" for c in vireo_old_cols],
            axis=1,
            copy=False,
        )

        # Verify the change actually happened on this object
        assert any(
            o != n for o, n in zip(vireo_old_cols, vireo_df.columns)
        ), "Rename did not apply to vireo_df!"

        adata.obs = pd.merge(
            adata.obs, vireo_df, left_on="barcode", right_on="cell", how="left"
        )
    except FileNotFoundError:
        print(f"Vireo output not found at {vireo_path}, skipping vireo merge.")

    # Read dropsift output
    # left join dropsift df to obs on barcode
    dropsift_df = pd.read_csv(dropsift_path)
    df_decode_bytes_inplace(dropsift_df)

    # dropsift has this column training_label_is_cell which is boolean with NA for unlabelled cells
    # must cast to string to make write_h5ad happy
    dropsift_df.training_label_is_cell = dropsift_df.training_label_is_cell.fillna(
        "Test"
    ).astype("str")

    dropsift_keep = {"cell_barcode"}
    dropsift_old_cols = list(dropsift_df.columns)
    dropsift_df = dropsift_df.set_axis(
        [c if c in dropsift_keep else f"dropsift_{c}" for c in dropsift_old_cols],
        axis=1,
        copy=False,
    )
    assert any(
        o != n for o, n in zip(dropsift_old_cols, dropsift_df.columns)
    ), "Rename did not apply to dropsift_df!"

    adata.obs = pd.merge(
        adata.obs,
        dropsift_df,
        left_on="barcode",
        right_on="cell_barcode",
        how="left",
        suffixes=(None, "_dropsift"),
    )

    # drop duplicative columns.
    adata.obs.drop(
        columns=[
            "cell",
            "cell_barcode",
            "dropsift_pct_intronic",
            "dropsift_pct_mt",
            "dropsift_num_retained_transcripts",
            "dropsift_num_transcripts",
        ],
        inplace=True,
        errors='ignore'
    )

    adata.obs.nUMI = adata.X.sum(axis=1).A1

    # now remove cells from adata that have a NA in donor_id (i.e. not assigned by vireo)
    # only do this if vireo was run in the first place
    if "donor_id" in adata.obs.columns:
        adata = adata[~adata.obs["donor_id"].isna(), :]

    # Suppress only during write_h5ad, as these are known to be harmless (e.g., dtype conversion warnings).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata.write_h5ad(adata_out_path)

def load_libraries_from_path(
    basepath: str,
    bcl_index_file: str = None,
    bcl_index: str = None,
) -> None:
    """
    Example usage:
    load_multiple_libraries("/path/to/gene-expression")

    By default, loads all bcls found under basepath and all rna_indices under each bcl
    Optionally specify a bcl_index_file with lines of the form `<BCL>/<RNA_INDEX>` to only load specific libraries
    Or optionally specify a single bcl_index string of the form `<BCL>/<RNA_INDEX>` to only load that library

    Returns
    -------
    None
        This function does not return a value; one or more anndata objects are saved to disk
    """

    if not (bcl_index is None) and not (bcl_index_file is None):
        raise ValueError(
            "Cannot specify both bcl_index and bcl_index_file; please choose one."
        )
    if bcl_index_file is not None:
        # read lines of file
        with open(bcl_index_file, "r") as f:
            lines = f.readlines()
        for line in lines:
            bcl = line.strip().split("/")[0]
            rna_index = line.strip().split("/")[1]
            print(f"\nLoading library: {bcl} / {rna_index}\n")
            load_single_library(basepath, bcl, rna_index)
            print(f"\nFinished library: {bcl} / {rna_index}\n")
    elif bcl_index is not None:
        bcl = bcl_index.split("/")[0]
        rna_index = bcl_index.split("/")[1]
        print(f"\nLoading library: {bcl} / {rna_index}\n")
        load_single_library(basepath, bcl, rna_index)
        print(f"\nFinished library: {bcl} / {rna_index}\n")
    else:
        bcls = [
            d for d in os.listdir(basepath) if os.path.isdir(os.path.join(basepath, d))
        ]
        for bcl in bcls:
            rna_indices = [
                d
                for d in os.listdir(os.path.join(basepath, bcl))
                if os.path.isdir(os.path.join(basepath, bcl, d))
            ]
            for rna_index in rna_indices:
                print(f"\nLoading library: {bcl} / {rna_index}\n")
                load_single_library(basepath, bcl, rna_index)
                print(f"\nFinished library: {bcl} / {rna_index}\n")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Load libraries from basepath to anndata.h5ad files"
    )
    ap.add_argument(
        "--basepath", required=True, type=str, help="Path to gene-expression directory"
    )
    ap.add_argument(
        "--bcl-index-file",
        required=False,
        type=str,
        default=None,
        help="Optional path to file listing BCLs and RNA indices to process",
    )
    ap.add_argument(
        "--bcl-index",
        required=False,
        type=str,
        default=None,
        help="Optional single BCL/RNA index to process",
    )
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    load_libraries_from_path(
        basepath=args.basepath,
        bcl_index_file=args.bcl_index_file,
        bcl_index=args.bcl_index,
    )
