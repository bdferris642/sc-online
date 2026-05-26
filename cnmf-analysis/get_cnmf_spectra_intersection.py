#!/usr/bin/env python3
"""
get_cnmf_spectra_intersection: Transpose a cNMF spectra file, annotate with gene symbols,
cache the result, then run Fisher's exact test gene-set intersection on each GEP using
run_intersection.R. After all per-GEP jobs finish, produces cross-GEP summary heatmaps.

Usage:
    get_cnmf_spectra_intersection.py --path <path/to/spectra.txt> [--outdir DIR] [--top-n N]
                                     [--zscore K] [--jobs N] [--min-size M] [--max-size M]

Input requirements
------------------
--path  Path to a single cNMF *spectra* consensus file.

        Required format
        ---------------
        Tab-separated text file produced by cNMF, named like:
            <run>.spectra.k_<K>.dt_<dt>.consensus.txt
        or the z-scored variant:
            <run>.gene_spectra_score.k_<K>.dt_<dt>.txt

        Layout
        ------
        • Rows    = GEPs (gene expression programs), indexed 1 .. K.
                    The first column is an unnamed integer index column.
        • Columns = Ensembl gene IDs (e.g. ENSG00000187583).
        • Values  = numeric gene loading scores (float).

        This is the GEP × gene matrix — NOT any of the following:
          ✗ Manifest / pointer files (multi-experiment metadata TSVs)
          ✗ Usages / embedding files  (<run>.usages.k_*.consensus.txt)
              — those are cells × GEPs, not GEPs × genes
          ✗ top_genes CSVs or sample heatmap files

        Example correct path (one file per cNMF run):
            /path/to/cNMF/my-run/k12/my-run.spectra.k_12.dt_0_05.consensus.txt

        Gene selection (per GEP)
        ------------------------
        Genes are selected by z-score threshold (default: mean + 3·sd) across the
        GEP's loading distribution, then capped at --top-n (default 150) by
        highest loading.  Gene names are resolved via HGNC symbol; rows without a
        matching symbol are excluded from the gene-set background.

Other arguments
---------------
--outdir    Directory under which all outputs are written.  Defaults to
            <dirname(path)>/<stem_noext>/.  Useful when --path is on a
            read-only filesystem.
--top-n     Maximum number of genes selected per GEP after z-score filter
            (default: 150).
--zscore    Z-score threshold for gene selection, i.e. genes with loading >
            mean + zscore·sd are selected (default: 3.0).
--jobs      Parallel R jobs (default: one per GEP).
--min-size  Minimum gene-set size after background filtering (default: 15).
--max-size  Maximum gene-set size after background filtering (default: 250).
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR          = Path(__file__).parent
RSCRIPT             = SCRIPT_DIR.parent / "micromamba_root" / "envs" / "cnmf-analysis-env" / "bin" / "Rscript"
os.environ["LD_LIBRARY_PATH"] = (
    str(RSCRIPT.parent.parent / "lib")
    + (":" + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")
)
INTERSECTION_SCRIPT = SCRIPT_DIR / "run_intersection.R"
ENSG_TO_SYMBOL      = "/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv"

GENE_SET_COLLECTIONS = [
    "kegg_2026_human",
    "reactome_2024",
    "gwas_2025",
    "jensen_disease_2025",
    "trrust_2019",
    "jaspar_tf_2025",
    "go_bp_2025",
    "go_mf_2025",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run gene-set intersection analysis on each GEP from a cNMF spectra file."
    )
    parser.add_argument("--path", type=Path, required=True, help="Path to cNMF spectra TSV file")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: dirname(path)/stem_noext)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=150,
        help="Cap on selected genes per GEP (default: 150)",
    )
    parser.add_argument(
        "--zscore",
        type=float,
        default=3.0,
        help="Select genes with loading > mean + K*std (default: 3.0)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Number of parallel jobs (default: n_geps)",
    )
    parser.add_argument(
        "--min-size",
        type=int,
        default=15,
        help="Min gene set size (default: 15)",
    )
    parser.add_argument(
        "--max-size",
        type=int,
        default=250,
        help="Max gene set size (default: 250)",
    )
    return parser.parse_args()


def build_transpose_cache(path: Path, cache_file: Path):
    import pandas as pd

    print(f"Loading spectra file: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    # df: rows=GEPs, cols=Ensembl IDs

    n_geps = df.shape[0]
    width = len(str(n_geps))

    # Transpose: rows become genes, columns become GEPs
    df_t = df.T
    df_t.index.name = "gene_id"

    # Rename GEP columns to zero-padded names
    new_col_names = {
        old: f"GEP_{str(i + 1).zfill(width)}"
        for i, old in enumerate(df_t.columns)
    }
    df_t = df_t.rename(columns=new_col_names)

    # Left join with gene symbol table
    print(f"Loading gene symbol map: {ENSG_TO_SYMBOL}")
    sym_df = pd.read_csv(ENSG_TO_SYMBOL, usecols=["ensembl_gene_id", "hgnc_symbol"])
    sym_df = sym_df.rename(columns={"ensembl_gene_id": "gene_id", "hgnc_symbol": "gene_name"})
    sym_df = sym_df.drop_duplicates(subset="gene_id", keep="first")

    df_t = df_t.reset_index()  # gene_id becomes a column
    df_t = df_t.merge(sym_df, on="gene_id", how="left")

    gep_cols = [c for c in df_t.columns if c.startswith("GEP_")]
    df_t = df_t[["gene_name", "gene_id"] + gep_cols]

    print(f"Saving transpose cache: {cache_file}")
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    df_t.to_csv(cache_file, sep="\t", index=False)

    return df_t, gep_cols


def build_intersection_input(df_t, gep_cols: list, intersection_input: Path):
    import pandas as pd

    print(f"Saving intersection input TSV: {intersection_input}")
    intersection_input.parent.mkdir(parents=True, exist_ok=True)
    out_df = df_t[["gene_name", "gene_id"] + gep_cols]
    out_df.to_csv(intersection_input, sep="\t", index=False)


def run_one_gep(
    gep_col: str,
    intersection_input: Path,
    top_n: int,
    zscore: float,
    min_size: int,
    max_size: int,
):
    cmd = [
        str(RSCRIPT), "--no-save", "--no-restore", "--no-site-file", "--no-init-file",
        str(INTERSECTION_SCRIPT),
        "--path", str(intersection_input),
        "--rank-col", gep_col,
        "--top-n", str(top_n),
        "--zscore", str(zscore),
        "--min-size", str(min_size),
        "--max-size", str(max_size),
    ]
    print(f"  Starting intersection for {gep_col}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(
            f"  ERROR running intersection for {gep_col}:\n{result.stderr}",
            file=sys.stderr,
        )
        raise RuntimeError(
            f"run_intersection.R failed for {gep_col} (exit {result.returncode})"
        )
    print(f"  Finished intersection for {gep_col}")
    return gep_col


def build_cross_gep_heatmaps(
    intersection_outdir: Path,
    gep_cols: list,
    basename_noext: str,
    cross_gep_outdir: Path,
):
    """Build cross-GEP significance and odds ratio heatmaps for each gene set collection."""
    import pandas as pd
    import numpy as np

    cross_gep_outdir.mkdir(parents=True, exist_ok=True)

    for collection in GENE_SET_COLLECTIONS:
        print(f"  Building cross-GEP heatmaps for: {collection}")

        # Collect per-GEP result CSVs for this collection
        records = {}
        for gep_col in gep_cols:
            csv_path = (
                intersection_outdir
                / gep_col
                / f"{basename_noext}_{gep_col}_{collection}_intersection.csv"
            )
            if not csv_path.exists():
                continue
            df = pd.read_csv(csv_path)
            if df.empty:
                continue
            records[gep_col] = df

        if not records:
            print(f"    No results found for {collection}, skipping.")
            continue

        # Gather all pathways significant (padj < 0.05) in ≥1 GEP
        sig_pathways = set()
        for gep_col, df in records.items():
            sig = df.loc[df["padj"] < 0.05, "pathway"]
            sig_pathways.update(sig.tolist())

        if not sig_pathways:
            print(f"    No significant pathways for {collection}, skipping.")
            continue

        sig_pathways = sorted(sig_pathways)

        # Build matrices
        neg_log10_padj_mat = pd.DataFrame(index=sig_pathways, columns=gep_cols, dtype=float)
        log2_or_mat = pd.DataFrame(index=sig_pathways, columns=gep_cols, dtype=float)

        neg_log10_padj_mat[:] = 0.0
        log2_or_mat[:] = float("nan")

        for gep_col, df in records.items():
            df_sig = df[df["pathway"].isin(sig_pathways)].set_index("pathway")
            for pathway in sig_pathways:
                if pathway not in df_sig.index:
                    continue
                row = df_sig.loc[pathway]
                padj = row["padj"]
                or_val = row["odds_ratio"]
                nlp = min(-np.log10(padj) if padj > 0 else 10.0, 10.0)
                neg_log10_padj_mat.loc[pathway, gep_col] = nlp if padj < 0.05 else 0.0
                if padj < 0.05 and or_val > 0:
                    log2_or_mat.loc[pathway, gep_col] = np.log2(or_val) if or_val > 0 else float("nan")

        _save_heatmap(
            mat=neg_log10_padj_mat,
            outdir=cross_gep_outdir,
            filename_stem=f"{collection}_cross_gep_neg_log10_padj",
            title=f"{collection}\n-log10(padj)",
            cmap="YlOrRd",
            center=None,
            vmin=0,
            vmax=10,
            mask=None,
        )

        # OR heatmap: mask where padj >= 0.05 (NaN entries)
        or_mask = log2_or_mat.isna()
        log2_or_mat_filled = log2_or_mat.fillna(0.0)

        _save_heatmap(
            mat=log2_or_mat_filled,
            outdir=cross_gep_outdir,
            filename_stem=f"{collection}_cross_gep_log2_OR",
            title=f"{collection}\nlog2(Odds Ratio)",
            cmap="RdBu_r",
            center=0,
            vmin=None,
            vmax=None,
            mask=or_mask,
        )


def _save_heatmap(mat, outdir, filename_stem, title, cmap, center, vmin, vmax, mask):
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    svg_path = outdir / f"{filename_stem}.svg"
    csv_path = outdir / f"{filename_stem}.csv"

    # Save companion CSV first (always)
    mat.to_csv(csv_path)

    n_rows, n_cols = mat.shape
    if n_rows == 0 or n_cols == 0:
        print(f"    Empty matrix for {filename_stem}, skipping figure.")
        return

    # Cluster rows; keep columns in order
    row_vals = mat.values.astype(float)
    row_vals_clean = np.nan_to_num(row_vals, nan=0.0)

    if n_rows > 1:
        try:
            row_dist = pdist(row_vals_clean, metric="euclidean")
            row_linkage = linkage(row_dist, method="average")
            row_order = leaves_list(row_linkage)
        except Exception:
            row_order = list(range(n_rows))
    else:
        row_order = [0]

    mat_ordered = mat.iloc[row_order, :]

    fig_height = max(4, min(n_rows * 0.3 + 2, 60))
    fig_width = max(4, min(n_cols * 0.7 + 3, 40))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    mask_ordered = mask.iloc[row_order, :] if mask is not None else None

    kwargs = dict(
        data=mat_ordered.astype(float),
        cmap=cmap,
        ax=ax,
        linewidths=0.3,
        linecolor="white",
    )
    if center is not None:
        kwargs["center"] = center
    if vmin is not None:
        kwargs["vmin"] = vmin
    if vmax is not None:
        kwargs["vmax"] = vmax
    if mask_ordered is not None:
        kwargs["mask"] = mask_ordered

    sns.heatmap(**kwargs)
    ax.set_title(title, fontsize=10, pad=8)
    ax.set_xlabel("GEP", fontsize=9)
    ax.set_ylabel("Pathway", fontsize=9)
    ax.tick_params(axis="y", labelsize=6)
    ax.tick_params(axis="x", labelsize=8, rotation=45)

    plt.tight_layout()
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {svg_path}")
    print(f"    Saved: {csv_path}")


def main():
    args = parse_args()
    path = args.path.resolve()

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    # Derive cache paths
    basename_noext = path.stem  # strips one extension
    cache_dir = (args.outdir / basename_noext).resolve() if args.outdir else path.parent / basename_noext
    cache_file = cache_dir / f"{basename_noext}-transpose.txt"
    intersection_input = cache_dir / f"{basename_noext}-intersection-input.tsv"
    intersection_outdir = cache_dir / "intersection"
    cross_gep_outdir = intersection_outdir / "cross_gep"

    # Build transpose cache if missing (reused from GSEA wrapper if present)
    if cache_file.exists():
        print(f"Transpose cache found: {cache_file}")
        import pandas as pd
        df_t = pd.read_csv(cache_file, sep="\t")
        gep_cols = [c for c in df_t.columns if c.startswith("GEP_")]
    else:
        df_t, gep_cols = build_transpose_cache(path, cache_file)

    # Build intersection input TSV if missing
    if not intersection_input.exists():
        build_intersection_input(df_t, gep_cols, intersection_input)
    else:
        print(f"Intersection input TSV found: {intersection_input}")

    # Run run_intersection.R for each GEP in parallel
    n_geps = len(gep_cols)
    jobs = args.jobs if args.jobs is not None else n_geps
    print(f"\nRunning intersection on {n_geps} GEPs with {jobs} parallel jobs...\n")

    errors = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(
                run_one_gep,
                col,
                intersection_input,
                args.top_n,
                args.zscore,
                args.min_size,
                args.max_size,
            ): col
            for col in gep_cols
        }
        for fut in as_completed(futures):
            col = futures[fut]
            try:
                fut.result()
            except Exception as exc:
                errors.append((col, exc))

    if errors:
        for col, exc in errors:
            print(f"FAILED: {col}: {exc}", file=sys.stderr)
        raise RuntimeError(f"{len(errors)} GEP(s) failed intersection — see stderr above")

    print(f"\nAll intersection jobs complete.")

    # Build cross-GEP summary heatmaps
    print(f"\nBuilding cross-GEP heatmaps in: {cross_gep_outdir}\n")
    build_cross_gep_heatmaps(
        intersection_outdir=intersection_outdir,
        gep_cols=gep_cols,
        basename_noext=basename_noext,
        cross_gep_outdir=cross_gep_outdir,
    )

    print(f"\nDone. Results in: {intersection_outdir}")


if __name__ == "__main__":
    main()
