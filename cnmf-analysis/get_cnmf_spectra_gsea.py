#!/usr/bin/env python3
"""
get_cnmf_spectra_gsea: Transpose a cNMF spectra file, annotate with gene symbols,
cache the result, then run GSEA on each GEP using run-gsea.R.

Usage:
    get_cnmf_spectra_gsea.py --path <path/to/spectra.txt> [--outdir DIR] [--jobs N]

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

        GEP selection
        -------------
        All K GEPs in the file are processed.  If only a QC-passing subset
        is desired, filter the file beforehand or use --rank-col selection
        in the downstream R script.

Other arguments
---------------
--outdir    Directory under which all outputs are written.  Defaults to
            <dirname(path)>/<stem_noext>/.  Useful when --path is on a
            read-only filesystem.
--jobs      Number of parallel GSEA jobs (default: n_geps).
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).parent
RSCRIPT    = SCRIPT_DIR.parent / "micromamba_root" / "envs" / "cnmf-analysis-env" / "bin" / "Rscript"
os.environ["LD_LIBRARY_PATH"] = (
    str(RSCRIPT.parent.parent / "lib")
    + (":" + os.environ["LD_LIBRARY_PATH"] if os.environ.get("LD_LIBRARY_PATH") else "")
)
GSEA_SCRIPT = SCRIPT_DIR.parent / "scripts" / "run-gsea.R"
ENSG_TO_SYMBOL = "/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run GSEA on each GEP from a cNMF spectra file."
    )
    parser.add_argument("--path", type=Path, required=True, help="Path to cNMF spectra TSV file")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: dirname(path)/stem_noext)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Number of parallel GSEA jobs (default: n_geps)",
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
        old: f"GEP_{str(i+1).zfill(width)}"
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


def build_gsea_csv(df_t, gep_cols: list, gsea_csv: Path):
    import pandas as pd

    print(f"Saving GSEA input CSV: {gsea_csv}")
    gsea_df = df_t[["gene_id"] + gep_cols].rename(columns={"gene_id": "gene"})
    gsea_csv.parent.mkdir(parents=True, exist_ok=True)
    gsea_df.to_csv(gsea_csv, index=False)


def run_one_gep(gep_col: str, gsea_csv: Path):
    cmd = [str(RSCRIPT), "--no-save", "--no-restore", "--no-site-file", "--no-init-file",
           str(GSEA_SCRIPT), "--path", str(gsea_csv), "--rank-col", gep_col]
    print(f"  Starting GSEA for {gep_col}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR running GSEA for {gep_col}:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"run-gsea.R failed for {gep_col} (exit {result.returncode})")
    print(f"  Finished GSEA for {gep_col}")
    return gep_col


def main():
    args = parse_args()
    path = args.path.resolve()

    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    # Derive cache paths
    basename_noext = path.stem  # strips one extension
    cache_dir = (args.outdir / basename_noext).resolve() if args.outdir else path.parent / basename_noext
    cache_file = cache_dir / f"{basename_noext}-transpose.txt"
    gsea_csv = cache_dir / f"{basename_noext}-gsea-input.csv"

    # Build transpose cache if missing
    if cache_file.exists():
        print(f"Transpose cache found: {cache_file}")
        import pandas as pd
        df_t = pd.read_csv(cache_file, sep="\t")
        gep_cols = [c for c in df_t.columns if c.startswith("GEP_")]
    else:
        df_t, gep_cols = build_transpose_cache(path, cache_file)

    # Build companion GSEA CSV if missing
    if not gsea_csv.exists():
        build_gsea_csv(df_t, gep_cols, gsea_csv)
    else:
        print(f"GSEA input CSV found: {gsea_csv}")

    # Run run-gsea.R for each GEP in parallel
    n_geps = len(gep_cols)
    jobs = args.jobs if args.jobs is not None else n_geps
    print(f"\nRunning GSEA on {n_geps} GEPs with {jobs} parallel jobs...\n")

    errors = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(run_one_gep, col, gsea_csv): col for col in gep_cols}
        for fut in as_completed(futures):
            col = futures[fut]
            try:
                fut.result()
            except Exception as exc:
                errors.append((col, exc))

    if errors:
        for col, exc in errors:
            print(f"FAILED: {col}: {exc}", file=sys.stderr)
        raise RuntimeError(f"{len(errors)} GEP(s) failed GSEA — see stderr above")

    print(f"\nAll GSEA jobs complete. Results in: {cache_dir}/gsea/signed/")


if __name__ == "__main__":
    main()
