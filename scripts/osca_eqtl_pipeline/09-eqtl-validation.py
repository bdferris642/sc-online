#!/usr/bin/env python3
"""
09-eqtl-validation.py — eQTL validation pipeline (Step 9)

Self-activating: re-execs into osca-venv if not already running in it.

Runs 10 validation analyses per cell class:
  1  qq_lambda        — QQ plot + genomic inflation (λ_GC) by MAF bin
  2  pi1_replication  — π₁ replication rate vs GTEx SN eQTLs
  3  effect_size_maf  — Effect size |β| vs MAF
  4  tss_enrichment   — Enrichment of sig SNPs near gene TSS
  5  atac_enrichment  — Overlap with Corces 2020 brain ATAC peaks
  6  egene_chrom      — eGenes per chromosome (raw + normalized)
  7  go_enrichment    — GO-BP and GO-MF enrichment of eGenes
  8  cis_distance     — Cis distance of lead eQTL SNPs
  9  eqtl_per_gene    — eQTL count distribution per eGene
 10  top_snp_fraction — Fraction of sig SNPs at rank 1 within gene

v2 changes vs v1:
  - --expr-csv removed; replaced with --pb-output-dir
  - Expression CSV auto-discovered: {pb_output_dir}/{CC_safe}_expression_matrix_ds.csv
  - Mean expression (per ENSG) computed inline; saved as {cc_out}/mean_expr.csv
  - _load_gene_loc() reads gene_loc_v2.txt (ensg_id, chr, TSS, gene_symbol, strand)
  - RDS column ensg_id used throughout (was Gene/Probe with Entrez/symbol)
  - Parallelism moved to bash orchestrator: each process handles one CC (--cell-class)

Usage:
  python 09-eqtl-validation.py \\
    --eqtl-dir      /path/to/eqtl_final_outs/my_run \\
    --gene-loc      /path/to/gene_loc_v2.txt \\
    --pb-output-dir /path/to/pb_output_dir \\
    --out-dir       /path/to/validation \\
    [--cell-class   da_neuron] \\
    [--gtex-sn      /path/to/gtex_sn_signif_pairs.txt.gz] \\
    [--atac-bed     /path/to/corces_2020_da_atac_peaks.bed.gz] \\
    [--go-bp-gmt    /path/to/GO_Biological_Process_2025.gmt] \\
    [--go-mf-gmt    /path/to/GO_Molecular_Function_2025.gmt] \\
    [--padj-thresh  0.05]
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

# ── Self-activation ────────────────────────────────────────────────────────────
SANDBOX = Path(__file__).resolve().parents[1]
ENV_BIN = SANDBOX / "micromamba_root" / "envs" / "osca-venv" / "bin"
if not ENV_BIN.is_dir():
    ENV_BIN = Path("/home/ferris/sc-online/scripts/micromamba_root/envs/osca-venv/bin")
if not ENV_BIN.is_dir():
    sys.exit("ERROR: osca-venv not found. Run setup.sh first.")
if sys.executable != str(ENV_BIN / "python"):
    os.execv(str(ENV_BIN / "python"), [str(ENV_BIN / "python")] + sys.argv)

# ── Standard imports (env guaranteed from here) ────────────────────────────────
import argparse
import traceback

import pandas as pd
import pyreadr  # type: ignore

sys.path.insert(0, str(Path(__file__).parent))
from validation import (
    qq_lambda,
    pi1_replication,
    effect_size_maf,
    tss_enrichment,
    atac_enrichment,
    egene_chrom,
    go_enrichment,
    cis_distance,
    eqtl_per_gene,
    top_snp_fraction,
)

# v2 RDS column schema (from step 5):
# SNP, Chr, BP, A1, A2, Freq, ensg_id, gene_symbol, Probe_Chr, Probe_bp,
# Orientation, b, SE, p, p_bonf, padj_gene, padj_snp,
# is_significant_gene, is_significant_snp


def _load_rds(rds_path: Path) -> pd.DataFrame:
    result = pyreadr.read_r(str(rds_path))
    df = result[None] if None in result else next(iter(result.values()))
    if not isinstance(df, pd.DataFrame):
        raise ValueError(f"Expected DataFrame from {rds_path}, got {type(df)}")
    return df


def _load_gene_loc(gene_loc_path: Path) -> pd.DataFrame:
    """
    Load gene_loc_v2.txt.
    Columns (tab-delimited, with header): ensg_id, chr, TSS, gene_symbol, strand
    Returns DataFrame with these standardized column names.
    """
    df = pd.read_csv(gene_loc_path, sep="\t",
                     dtype={"ensg_id": str, "chr": str, "TSS": int,
                            "gene_symbol": str, "strand": str})
    return df


def _make_mean_expr_csv(pb_output_dir: Path, cc: str, cc_out: Path) -> Path | None:
    """
    Discover expression CSV for cell class `cc`, compute per-ENSG mean expression,
    and write {cc_out}/mean_expr.csv with columns gene_id, mean_expr.
    Returns path to the CSV, or None if not found.
    """
    expr_csv_path = pb_output_dir / f"{cc}_expression_matrix_ds.csv"
    if not expr_csv_path.exists():
        print(f"  [step9] Expression CSV not found: {expr_csv_path} — skipping expr QQ.")
        return None

    print(f"  [step9] Computing mean expression from {expr_csv_path} ...")
    # CSV: rows = participants, cols = participant_id + ENSG IDs
    expr = pd.read_csv(expr_csv_path, index_col=0)
    mean_expr = expr.mean(axis=0).reset_index()
    mean_expr.columns = ["gene_id", "mean_expr"]
    # keep only ENSG IDs
    mean_expr = mean_expr[mean_expr["gene_id"].str.startswith("ENSG")]
    out_path = cc_out / "mean_expr.csv"
    mean_expr.to_csv(out_path, index=False)
    print(f"  [step9] Mean expression saved → {out_path} ({len(mean_expr)} genes)")
    return out_path


def _run_all(cc: str, rds_path: Path, out_dir: Path, gene_loc_df: pd.DataFrame,
             args: argparse.Namespace) -> None:
    print(f"\n{'='*60}")
    print(f" Cell class: {cc}")
    print(f"{'='*60}")
    cc_out = out_dir / cc
    cc_out.mkdir(parents=True, exist_ok=True)

    print(f"  Loading {rds_path} …")
    df = _load_rds(rds_path)
    print(f"  Loaded {len(df):,} rows × {len(df.columns)} columns")

    # Auto-discover expression CSV for expression-stratified QQ
    expr_csv: Path | None = None
    if args.pb_output_dir is not None:
        pb_dir = Path(args.pb_output_dir)
        expr_csv = _make_mean_expr_csv(pb_dir, cc, cc_out)

    shared_kwargs = dict(
        gene_loc_df=gene_loc_df,
        padj_thresh=args.padj_thresh,
        gtex_sn=args.gtex_sn,
        atac_bed=args.atac_bed,
        go_bp_gmt=args.go_bp_gmt,
        go_mf_gmt=args.go_mf_gmt,
        expr_csv=expr_csv,
    )

    analyses = [
        # Fast analyses first (no external resources)
        ("qq_lambda",        qq_lambda),
        ("effect_size_maf",  effect_size_maf),
        ("egene_chrom",      egene_chrom),
        ("cis_distance",     cis_distance),
        ("eqtl_per_gene",    eqtl_per_gene),
        ("top_snp_fraction", top_snp_fraction),
        ("tss_enrichment",   tss_enrichment),
        ("go_enrichment",    go_enrichment),
        # Slow external-resource analyses last
        ("pi1_replication",  pi1_replication),
        ("atac_enrichment",  atac_enrichment),
    ]

    for name, mod in analyses:
        print(f"\n  ── {name} ──")
        try:
            mod.run(df, cc_out, **shared_kwargs)
        except Exception:
            print(f"  ERROR in {name}:")
            traceback.print_exc()


def main() -> None:
    ap = argparse.ArgumentParser(description="eQTL validation pipeline (step 9)")
    ap.add_argument("--eqtl-dir",      required=True,  help="Dir with eqtl_{cc}.rds files (step 5 output)")
    ap.add_argument("--gene-loc",      required=True,  help="gene_loc_v2.txt path (ensg_id-keyed)")
    ap.add_argument("--pb-output-dir", default=None,   help="Pseudobulk output dir; used to auto-discover per-CC expression CSV for QQ stratification")
    ap.add_argument("--gtex-sn",       default=None,   help="GTEx SN significant pairs .txt.gz")
    ap.add_argument("--atac-bed",      default=None,   help="Corces 2020 ATAC peaks .bed.gz")
    ap.add_argument("--go-bp-gmt",     default=None,   help="GO Biological Process GMT")
    ap.add_argument("--go-mf-gmt",     default=None,   help="GO Molecular Function GMT")
    ap.add_argument("--cell-class",    default=None,   help="Process only this cell class (default: all discovered)")
    ap.add_argument("--out-dir",       required=True,  help="Validation output root directory")
    ap.add_argument("--padj-thresh",   type=float, default=0.05)
    args = ap.parse_args()

    eqtl_dir = Path(args.eqtl_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover cell classes from eqtl_{cc}.rds (exclude _sig and present_in_all)
    rds_files = sorted([
        f for f in eqtl_dir.glob("eqtl_*.rds")
        if "_sig" not in f.name and "present_in_all" not in f.name
    ])
    if not rds_files:
        sys.exit(f"ERROR: No eqtl_*.rds files found in {eqtl_dir}")

    cell_classes = {}
    for f in rds_files:
        cc = f.stem.replace("eqtl_", "", 1)
        cell_classes[cc] = f

    if args.cell_class:
        if args.cell_class not in cell_classes:
            sys.exit(f"ERROR: --cell-class '{args.cell_class}' not found. Available: {list(cell_classes.keys())}")
        cell_classes = {args.cell_class: cell_classes[args.cell_class]}

    print(f"Found {len(cell_classes)} cell class(es): {list(cell_classes.keys())}")

    gene_loc_df = _load_gene_loc(Path(args.gene_loc))
    print(f"Loaded gene_loc_v2: {len(gene_loc_df)} genes")

    # Run analyses sequentially within this process.
    # Parallelism across cell classes is handled externally by the bash orchestrator,
    # which launches one process per CC (each with --cell-class).
    for cc, rds_path in cell_classes.items():
        _run_all(cc, rds_path, out_dir, gene_loc_df, args)

    print(f"\n{'='*60}")
    print(" Step 9 eQTL validation complete.")
    print(f" Outputs → {out_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
