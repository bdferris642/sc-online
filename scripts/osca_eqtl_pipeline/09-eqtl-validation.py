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

Usage:
  python 09-eqtl-validation.py \\
    --eqtl-dir  /path/to/eqtl_final_outs/my_run \\
    --gene-loc  /path/to/gene_loc_new.txt \\
    --out-dir   /path/to/validation \\
    [--gtex-sn  /path/to/gtex_sn_signif_pairs.txt.gz] \\
    [--atac-bed /path/to/corces_2020_da_atac_peaks.bed.gz] \\
    [--go-bp-gmt /path/to/GO_Biological_Process_2025.gmt] \\
    [--go-mf-gmt /path/to/GO_Molecular_Function_2025.gmt] \\
    [--padj-thresh 0.05] \\
    [--threads 4]
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
from concurrent.futures import ProcessPoolExecutor, as_completed

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

# RDS column schema (from step 5):
# SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr, Probe_bp,
# Gene, Orientation, b, SE, p, p_bonf, padj_gene, padj_snp,
# is_significant_gene, is_significant_snp


def _load_rds(rds_path: Path) -> pd.DataFrame:
    result = pyreadr.read_r(str(rds_path))
    # pyreadr returns OrderedDict; None key = unnamed R object
    df = result[None] if None in result else next(iter(result.values()))
    if not isinstance(df, pd.DataFrame):
        raise ValueError(f"Expected DataFrame from {rds_path}, got {type(df)}")
    return df


def _load_gene_loc(gene_loc_path: Path) -> pd.DataFrame:
    """
    Load gene_loc_new.txt.
    Actual columns (tab-delimited, with header): probe, chr, TSS, NAME, strand
    where probe = Entrez/probe ID, NAME = gene symbol.
    Returns DataFrame with standardized names: entrez, chr, TSS, symbol, strand.
    """
    df = pd.read_csv(gene_loc_path, sep="\t",
                     dtype={"probe": str, "chr": str, "TSS": int, "NAME": str, "strand": str})
    df = df.rename(columns={"probe": "entrez", "NAME": "symbol"})
    return df


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

    shared_kwargs = dict(
        gene_loc_df=gene_loc_df,
        padj_thresh=args.padj_thresh,
        gtex_sn=args.gtex_sn,
        atac_bed=args.atac_bed,
        go_bp_gmt=args.go_bp_gmt,
        go_mf_gmt=args.go_mf_gmt,
        expr_csv=args.expr_csv,
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
    ap.add_argument("--gene-loc",      required=True,  help="gene_loc_new.txt path")
    ap.add_argument("--gtex-sn",       default=None,   help="GTEx SN significant pairs .txt.gz")
    ap.add_argument("--atac-bed",      default=None,   help="Corces 2020 ATAC peaks .bed.gz")
    ap.add_argument("--go-bp-gmt",     default=None,   help="GO Biological Process GMT")
    ap.add_argument("--go-mf-gmt",     default=None,   help="GO Molecular Function GMT")
    ap.add_argument("--expr-csv",      default=None,   help="CSV with gene_id,mean_expr for expression-stratified QQ")
    ap.add_argument("--cell-class",    default=None,   help="Process only this cell class (default: all discovered)")
    ap.add_argument("--out-dir",       required=True,  help="Validation output root directory")
    ap.add_argument("--padj-thresh",   type=float, default=0.05)
    ap.add_argument("--threads",       type=int,   default=4,
                    help="Parallel workers (one per cell class)")
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
    print(f"Loaded gene_loc: {len(gene_loc_df)} genes")

    # Run analyses in parallel (one worker per cell class)
    n_workers = min(args.threads, len(cell_classes))
    if n_workers <= 1:
        for cc, rds_path in cell_classes.items():
            _run_all(cc, rds_path, out_dir, gene_loc_df, args)
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as exe:
            futs = {
                exe.submit(_run_all, cc, rds_path, out_dir, gene_loc_df, args): cc
                for cc, rds_path in cell_classes.items()
            }
            for fut in as_completed(futs):
                cc = futs[fut]
                try:
                    fut.result()
                except Exception:
                    print(f"ERROR processing cell class {cc}:")
                    traceback.print_exc()

    print(f"\n{'='*60}")
    print(" Step 9 eQTL validation complete.")
    print(f" Outputs → {out_dir}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
