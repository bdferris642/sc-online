"""
egene_chrom.py
--------------
eGenes per chromosome: raw count + normalized (sig/total tested per chrom).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


# Chromosomes to annotate with labels
CHROM_LABELS = {"6": "chr6\n(MHC)", "17": "chr17\n(MAPT/H1H2)"}


def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None, padj_thresh: float = 0.05, **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_work = df.copy()
    df_work["Chr"] = df_work["Chr"].astype(str)
    df_work["Gene"] = df_work["Gene"].astype(str)

    # Total tested genes per chrom (one entry per unique Gene-Chrom pair)
    per_gene = df_work.drop_duplicates(subset=["Gene", "Chr"])
    total_per_chrom = per_gene.groupby("Chr")["Gene"].nunique()

    # Significant eGenes per chrom
    padj_gene = df_work.get("padj_gene", pd.Series(np.nan, index=df_work.index))
    sig_genes_df = df_work[padj_gene < padj_thresh].drop_duplicates(subset=["Gene", "Chr"])
    sig_per_chrom = sig_genes_df.groupby("Chr")["Gene"].nunique()

    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    chroms = [c for c in chrom_order if c in total_per_chrom.index]

    csv_df = pd.DataFrame({
        "chr": chroms,
        "n_sig_egenes": [sig_per_chrom.get(c, 0) for c in chroms],
        "n_total_genes": [total_per_chrom.get(c, 0) for c in chroms],
    })
    csv_df["frac_sig"] = csv_df["n_sig_egenes"] / csv_df["n_total_genes"].replace(0, np.nan)
    csv_df.to_csv(out_dir / "egene_per_chrom.csv", index=False)

    xlab = [CHROM_LABELS.get(c, c) for c in chroms]

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    ax = axes[0]
    ax.bar(range(len(chroms)), csv_df["n_sig_egenes"], color="steelblue", edgecolor="black")
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(xlab, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("# significant eGenes")
    ax.set_title("eGenes per chromosome (raw count)")
    for i, c in enumerate(chroms):
        if c in CHROM_LABELS:
            ax.axvline(i, color="red", linestyle="--", alpha=0.5, lw=1)

    ax = axes[1]
    ax.bar(range(len(chroms)), csv_df["frac_sig"].fillna(0), color="darkorange", edgecolor="black")
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(xlab, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Fraction (sig eGenes / total tested genes)")
    ax.set_title("eGenes per chromosome (normalized)")
    for i, c in enumerate(chroms):
        if c in CHROM_LABELS:
            ax.axvline(i, color="red", linestyle="--", alpha=0.5, lw=1)

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"egene_per_chrom.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [egene_chrom] Outputs → {out_dir}")
