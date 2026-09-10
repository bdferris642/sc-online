"""
eqtl_per_gene.py
----------------
Distribution of the number of significant eSNPs per eGene.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None, padj_thresh: float = 0.05, **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    padj_snp = df.get("padj_snp", pd.Series(np.nan, index=df.index))
    sig_df = df[padj_snp < padj_thresh]

    counts = sig_df.groupby("ensg_id").size().reset_index(name="n_sig_snps")
    # add gene_symbol column if available
    if "gene_symbol" in sig_df.columns:
        sym_map = sig_df.drop_duplicates("ensg_id").set_index("ensg_id")["gene_symbol"]
        counts["gene_symbol"] = counts["ensg_id"].map(sym_map)
    counts.to_csv(out_dir / "eqtls_per_gene.csv", index=False)

    if counts.empty:
        print("  [eqtl_per_gene] No significant SNPs found.")
        return

    med = counts["n_sig_snps"].median()
    p90 = counts["n_sig_snps"].quantile(0.9)
    mx = counts["n_sig_snps"].max()

    log_counts = np.log10(counts["n_sig_snps"].values + 1)

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(log_counts, bins=40, color="steelblue", edgecolor="black", alpha=0.8)
    ax.set_xlabel("log₁₀(# significant eSNPs per eGene + 1)")
    ax.set_ylabel("Count (genes)")
    ax.set_title(f"eQTLs per eGene\nMedian={med:.0f}, 90th pctile={p90:.0f}, Max={mx}")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"eqtls_per_gene.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [eqtl_per_gene] {len(counts)} eGenes, median={med:.0f} eSNPs → {out_dir}")
