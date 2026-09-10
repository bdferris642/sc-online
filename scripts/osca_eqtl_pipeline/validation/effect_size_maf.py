"""
effect_size_maf.py
------------------
Effect size |β| vs MAF: hexbin scatter + binned means for all tested vs significant.
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

    maf = df["Freq"].clip(0, 1).apply(lambda x: min(x, 1 - x))
    abs_b = df["b"].abs()

    df_work = pd.DataFrame({"maf": maf, "abs_b": abs_b, "padj_snp": df.get("padj_snp", np.nan)}).dropna(subset=["maf", "abs_b"])
    df_sig = df_work[df_work["padj_snp"] < padj_thresh] if "padj_snp" in df_work.columns else pd.DataFrame()

    # MAF bins (10 bins 0–0.5)
    bin_edges = np.linspace(0, 0.5, 11)
    bin_labels = [f"{bin_edges[i]:.2f}–{bin_edges[i+1]:.2f}" for i in range(10)]

    def bin_stats(sub):
        cats = pd.cut(sub["maf"], bins=bin_edges, labels=bin_labels, include_lowest=True)
        return sub.groupby(cats, observed=True)["abs_b"].agg(["mean", "sem", "count"]).reset_index()

    stats_all = bin_stats(df_work)
    stats_all["group"] = "all"
    stats_sig = bin_stats(df_sig) if not df_sig.empty else pd.DataFrame()
    if not stats_sig.empty:
        stats_sig["group"] = "significant"

    csv_df = pd.concat([stats_all, stats_sig], ignore_index=True)
    csv_df.to_csv(out_dir / "effect_size_vs_maf.csv", index=False)

    # Downsample to 500k for hexbin
    n_max = 500_000
    sub_plot = df_work.sample(min(n_max, len(df_work)), random_state=42)

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Hexbin
    ax = axes[0]
    hb = ax.hexbin(sub_plot["maf"], sub_plot["abs_b"], gridsize=60, cmap="YlOrRd", mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax, label="log10(count)")
    ax.set_xlabel("MAF")
    ax.set_ylabel("|β|")
    ax.set_title("Effect size vs MAF (hexbin)")

    # Binned means
    ax = axes[1]
    bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(10)]
    ax.errorbar(bin_centers, stats_all["mean"], yerr=stats_all["sem"], fmt="o-", label="All tested", color="steelblue")
    if not stats_sig.empty:
        ax.errorbar(bin_centers, stats_sig["mean"], yerr=stats_sig["sem"], fmt="s--", label="Significant", color="darkorange")
    ax.set_xlabel("MAF bin (center)")
    ax.set_ylabel("Mean |β| ± SE")
    ax.set_title("Binned mean effect size vs MAF")
    ax.legend()

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"effect_size_vs_maf.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [effect_size_maf] Outputs → {out_dir}")
