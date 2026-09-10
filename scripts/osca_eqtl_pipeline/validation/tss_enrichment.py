"""
tss_enrichment.py
-----------------
Enrichment of significant eQTL SNPs near TSS relative to background.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde


def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None, padj_thresh: float = 0.05, **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Join Probe_bp from gene_loc_v2 (ensg_id-keyed) if needed; fall back to df["Probe_bp"] directly
    if "Probe_bp" in df.columns:
        probe_bp = df["Probe_bp"]
    elif gene_loc_df is not None:
        # gene_loc_v2 cols: ensg_id, chr, TSS, gene_symbol, strand
        tss_map = dict(zip(gene_loc_df["ensg_id"].astype(str), gene_loc_df["TSS"]))
        probe_bp = df["ensg_id"].astype(str).map(tss_map)
    else:
        print("  [tss_enrichment] No Probe_bp or gene_loc — skipping.")
        return

    bp = pd.to_numeric(df["BP"], errors="coerce")
    tss = pd.to_numeric(probe_bp, errors="coerce")
    dist = (bp - tss).abs()

    padj_snp = df.get("padj_snp", pd.Series(np.nan, index=df.index))
    sig_mask = padj_snp < padj_thresh

    dist_sig = dist[sig_mask].dropna().values
    dist_all = dist.dropna().values
    # Downsample background to 1M
    if len(dist_all) > 1_000_000:
        rng = np.random.default_rng(42)
        dist_all = rng.choice(dist_all, 1_000_000, replace=False)

    dist_bins = [0, 10_000, 100_000, 500_000, 1_000_000]
    bin_labels = ["<10kb", "10–100kb", "100–500kb", "500kb–1Mb"]

    # Per-bin fold enrichment
    n_sig_total = len(dist_sig)
    n_all_total = len(dist_all)
    records = []
    for lo, hi, label in zip(dist_bins[:-1], dist_bins[1:], bin_labels):
        n_sig_bin = np.sum((dist_sig >= lo) & (dist_sig < hi))
        n_all_bin = np.sum((dist_all >= lo) & (dist_all < hi))
        frac_sig = n_sig_bin / n_sig_total if n_sig_total > 0 else 0
        frac_all = n_all_bin / n_all_total if n_all_total > 0 else 0
        fe = frac_sig / frac_all if frac_all > 0 else np.nan
        records.append({"dist_bin": label, "n_sig": n_sig_bin, "n_background": n_all_bin,
                        "frac_sig": frac_sig, "frac_background": frac_all, "fold_enrichment": fe})

    csv_df = pd.DataFrame(records)
    csv_df.to_csv(out_dir / "tss_enrichment.csv", index=False)

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram + KDE on log10 distance
    ax = axes[0]
    log_sig = np.log10(dist_sig + 1)
    log_all = np.log10(dist_all + 1)
    ax.hist(log_all, bins=60, density=True, alpha=0.4, color="steelblue", label="All tested")
    ax.hist(log_sig, bins=60, density=True, alpha=0.6, color="darkorange", label="Significant")
    for data, color in [(log_all, "steelblue"), (log_sig, "darkorange")]:
        if len(data) > 1:
            kde = gaussian_kde(data, bw_method=0.3)
            x = np.linspace(data.min(), data.max(), 200)
            ax.plot(x, kde(x), color=color, lw=2)
    ax.set_xlabel("log₁₀(|SNP BP – TSS| + 1)")
    ax.set_ylabel("Density")
    ax.set_title("TSS distance distribution")
    ax.legend()

    # Fold enrichment bar chart
    ax = axes[1]
    fe_vals = csv_df["fold_enrichment"].fillna(0).values
    ax.bar(bin_labels, fe_vals, color="teal", edgecolor="black")
    ax.axhline(1, color="red", linestyle="--", lw=1)
    ax.set_xlabel("Distance to TSS bin")
    ax.set_ylabel("Fold enrichment (sig / background)")
    ax.set_title("Enrichment of significant SNPs by TSS distance")

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"tss_enrichment.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [tss_enrichment] Outputs → {out_dir}")
