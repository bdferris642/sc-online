"""
cis_distance.py
---------------
Distribution of cis distance for lead eQTL SNPs (min-p per eGene).
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

    # Probe_bp from df or gene_loc
    if "Probe_bp" in df.columns:
        probe_bp = df["Probe_bp"]
    elif gene_loc_df is not None:
        # gene_loc named cols: entrez, chr, TSS, symbol, strand
        tss_map = dict(zip(gene_loc_df["entrez"].astype(str), gene_loc_df["TSS"]))
        probe_bp = df["Gene"].astype(str).map(tss_map)
    else:
        print("  [cis_distance] No Probe_bp — skipping.")
        return

    df_work = df.copy()
    df_work["BP"] = pd.to_numeric(df_work["BP"], errors="coerce")
    df_work["Probe_bp_val"] = pd.to_numeric(probe_bp, errors="coerce")
    df_work["cis_dist"] = (df_work["BP"] - df_work["Probe_bp_val"]).abs()
    df_work = df_work.dropna(subset=["cis_dist"])

    # Lead SNP per eGene = min p per gene
    padj_gene = df_work.get("padj_gene", pd.Series(np.nan, index=df_work.index))
    sig_genes = df_work[padj_gene < padj_thresh]["Gene"].unique()
    df_sig = df_work[df_work["Gene"].isin(sig_genes)]

    lead_idx = df_sig.groupby("Gene")["p"].idxmin()
    lead_df = df_sig.loc[lead_idx].copy()

    csv_df = lead_df[["Gene", "SNP", "BP", "Probe_bp_val", "cis_dist", "p"]].reset_index(drop=True)
    csv_df.to_csv(out_dir / "cis_distance.csv", index=False)

    dist_sig = lead_df["cis_dist"].values
    dist_all = df_work["cis_dist"].values
    if len(dist_all) > 1_000_000:
        rng = np.random.default_rng(42)
        dist_all = rng.choice(dist_all, 1_000_000, replace=False)

    med = np.median(dist_sig) if len(dist_sig) > 0 else np.nan
    p90 = np.percentile(dist_sig, 90) if len(dist_sig) > 0 else np.nan

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    log_sig = np.log10(dist_sig + 1)
    log_all = np.log10(dist_all + 1)
    ax.hist(log_all, bins=60, density=True, alpha=0.4, color="steelblue", label="All SNP-gene pairs")
    ax.hist(log_sig, bins=60, density=True, alpha=0.6, color="darkorange", label="Lead eQTL SNPs")
    for data, color in [(log_all, "steelblue"), (log_sig, "darkorange")]:
        if len(data) > 1:
            kde = gaussian_kde(data, bw_method=0.3)
            x = np.linspace(data.min(), data.max(), 200)
            ax.plot(x, kde(x), color=color, lw=2)
    ax.set_xlabel("log₁₀(cis distance + 1) [bp]")
    ax.set_ylabel("Density")
    ax.set_title(f"Cis distance of lead eQTL SNPs\nMedian={med:.0f} bp, 90th pctile={p90:.0f} bp")
    ax.legend()
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"cis_distance.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [cis_distance] {len(dist_sig)} lead SNPs, median={med:.0f} bp → {out_dir}")
