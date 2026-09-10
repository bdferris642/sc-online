"""
atac_enrichment.py
------------------
Overlap of significant eQTL SNPs with Corces 2020 DA neuron ATAC-seq peaks,
per brain cell type. Uses fast interval overlap via numpy.searchsorted.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact


def _load_peaks(atac_bed: Path) -> pd.DataFrame:
    """Load BED gz; return DataFrame with columns: chr, start, end, cell_type."""
    peaks = pd.read_csv(atac_bed, sep="\t", header=None, comment="#",
                        dtype={0: str, 1: int, 2: int})
    if peaks.shape[1] >= 4:
        peaks.columns = list(peaks.columns[:4]) + list(peaks.columns[4:])
        peaks = peaks.rename(columns={0: "chr", 1: "start", 2: "end", 3: "cell_type"})
    else:
        peaks.columns = ["chr", "start", "end"] + list(range(3, peaks.shape[1]))
        peaks["cell_type"] = "all"
    peaks["chr"] = peaks["chr"].str.replace("chr", "", regex=False)
    return peaks


def _overlap_fraction(snp_chr: np.ndarray, snp_pos: np.ndarray,
                      peak_chr: np.ndarray, peak_start: np.ndarray, peak_end: np.ndarray) -> float:
    """Fraction of SNPs overlapping any peak (binary, per-SNP)."""
    if len(snp_chr) == 0:
        return 0.0
    count = 0
    for chrom in np.unique(snp_chr):
        c_mask = peak_chr == chrom
        if not c_mask.any():
            continue
        ps = peak_start[c_mask]
        pe = peak_end[c_mask]
        order = np.argsort(ps)
        ps_s, pe_s = ps[order], pe[order]
        pos = snp_pos[snp_chr == chrom]
        idx = np.searchsorted(ps_s, pos, side="right") - 1
        in_peak = (idx >= 0) & (pos < pe_s[np.clip(idx, 0, len(pe_s) - 1)])
        count += int(in_peak.sum())
    return count / len(snp_chr)


def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None,
        atac_bed=None, padj_thresh: float = 0.05, **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if atac_bed is None or not Path(atac_bed).exists():
        print(f"  [atac_enrichment] ATAC BED not found — skipping.")
        return

    print(f"  [atac_enrichment] Loading ATAC peaks from {atac_bed} …")
    peaks = _load_peaks(Path(atac_bed))
    cell_types = peaks["cell_type"].unique().tolist()
    print(f"  [atac_enrichment] Cell types in BED: {cell_types}")

    # SNP positions
    df_snp = df[["Chr", "BP"]].copy()
    df_snp["Chr"] = df_snp["Chr"].astype(str)
    df_snp["BP"] = pd.to_numeric(df_snp["BP"], errors="coerce")
    df_snp = df_snp.dropna(subset=["BP"])

    padj_snp = df.get("padj_snp", pd.Series(np.nan, index=df.index))
    sig_mask = (padj_snp < padj_thresh).reindex(df_snp.index).fillna(False)

    # Downsample background to 200k
    bg_idx = df_snp.index
    if len(bg_idx) > 200_000:
        rng = np.random.default_rng(42)
        bg_idx = rng.choice(bg_idx, 200_000, replace=False)

    snp_chr_sig = df_snp.loc[sig_mask, "Chr"].values
    snp_pos_sig = df_snp.loc[sig_mask, "BP"].values.astype(int)
    snp_chr_bg = df_snp.loc[bg_idx, "Chr"].values
    snp_pos_bg = df_snp.loc[bg_idx, "BP"].values.astype(int)

    n_sig = len(snp_chr_sig)
    n_bg = len(snp_chr_bg)

    records = []
    for ct in cell_types:
        ct_peaks = peaks[peaks["cell_type"] == ct]
        pc = ct_peaks["chr"].values
        ps = ct_peaks["start"].values.astype(int)
        pe = ct_peaks["end"].values.astype(int)

        frac_sig = _overlap_fraction(snp_chr_sig, snp_pos_sig, pc, ps, pe)
        frac_bg = _overlap_fraction(snp_chr_bg, snp_pos_bg, pc, ps, pe)

        a = int(frac_sig * n_sig)
        b = n_sig - a
        c = int(frac_bg * n_bg)
        d = n_bg - c
        or_val, pval = fisher_exact([[a, b], [c, d]])
        # 95% CI on OR (Woolf's method)
        if a > 0 and b > 0 and c > 0 and d > 0:
            log_or = np.log(or_val)
            se_log = np.sqrt(1/a + 1/b + 1/c + 1/d)
            ci_lo = np.exp(log_or - 1.96 * se_log)
            ci_hi = np.exp(log_or + 1.96 * se_log)
        else:
            ci_lo, ci_hi = np.nan, np.nan
        records.append({"cell_type": ct, "n_sig_overlap": a, "n_sig_total": n_sig,
                        "n_bg_overlap": c, "n_bg_total": n_bg,
                        "frac_sig": frac_sig, "frac_background": frac_bg,
                        "OR": or_val, "pval": pval, "CI_lo": ci_lo, "CI_hi": ci_hi})

    csv_df = pd.DataFrame(records)
    csv_df.to_csv(out_dir / "atac_enrichment.csv", index=False)

    # Multi-bar chart: OR per cell type
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(max(8, len(cell_types) * 1.5), 5))
    x = np.arange(len(records))
    ors = csv_df["OR"].values
    ci_lo = csv_df["CI_lo"].values
    ci_hi = csv_df["CI_hi"].values
    ax.bar(x, ors, color="mediumseagreen", edgecolor="black", zorder=3)
    yerr_lo = np.where(np.isnan(ci_lo), 0, ors - ci_lo)
    yerr_hi = np.where(np.isnan(ci_hi), 0, ci_hi - ors)
    ax.errorbar(x, ors, yerr=[yerr_lo, yerr_hi], fmt="none", color="black", capsize=4, zorder=4)
    ax.axhline(1, color="red", linestyle="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(csv_df["cell_type"].values, rotation=30, ha="right")
    ax.set_ylabel("Odds Ratio (sig eQTL SNPs vs background)")
    ax.set_title("ATAC-seq peak overlap enrichment by brain cell type")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"atac_enrichment.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [atac_enrichment] Outputs → {out_dir}")
