"""
qq_lambda.py
------------
QQ plot and genomic inflation factor (lambda_GC) stratified by:
  1. MAF bin (3 bins)
  2. Expression quartile (4 quartiles, optional -- requires expr_csv kwarg)

Uses 1000 quantiles so the CSV/plot is always compact regardless of input size.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2

N_QUANTILES = 1000


def _quantile_qq(p_bin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (expected, observed) -log10 p-values at N_QUANTILES quantile points.

    q goes 0 -> 1, so:
      np.quantile(p, q)  = smallest -> largest p   (obs: large -> 0 in -log10)
      -log10(q)          = largest -> 0             (exp: large -> 0 in -log10)
    Both axes decrease together -> positive slope, diagonal when p ~ Uniform(0,1).
    Bug note: using -log10(1-q) would invert the expected axis, giving negative slope.
    """
    q = np.linspace(0, 1, N_QUANTILES + 2)[1:-1]   # interior quantiles, never 0 or 1
    obs = -np.log10(np.clip(np.quantile(p_bin, q), 1e-300, 1.0))
    exp = -np.log10(np.clip(q, 1e-300, 1.0))        # q-th quantile of Uniform(0,1) is q
    return exp, obs


def _lambda_gc(p: np.ndarray) -> float:
    return float(np.median(chi2.ppf(1 - np.clip(p, 1e-300, 1.0), 1)) / chi2.ppf(0.5, 1))


def _qq_panel(ax, df_plot: pd.DataFrame, group_col: str, groups: list,
              colors: list, lambdas: dict) -> None:
    """Scatter quantile QQ onto ax, one series per group."""
    max_val = 0.0
    for grp, color in zip(groups, colors):
        sub = df_plot[df_plot[group_col] == grp]
        if sub.empty:
            continue
        lam = lambdas.get(grp, np.nan)
        lam_str = f"lambda={lam:.3f}" if not np.isnan(lam) else "lambda=N/A"
        ax.scatter(sub["expected_neg_log10_p"], sub["observed_neg_log10_p"],
                   s=8, alpha=0.6, color=color, label=f"{grp} ({lam_str})", rasterized=True)
        max_val = max(max_val, sub["expected_neg_log10_p"].max())
    ax.plot([0, max_val or 8], [0, max_val or 8], "k--", lw=1)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")


def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None, expr_csv=None, **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── MAF-stratified QQ ─────────────────────────────────────────────────────
    maf = df["Freq"].clip(0, 1).apply(lambda x: min(x, 1 - x))

    maf_bins   = [(0, 0.1), (0.1, 0.25), (0.25, 0.5)]
    maf_labels = ["MAF<0.10", "0.10<=MAF<0.25", "0.25<=MAF<0.50"]
    maf_colors = sns.color_palette("tab10", len(maf_bins))

    maf_records = []
    maf_lambdas = {}
    for (lo, hi), label in zip(maf_bins, maf_labels):
        idx = df.index[df["Freq"].notna() & (maf >= lo) & (maf < hi)]
        p_bin = df.loc[idx, "p"].dropna().values
        if len(p_bin) < 10:
            maf_lambdas[label] = np.nan
            continue
        maf_lambdas[label] = _lambda_gc(p_bin)
        exp, obs = _quantile_qq(p_bin)
        for e, o in zip(exp, obs):
            maf_records.append({"maf_bin": label,
                                 "expected_neg_log10_p": e,
                                 "observed_neg_log10_p": o})

    pd.DataFrame(maf_records).to_csv(out_dir / "qq_by_maf_bin.csv", index=False)
    df_maf = pd.DataFrame(maf_records)

    sns.set_style("whitegrid")

    # Combined MAF QQ
    fig, ax = plt.subplots(figsize=(10, 6))
    _qq_panel(ax, df_maf, "maf_bin", maf_labels, maf_colors, maf_lambdas)
    ax.set_title(f"QQ plot by MAF bin  (n={N_QUANTILES} quantiles)")
    ax.legend(markerscale=2)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"qq_by_maf_bin.{ext}", dpi=150)
    plt.close(fig)

    # 4-panel MAF QQ
    fig, axes = plt.subplots(1, len(maf_bins), figsize=(14, 5))
    for ax, (lo, hi), label, color in zip(axes, maf_bins, maf_labels, maf_colors):
        sub = df_maf[df_maf["maf_bin"] == label]
        lam = maf_lambdas.get(label, np.nan)
        lam_str = f"lambda={lam:.3f}" if not np.isnan(lam) else "lambda=N/A"
        if not sub.empty:
            ax.scatter(sub["expected_neg_log10_p"], sub["observed_neg_log10_p"],
                       s=8, alpha=0.6, color=color, rasterized=True)
            mv = sub["expected_neg_log10_p"].max()
            ax.plot([0, mv], [0, mv], "k--", lw=1)
        ax.set_title(f"{label}\n{lam_str}")
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"qq_by_maf_bin_4panel.{ext}", dpi=150)
    plt.close(fig)

    print(f"  [qq_lambda] MAF lambda_GC: {maf_lambdas}")

    # ── Expression-stratified QQ (optional) ───────────────────────────────────
    if expr_csv is None or not Path(str(expr_csv)).exists():
        print(f"  [qq_lambda] No expr_csv provided -- skipping expression QQ.")
        print(f"  [qq_lambda] Outputs -> {out_dir}")
        return

    print(f"  [qq_lambda] Loading expression from {expr_csv} ...")
    expr_df = pd.read_csv(expr_csv)
    if "gene_id" not in expr_df.columns or "mean_expr" not in expr_df.columns:
        raise ValueError("expr_csv must have columns: gene_id, mean_expr")
    expr_map = dict(zip(expr_df["gene_id"].astype(str), expr_df["mean_expr"]))

    # Join on Probe (gene symbol) — Gene column is Entrez int, Probe is symbol
    probe_col = "Probe" if "Probe" in df.columns else "Gene"
    gene_expr = df[probe_col].astype(str).map(expr_map)
    n_mapped = gene_expr.notna().sum()
    print(f"  [qq_lambda] Mapped expression for {n_mapped:,}/{len(df):,} rows.")
    if n_mapped < 100:
        print("  [qq_lambda] Too few genes with expression -- skipping expression QQ.")
        print(f"  [qq_lambda] Outputs -> {out_dir}")
        return

    gene_expr_q = pd.qcut(gene_expr, 4, labels=["Q1 (low)", "Q2", "Q3", "Q4 (high)"])
    expr_labels = ["Q1 (low)", "Q2", "Q3", "Q4 (high)"]
    expr_colors = sns.color_palette("coolwarm", 4)

    expr_records = []
    expr_lambdas = {}
    for label in expr_labels:
        idx = df.index[gene_expr_q == label]
        p_bin = df.loc[idx, "p"].dropna().values
        if len(p_bin) < 10:
            expr_lambdas[label] = np.nan
            continue
        expr_lambdas[label] = _lambda_gc(p_bin)
        exp, obs = _quantile_qq(p_bin)
        for e, o in zip(exp, obs):
            expr_records.append({"expr_bin": label,
                                  "expected_neg_log10_p": e,
                                  "observed_neg_log10_p": o})

    pd.DataFrame(expr_records).to_csv(out_dir / "qq_by_expr_quartile.csv", index=False)
    df_expr = pd.DataFrame(expr_records)

    # Combined expression QQ
    fig, ax = plt.subplots(figsize=(10, 6))
    _qq_panel(ax, df_expr, "expr_bin", expr_labels, expr_colors, expr_lambdas)
    ax.set_title(f"QQ plot by expression quartile  (n={N_QUANTILES} quantiles)")
    ax.legend(markerscale=2)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"qq_by_expr_quartile.{ext}", dpi=150)
    plt.close(fig)

    # 4-panel expression QQ
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    for ax, label, color in zip(axes, expr_labels, expr_colors):
        sub = df_expr[df_expr["expr_bin"] == label]
        lam = expr_lambdas.get(label, np.nan)
        lam_str = f"lambda={lam:.3f}" if not np.isnan(lam) else "lambda=N/A"
        if not sub.empty:
            ax.scatter(sub["expected_neg_log10_p"], sub["observed_neg_log10_p"],
                       s=8, alpha=0.6, color=color, rasterized=True)
            mv = sub["expected_neg_log10_p"].max()
            ax.plot([0, mv], [0, mv], "k--", lw=1)
        ax.set_title(f"{label}\n{lam_str}")
        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"qq_by_expr_quartile_4panel.{ext}", dpi=150)
    plt.close(fig)

    print(f"  [qq_lambda] Expression lambda_GC: {expr_lambdas}")
    print(f"  [qq_lambda] Outputs -> {out_dir}")
