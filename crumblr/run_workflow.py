#!/usr/bin/env python3
"""
run_workflow.py  –  crumblr proportion-change workflow

Aggregates cell-level data to sample-level counts, calls the R crumblr
analysis, then generates comprehensive Python-side plots.

Environment Setup (one-time)
-----------------------------
Run setup_env.sh once to install all dependencies into a self-contained
directory alongside this script:

    bash /path/to/crumblr/setup_env.sh

This creates:
  - A micromamba conda environment (crumblr) with R 4.4 + crumblr /
    Bioconductor packages, located at:
        <sandbox>/micromamba_root/envs/crumblr/
  - A Python virtual environment with numpy, pandas, matplotlib, seaborn,
    and scipy, located at:
        <sandbox>/crumblr-venv/

where <sandbox> is the parent directory of the crumblr/ folder
(i.e., the directory containing both crumblr/ and crumblr-venv/).

Running the Script
------------------
No manual environment activation is needed. The script automatically
re-execs itself with the crumblr-venv Python if it is not already running
inside it:

    python /path/to/crumblr/run_workflow.py --csv ... --output-dir results/

If the venv is not found, the script exits immediately with a clear message
pointing to setup_env.sh. To activate manually (e.g. for interactive use):

    source /path/to/cc-sandbox/crumblr-venv/bin/activate

The Rscript binary and R library path are resolved automatically relative
to this script's location — MAMBA_ROOT_PREFIX does not need to be set.

Path Resolution (in priority order)
-------------------------------------
Rscript binary:
  1. $MAMBA_ROOT_PREFIX/envs/crumblr/bin/Rscript  (if env var is set)
  2. <sandbox>/micromamba_root/envs/crumblr/bin/Rscript  (auto-detected)
  3. Rscript on PATH  (system fallback)

--r-script:
  Defaults to <dir of this script>/crumblr_analysis.R

--r-libs:
  Defaults to <sandbox>/R_libs if that directory exists (used when running
  with a system R instead of the conda env); ignored otherwise.

Usage
-----
python run_workflow.py \\
    --csv            cells.csv          \\
    --proportion-col cell_type          \\
    --sample-col     sample_id          \\
    --donor-col      donor_id           \\
    --contrast       case_control       \\
    --fixed-effects  age sex region_SN region_VTA dapi_nurr brain_bank study \\
    --random-effects donor_id           \\
    --output-dir     results/           \\
    [--ref-level     ctr]               \\
    [--r-libs        /path/to/R_libs]   \\
    [--r-script      /path/to/crumblr_analysis.R]
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# ── Self-activating: re-exec with crumblr-venv Python if not already using it ─
_VENV_PYTHON = Path(__file__).parent.parent / "crumblr-venv" / "bin" / "python3"
if Path(sys.executable).resolve() != _VENV_PYTHON.resolve():
    if not _VENV_PYTHON.is_file():
        sys.exit(
            f"ERROR: crumblr-venv not found at {_VENV_PYTHON.parent.parent}.\n"
            f"Run setup first:  bash {Path(__file__).parent / 'setup_env.sh'}"
        )
    os.execv(str(_VENV_PYTHON), [str(_VENV_PYTHON)] + sys.argv)
# ──────────────────────────────────────────────────────────────────────────────

import argparse
import json
import subprocess
import tempfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist

# ── Helpers ───────────────────────────────────────────────────────────────────

PALETTE = {
    "pd":  "#e74c3c",
    "ctr": "#2980b9",
    "PD":  "#e74c3c",
    "CTR": "#2980b9",
}


def _safe_name(s: str) -> str:
    return s.replace(" ", "_").replace("/", "_").replace("|", "_")


# ── Step 1: aggregate cells → sample-level counts ────────────────────────────

def aggregate_counts(
    df: pd.DataFrame,
    proportion_col: str,
    sample_col: str,
) -> pd.DataFrame:
    """Return a (samples × cell_types) count DataFrame."""
    counts = (
        df.groupby([sample_col, proportion_col])
          .size()
          .unstack(fill_value=0)
    )
    # Remove zero-count cell types
    counts = counts.loc[:, counts.sum() > 0]
    return counts


def build_metadata(
    df: pd.DataFrame,
    sample_col: str,
    covariate_cols: list[str],
) -> pd.DataFrame:
    """Return one row per sample with covariate values (first-occurrence).

    sample_col is excluded from covariate_cols because it is the groupby key
    and becomes the index automatically — this handles the common case where
    --sample-col, --donor-col, and a random-effect term all refer to the same
    column.
    """
    cols = [c for c in covariate_cols if c != sample_col]
    return df.groupby(sample_col)[cols].first()


# ── Step 1b: validate random effects ─────────────────────────────────────────

def validate_random_effects(
    metadata: pd.DataFrame,
    sample_col: str,
    random_effects: list[str],
) -> None:
    """Fail fast if any random effect has as many levels as observations.

    dream() / lme4 require strictly more observations than random-effect
    levels.  The most common mistake is specifying --random-effects with the
    same column used as --sample-col (or --donor-col), which produces one row
    per level after aggregation.
    """
    n_obs = len(metadata)
    for re_col in random_effects:
        if re_col in metadata.columns:
            n_levels = metadata[re_col].nunique()
        elif re_col == metadata.index.name:
            n_levels = metadata.index.nunique()
        else:
            continue  # R will report a missing-column error

        if n_levels >= n_obs:
            raise ValueError(
                f"Random effect '{re_col}' has {n_levels} unique levels for "
                f"{n_obs} observations — dream() requires n_levels < n_obs.\n\n"
                f"This happens when the column used for sample aggregation "
                f"('{sample_col}') is the same as a random-effect column, "
                f"leaving only one observation per level after groupby.\n\n"
                f"Options:\n"
                f"  1. Use a finer-grained sample column (e.g. participant+region)\n"
                f"     so multiple samples exist per '{re_col}' level.\n"
                f"  2. Remove '{re_col}' from --random-effects if each row is\n"
                f"     already an independent sample."
            )


# ── Step 2: relevel the contrast column ──────────────────────────────────────

def relevel(meta: pd.DataFrame, contrast: str, ref: str | None) -> pd.DataFrame:
    if ref is not None and ref in meta[contrast].values:
        cats = [ref] + [c for c in meta[contrast].unique() if c != ref]
        meta[contrast] = pd.Categorical(meta[contrast], categories=cats)
    return meta


# ── Step 3: build R formula string ───────────────────────────────────────────

def build_formula(
    contrast: str,
    fixed_effects: list[str],
    random_effects: list[str],
) -> str:
    fixed_terms = [contrast] + [f for f in fixed_effects if f != contrast]
    random_terms = [f"(1|{r})" for r in random_effects]
    rhs = " + ".join(fixed_terms + random_terms)
    return f"~ {rhs}"


# ── Step 4: call R ────────────────────────────────────────────────────────────

def run_r_analysis(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    formula: str,
    contrast: str,
    output_dir: Path,
    r_script: Path,
    r_libs: Path | None,
) -> Path:
    """Write temp CSVs, call crumblr_analysis.R, return output_dir."""
    with tempfile.TemporaryDirectory() as tmp:
        counts_csv   = Path(tmp) / "counts.csv"
        metadata_csv = Path(tmp) / "metadata.csv"

        # Align rows
        shared = counts.index.intersection(metadata.index)
        if len(shared) == 0:
            raise ValueError("No shared samples between counts and metadata.")
        counts   = counts.loc[shared]
        metadata = metadata.loc[shared]

        counts.to_csv(counts_csv)

        # Write metadata so R can use the sample column as both row.names
        # (for alignment with the counts matrix) and as a formula variable.
        # When sample_col == donor_col or a random-effect column, it lives
        # only in the pandas index after build_metadata(), making it invisible
        # to dream().  Fix: promote it to a data column and drop the index name
        # so R's read.csv(..., row.names=1) still picks it up as row names.
        meta_for_r = metadata.copy()
        idx_name = meta_for_r.index.name
        if idx_name and idx_name not in meta_for_r.columns:
            meta_for_r.insert(0, idx_name, meta_for_r.index)
        meta_for_r.index.name = None
        meta_for_r.to_csv(metadata_csv)

        env = os.environ.copy()
        if r_libs is not None and r_libs.is_dir():
            env["R_LIBS_CUSTOM"] = str(r_libs)

        # Resolve Rscript: env var → auto-detect relative to this script → system
        mamba_root = env.get("MAMBA_ROOT_PREFIX", "")
        if not mamba_root:
            mamba_root = str(Path(__file__).parent.parent / "micromamba_root")
        conda_rscript = Path(mamba_root) / "envs" / "crumblr" / "bin" / "Rscript"
        rscript_bin = str(conda_rscript) if conda_rscript.is_file() else "Rscript"

        cmd = [
            rscript_bin,
            "--no-save", "--no-restore", "--no-site-file", "--no-init-file",
            str(r_script),
            str(counts_csv),
            str(metadata_csv),
            formula,
            contrast,
            str(output_dir),
        ]
        print(f"\n[R] {rscript_bin} … formula='{formula}'")
        result = subprocess.run(cmd, env=env, capture_output=False)
        if result.returncode != 0:
            raise RuntimeError(
                f"crumblr_analysis.R exited with code {result.returncode}."
            )

    return output_dir


# ── Step 5: plots ─────────────────────────────────────────────────────────────

def _cell_class(ct: str) -> str:
    """Map a cell type name to its broad class.

    Handles both fine-grained subtypes (Da_SOX6_EYA4 → Dopaminergic) and the
    case where --proportion-col is already the class column (Dopaminergic →
    Dopaminergic), avoiding everything defaulting to Inhibitory/purple.
    """
    if ct in CLASS_COLORS:
        return ct
    if ct.startswith("Da_"):
        return "Dopaminergic"
    if ct.startswith("Exc_"):
        return "Excitatory"
    if ct.startswith("Inh_"):
        return "Inhibitory"
    if ct.startswith("Astro"):
        return "Astro"
    if ct.startswith("Oligo"):
        return "Oligo"
    if ct.startswith("Micro") or ct.startswith("Macro") or ct.startswith("Mono"):
        return "MG"
    if ct.startswith("Endo") or ct.startswith("Epen"):
        return "Endo"
    if ct.startswith("OPC"):
        return "OPC"
    if ct.startswith("Immune"):
        return "Immune"
    return "Other"


CLASS_COLORS = {
    "Dopaminergic": "#e67e22", 
    "Excitatory":   "#27ae60",
    "Inhibitory":   "#8e44ad",
    "Oligo":        "#39345e", 
    "Astro":        "#0e720a", 
    "MG":           "#880e00", 
    "Endo":         "#1aa6bc", 
    "OPC":          "#e84393", 
    "Immune":       "#facc12", 
    "Other":        "#7f8c8d",
}


def plot_composition_bar(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: str,
    output_dir: Path,
) -> None:
    """Stacked bar: mean proportions per contrast group."""
    props = counts.div(counts.sum(axis=1), axis=0)
    props[contrast] = metadata[contrast].astype(str)
    mean_props = props.groupby(contrast)[counts.columns.tolist()].mean()

    # Sort cell types by class then name
    ordered_cts = sorted(counts.columns, key=lambda c: (_cell_class(c), c))
    mean_props   = mean_props[ordered_cts]

    class_palette = {ct: CLASS_COLORS[_cell_class(ct)] for ct in ordered_cts}
    # Shade individual subtypes within class
    from matplotlib.colors import to_rgb, to_hex
    import colorsys

    def shade_variants(base_hex: str, n: int) -> list[str]:
        r, g, b = to_rgb(base_hex)
        h, s, v = colorsys.rgb_to_hsv(r, g, b)
        return [
            to_hex(colorsys.hsv_to_rgb(h, s * (0.5 + 0.5 * i / max(n - 1, 1)),
                                        v * (0.7 + 0.3 * i / max(n - 1, 1))))
            for i in range(n)
        ]

    colors = []
    from itertools import groupby
    for cls, group in groupby(ordered_cts, key=_cell_class):
        grp = list(group)
        shades = shade_variants(CLASS_COLORS[cls], len(grp))
        colors.extend(shades)

    fig, ax = plt.subplots(figsize=(7, 6))
    bottom = np.zeros(len(mean_props))
    for i, ct in enumerate(ordered_cts):
        vals = mean_props[ct].values
        ax.bar(mean_props.index, vals, bottom=bottom, color=colors[i],
               label=ct, width=0.55, edgecolor="white", linewidth=0.3)
        bottom += vals

    # Class legend patches
    patches = [mpatches.Patch(color=c, label=cls)
               for cls, c in CLASS_COLORS.items()]
    ax.legend(handles=patches, title="Cell class", bbox_to_anchor=(1.02, 1),
              loc="upper left", frameon=False)
    ax.set_ylabel("Mean proportion")
    ax.set_title(f"Cell type composition by {contrast}")
    ax.set_ylim(0, 1)
    plt.tight_layout()
    fig.savefig(output_dir / "composition_bar.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "composition_bar.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: composition_bar.{pdf,png}")


def plot_proportion_boxplots(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: str,
    results: pd.DataFrame,
    output_dir: Path,
    n_top: int = 12,
) -> None:
    """Strip + box plots for the top-N cell types, balanced across directions.

    Selects up to n_top // 2 enriched (logFC > 0) and n_top // 2 depleted
    (logFC < 0) by FDR, so the panel always represents both directions when
    hits exist in both.  Leftover budget is filled from whichever direction
    has more results.  Panels are ordered depleted → enriched by logFC.
    """
    props = counts.div(counts.sum(axis=1), axis=0)
    props[contrast] = metadata.loc[props.index, contrast].astype(str)
    groups = sorted(props[contrast].unique())

    ranked = (results.copy()
              .sort_values("adj.P.Val")
              .drop_duplicates("cell_type"))
    enriched = ranked[ranked["logFC"] > 0]["cell_type"].tolist()
    depleted = ranked[ranked["logFC"] < 0]["cell_type"].tolist()

    n_each = n_top // 2
    picked_enr = enriched[:n_each]
    picked_dep = depleted[:n_each]

    # Redistribute unused budget to the other direction
    slack = n_top - len(picked_enr) - len(picked_dep)
    if slack > 0:
        if len(picked_enr) < n_each:
            picked_dep = depleted[:n_each + slack]
        else:
            picked_enr = enriched[:n_each + slack]

    # Order: most depleted first, then most enriched
    all_ranked = ranked.set_index("cell_type")
    dep_sorted = sorted(picked_dep,
                        key=lambda c: all_ranked.loc[c, "logFC"] if c in all_ranked.index else 0)
    enr_sorted = sorted(picked_enr,
                        key=lambda c: all_ranked.loc[c, "logFC"] if c in all_ranked.index else 0,
                        reverse=True)
    top_cts = dep_sorted + enr_sorted

    ncols = 4
    nrows = int(np.ceil(len(top_cts) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.2, nrows * 3.0))
    axes = axes.flatten()

    for ax, ct in zip(axes, top_cts):
        row = ranked[ranked["cell_type"] == ct]
        fdr   = row["adj.P.Val"].values[0] if len(row) else np.nan
        lfc   = row["logFC"].values[0] if len(row) else np.nan

        data   = [props.loc[props[contrast] == g, ct].values * 100
                  for g in groups]
        colors = [PALETTE.get(g, "#95a5a6") for g in groups]

        bp = ax.boxplot(data, patch_artist=True, widths=0.45, showfliers=False,
                        medianprops=dict(color="black", linewidth=1.5))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        for j, (g, d) in enumerate(zip(groups, data)):
            jitter = np.random.default_rng(abs(hash(ct)) % (2**31)).uniform(
                -0.15, 0.15, size=len(d))
            ax.scatter(np.full(len(d), j + 1) + jitter, d,
                       color=PALETTE.get(g, "#95a5a6"),
                       alpha=0.75, s=18, zorder=3)

        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups, fontsize=9)
        ax.set_ylabel("Proportion (%)", fontsize=8)
        cls = _cell_class(ct)
        ax.set_title(
            f"{ct}\n[{cls[:3]}]  logFC={lfc:.2f}  FDR={fdr:.2g}",
            fontsize=7.5
        )

    for ax in axes[len(top_cts):]:
        ax.set_visible(False)

    fig.suptitle(f"Top {len(top_cts)} cell types by crumblr FDR", y=1.01)
    plt.tight_layout()
    fig.savefig(output_dir / "proportion_boxplots.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "proportion_boxplots.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("  Saved: proportion_boxplots.{pdf,png}")


def plot_scatter_continuous(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: str,
    results: pd.DataFrame,
    output_dir: Path,
    n_top: int = 12,
    color_by: str | None = None,
) -> None:
    """Scatter: proportion % vs continuous contrast, OLS regression line ± 1 SE.

    Points are colored by color_by (using PALETTE) when supplied.
    Cell types are selected balanced by sign of dream() slope, same logic as
    proportion_boxplots.  logFC / FDR in titles come from dream() (CLR scale,
    covariate-adjusted); the regression line is marginal OLS on proportion %.
    """
    from scipy import stats as sp_stats

    props  = counts.div(counts.sum(axis=1), axis=0)
    x_meta = metadata.loc[props.index, contrast].astype(float)

    # ── Point colors ──────────────────────────────────────────────────────────
    if color_by and color_by in metadata.columns:
        color_series   = metadata.loc[props.index, color_by].astype(str)
        pt_colors      = [PALETTE.get(v, "#95a5a6") for v in color_series]
        legend_handles = [
            mpatches.Patch(color=PALETTE.get(v, "#95a5a6"), label=v)
            for v in sorted(color_series.unique())
        ]
    else:
        pt_colors      = ["#95a5a6"] * len(props)
        legend_handles = []

    # ── Select top N balanced by direction ────────────────────────────────────
    ranked   = results.copy().sort_values("adj.P.Val").drop_duplicates("cell_type")
    enriched = ranked[ranked["logFC"] > 0]["cell_type"].tolist()
    depleted = ranked[ranked["logFC"] < 0]["cell_type"].tolist()
    n_each   = n_top // 2
    picked_enr = enriched[:n_each]
    picked_dep = depleted[:n_each]
    slack = n_top - len(picked_enr) - len(picked_dep)
    if slack > 0:
        if len(picked_enr) < n_each:
            picked_dep = depleted[:n_each + slack]
        else:
            picked_enr = enriched[:n_each + slack]
    all_ranked = ranked.set_index("cell_type")
    dep_sorted = sorted(picked_dep,
                        key=lambda c: all_ranked.loc[c, "logFC"] if c in all_ranked.index else 0)
    enr_sorted = sorted(picked_enr,
                        key=lambda c: all_ranked.loc[c, "logFC"] if c in all_ranked.index else 0,
                        reverse=True)
    top_cts = dep_sorted + enr_sorted

    # Reserve one extra panel for the legend if needed
    n_panels = len(top_cts) + (1 if legend_handles else 0)
    ncols    = 4
    nrows    = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.2))
    axes = np.array(axes).flatten()

    for ax, ct in zip(axes, top_cts):
        row = all_ranked.loc[ct] if ct in all_ranked.index else None
        fdr = row["adj.P.Val"] if row is not None else np.nan
        lfc = row["logFC"]     if row is not None else np.nan

        y_pct  = props[ct] * 100
        x_arr  = x_meta.values.astype(float)
        y_arr  = y_pct.values.astype(float)
        c_arr  = np.array(pt_colors)

        mask   = np.isfinite(x_arr) & np.isfinite(y_arr)
        xm, ym, cm = x_arr[mask], y_arr[mask], c_arr[mask]

        ax.scatter(xm, ym, c=cm, alpha=0.70, s=18, zorder=3)

        if len(xm) > 2:
            sl, ic, _, _, _ = sp_stats.linregress(xm, ym)
            xl   = np.linspace(xm.min(), xm.max(), 120)
            yl   = sl * xl + ic
            # SE of regression line (not prediction interval)
            resid   = ym - (sl * xm + ic)
            s_res   = np.sqrt(np.sum(resid ** 2) / (len(xm) - 2))
            ss_x    = np.sum((xm - xm.mean()) ** 2)
            se_line = s_res * np.sqrt(1 / len(xm) + (xl - xm.mean()) ** 2 / ss_x)
            ax.plot(xl, yl, color="black", linewidth=1.5, zorder=4)
            ax.fill_between(xl, yl - se_line, yl + se_line,
                            color="black", alpha=0.12, zorder=2)

        cls = _cell_class(ct)
        ax.set_title(f"{ct}\n[{cls[:3]}]  logFC={lfc:.3f}  FDR={fdr:.2g}",
                     fontsize=7.5)
        ax.set_xlabel(contrast, fontsize=8)
        ax.set_ylabel("Proportion (%)", fontsize=8)

    # ── Legend panel ──────────────────────────────────────────────────────────
    if legend_handles:
        ax_leg = axes[len(top_cts)]
        ax_leg.axis("off")
        ax_leg.legend(handles=legend_handles, title=color_by or "",
                      loc="center", frameon=False, fontsize=9)
        for ax in axes[len(top_cts) + 1:]:
            ax.set_visible(False)
    else:
        for ax in axes[len(top_cts):]:
            ax.set_visible(False)

    fig.suptitle(
        f"Top {len(top_cts)} cell types by crumblr FDR  ({contrast})\n"
        f"logFC/FDR: dream() covariate-adjusted [CLR scale] — "
        f"line: marginal OLS ± 1 SE",
        fontsize=9, y=1.02,
    )
    plt.tight_layout()
    fig.savefig(output_dir / "proportion_scatter.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "proportion_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: proportion_scatter.{pdf,png}")


def plot_volcano(results: pd.DataFrame, contrast: str, output_dir: Path) -> None:
    """Volcano: log2FC vs –log10(BH-adjusted p), coloured by FDR significance."""
    for coeff, sub in results.groupby("coefficient"):
        sub = sub.copy()
        sub["neg_log10_fdr"] = -np.log10(sub["adj.P.Val"].clip(lower=1e-300))
        sub["significant"]   = sub["adj.P.Val"] < 0.05
        sub["cell_class"]    = sub["cell_type"].map(_cell_class)

        fig, ax = plt.subplots(figsize=(7, 5))
        for cls, grp in sub.groupby("cell_class"):
            marker = {"Dopaminergic": "o", "Excitatory": "s",
                      "Inhibitory": "^"}.get(cls, "o")
            ax.scatter(
                grp["logFC"], grp["neg_log10_fdr"],
                c=[CLASS_COLORS[cls]] * len(grp),
                alpha=0.8, s=55, marker=marker,
                edgecolors=["black" if s else "none" for s in grp["significant"]],
                linewidths=[0.8 if s else 0.0 for s in grp["significant"]],
                label=cls,
            )

        # FDR = 0.05 threshold line
        fdr_line = -np.log10(0.05)
        ax.axhline(fdr_line, color="grey", linestyle="--", linewidth=0.8)
        ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)

        # Label significant hits
        sig = sub[sub["significant"]]
        for _, r in sig.iterrows():
            ax.text(r["logFC"], r["neg_log10_fdr"] + 0.12,
                    r["cell_type"], fontsize=7.5, ha="center", va="bottom")

        ax.set_xlabel("log2 fold-change (CLR scale)", fontsize=11)
        ax.set_ylabel("–log10(BH-adjusted p)", fontsize=11)
        ax.set_title(f"Volcano – {coeff}", fontsize=12)
        ax.legend(title="Cell class", frameon=False, fontsize=9)
        plt.tight_layout()

        safe = _safe_name(coeff)
        fig.savefig(output_dir / f"volcano_{safe}.pdf", bbox_inches="tight")
        fig.savefig(output_dir / f"volcano_{safe}.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: volcano_{safe}.{{pdf,png}}")


def plot_clr_heatmap(
    clr_df: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: str,
    output_dir: Path,
    color_by: str | None = None,
) -> None:
    """Hierarchically-clustered heatmap of CLR values.

    Samples are sorted by contrast (works for both numeric and categorical).
    Row colour annotation uses color_by when supplied (e.g. case_control when
    contrast is a continuous variable like age); defaults to contrast itself.
    """
    ann_col = color_by if (color_by and color_by in metadata.columns) else contrast

    shared = clr_df.index.intersection(metadata.index)
    clr    = clr_df.loc[shared]
    meta   = metadata.loc[shared]

    # Sort samples by contrast value (numeric → ascending age; categorical → alpha)
    try:
        order = meta[contrast].argsort()
    except TypeError:
        order = meta[contrast].astype(str).argsort()
    clr  = clr.iloc[order]
    meta = meta.iloc[order]

    # Cluster cell types
    ct_dist  = pdist(clr.T, metric="correlation")
    ct_link  = hierarchy.linkage(ct_dist, method="average")
    ct_order = hierarchy.leaves_list(ct_link)
    clr      = clr.iloc[:, ct_order]

    row_colors = meta[ann_col].astype(str).map(PALETTE).fillna("#95a5a6")

    g = sns.clustermap(
        clr,
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        cmap="RdBu_r",
        center=0,
        figsize=(max(14, len(clr.columns) * 0.38),
                 max(8,  len(clr)          * 0.18)),
        yticklabels=False,
        xticklabels=True,
        cbar_kws={"label": "CLR-transformed proportion"},
    )
    g.ax_heatmap.set_xlabel("Cell type")
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7
    )
    handles = [mpatches.Patch(color=v, label=k) for k, v in PALETTE.items()
               if k in meta[ann_col].astype(str).unique()]
    g.ax_heatmap.legend(
        handles=handles, title=ann_col,
        bbox_to_anchor=(1.25, 1), loc="upper left", frameon=False
    )
    sort_label = f"sorted by {contrast}" + (
        f", annotated by {ann_col}" if ann_col != contrast else ""
    )
    g.figure.suptitle(f"CLR proportion heatmap (samples × cell types)\n{sort_label}",
                      y=1.01)

    g.figure.savefig(output_dir / "clr_heatmap.pdf", bbox_inches="tight")
    g.figure.savefig(output_dir / "clr_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close("all")
    print("  Saved: clr_heatmap.{pdf,png}")


def plot_da_class_summary(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    contrast: str,
    output_dir: Path,
) -> None:
    """Cell-class-level proportion summary across contrast groups."""
    props = counts.div(counts.sum(axis=1), axis=0)
    props.columns = [_cell_class(c) for c in props.columns]
    # Sum within class
    class_props = props.T.groupby(level=0).sum().T
    class_props[contrast] = metadata.loc[class_props.index, contrast].astype(str)

    class_long = class_props.melt(
        id_vars=contrast, var_name="cell_class", value_name="proportion"
    )

    groups = sorted(class_long[contrast].unique())
    classes = ["Dopaminergic", "Excitatory", "Inhibitory"]

    fig, axes = plt.subplots(1, 3, figsize=(11, 4), sharey=False)
    for ax, cls in zip(axes, classes):
        sub = class_long[class_long["cell_class"] == cls]
        data   = [sub.loc[sub[contrast] == g, "proportion"].values * 100
                  for g in groups]
        colors = [PALETTE.get(g, "#95a5a6") for g in groups]

        bp = ax.boxplot(data, patch_artist=True, widths=0.5, showfliers=False,
                        medianprops=dict(color="black", linewidth=2))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)

        for j, (g, d) in enumerate(zip(groups, data)):
            jitter = np.random.default_rng(abs(hash(cls + g)) % (2**31)).uniform(
                -0.15, 0.15, size=len(d))
            ax.scatter(np.full(len(d), j + 1) + jitter, d,
                       color=PALETTE.get(g, "#95a5a6"),
                       alpha=0.8, s=30, zorder=3)

        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups)
        ax.set_ylabel("Proportion (%)")
        ax.set_title(f"{cls} Neurons", color=CLASS_COLORS.get(cls, "black"),
                     fontsize=11, fontweight="bold")

    fig.suptitle(f"Broad cell-class proportions by {contrast}", fontsize=12)
    plt.tight_layout()
    fig.savefig(output_dir / "cell_class_summary.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "cell_class_summary.png", dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("  Saved: cell_class_summary.{pdf,png}")


def plot_logfc_bar(results: pd.DataFrame, output_dir: Path) -> None:
    """Horizontal bar chart of log2FC for all cell types, coloured by class."""
    for coeff, sub in results.groupby("coefficient"):
        sub = sub.copy().sort_values("logFC")
        sub["cell_class"] = sub["cell_type"].map(_cell_class)
        sub["color"]      = sub["cell_class"].map(CLASS_COLORS)
        sub["edge"]       = np.where(sub["adj.P.Val"] < 0.05, "black", "none")
        sub["lw"]         = np.where(sub["adj.P.Val"] < 0.05, 1.2, 0)
        # SE not returned by topTable; derive from logFC / t
        sub["SE"] = np.where(
            sub["t"].abs() > 1e-10,
            (1.96 * sub["logFC"] / sub["t"]).abs(),
            np.nan,
        )

        fig, ax = plt.subplots(figsize=(7, max(5, len(sub) * 0.30)))
        ax.barh(
            range(len(sub)), sub["logFC"],
            color=sub["color"], alpha=0.80,
            edgecolor=sub["edge"], linewidth=sub["lw"]
        )
        ax.errorbar(
            sub["logFC"], range(len(sub)),
            xerr=sub["SE"], fmt="none",
            color="black", linewidth=0.8, capsize=2
        )
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["cell_type"], fontsize=8)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("log2 fold-change (CLR scale)")
        ax.set_title(f"crumblr log2FC – {coeff}\n(black border = FDR < 0.05)")

        # Class legend
        handles = [mpatches.Patch(color=c, label=cls)
                   for cls, c in CLASS_COLORS.items()]
        ax.legend(handles=handles, title="Cell class",
                  bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

        plt.tight_layout()
        safe = _safe_name(coeff)
        fig.savefig(output_dir / f"logfc_bar_{safe}.pdf", bbox_inches="tight")
        fig.savefig(output_dir / f"logfc_bar_{safe}.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: logfc_bar_{safe}.{{pdf,png}}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--csv",            required=True,       help="Cell-level CSV path")
    p.add_argument("--proportion-col", required=True,       help="Column to compute proportions of (e.g. cell_type)")
    p.add_argument("--sample-col",     default=None,        help="Column identifying each sample. Defaults to --donor-col when omitted (for datasets where donor IS the sample unit).")
    p.add_argument("--donor-col",      default=None,        help="Donor/subject column (for random effects bookkeeping)")
    p.add_argument("--contrast",       required=True,       help="Main contrast variable (e.g. case_control)")
    p.add_argument("--fixed-effects",  nargs="+", default=[], help="Fixed-effect covariate column names")
    p.add_argument("--random-effects", nargs="*", default=[], help="Random-effect column names (will be wrapped in (1|...))")
    p.add_argument("--ref-level",      default=None,        help="Reference level for contrast variable")
    p.add_argument("--output-dir",     required=True,       help="Output directory")
    p.add_argument("--r-libs",   default=None,
                   help="Path to custom R library (sets R_LIBS_CUSTOM). "
                        "Defaults to <sandbox>/R_libs if that directory exists.")
    p.add_argument("--r-script", default=None,
                   help="Path to crumblr_analysis.R. "
                        "Defaults to <dir of this script>/crumblr_analysis.R.")
    p.add_argument("--n-top",          type=int, default=60,  help="Total cell types shown in boxplots/scatter; budget split evenly between enriched and depleted (default: 60)")
    p.add_argument("--min-n-cells",    type=int, default=100, help="Minimum total cells per sample; samples below this threshold are dropped before analysis (default: 100)")
    p.add_argument("--min-n-cells-category",    type=int, default=150, help="Minimum total cells per category per category; categories below this threshold are dropped before analysis (default: 150)")
    p.add_argument("--color-by",       default=None,         help="Column used to colour points/heatmap rows when --contrast is continuous (e.g. case_control). Ignored for categorical contrasts.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    args_dict = {k: (str(v) if isinstance(v, Path) else v)
                 for k, v in vars(args).items()}
    with open(out_dir / "run_args.json", "w") as _f:
        json.dump(args_dict, _f, indent=2)

    r_script = Path(args.r_script) if args.r_script else \
               Path(__file__).parent / "crumblr_analysis.R"
    if not r_script.exists():
        raise FileNotFoundError(f"R script not found: {r_script}")

    if args.r_libs:
        r_libs = Path(args.r_libs)
    else:
        default_r_libs = Path(__file__).parent.parent / "R_libs"
        r_libs = default_r_libs if default_r_libs.is_dir() else None

    # ── Resolve sample column ─────────────────────────────────────────────────
    sample_col = args.sample_col or args.donor_col
    if not sample_col:
        raise ValueError(
            "Specify --donor-col (or --sample-col) to identify the sample unit."
        )

    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"[1/6] Loading {args.csv} …")
    df = pd.read_csv(args.csv, low_memory=False)

    required_cols = {args.proportion_col, sample_col, args.contrast}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    if args.min_n_cells_category > 0:
        category_counts = df[args.proportion_col].value_counts()
        df = df[df[args.proportion_col].isin(category_counts[category_counts >= args.min_n_cells_category].index)]

    # ── Aggregate to counts ───────────────────────────────────────────────────
    print("[2/6] Aggregating cell counts …")
    counts = aggregate_counts(df, args.proportion_col, sample_col)

    n_cells = counts.sum(axis=1)
    low = n_cells[n_cells < args.min_n_cells]
    if len(low):
        print(f"  Dropping {len(low)} sample(s) with < {args.min_n_cells} cells: "
              f"{low.index.tolist()}")
        counts = counts.loc[n_cells >= args.min_n_cells]
    if len(counts) == 0:
        raise ValueError(
            f"No samples remain after applying --min-n-cells {args.min_n_cells}."
        )

    all_meta_cols = list({args.contrast} |
                         set(args.fixed_effects) |
                         set(args.random_effects) |
                         ({args.donor_col} if args.donor_col else set()))
    all_meta_cols = [c for c in all_meta_cols if c in df.columns]

    metadata = build_metadata(df, sample_col, all_meta_cols)
    metadata = relevel(metadata, args.contrast, args.ref_level)
    validate_random_effects(metadata, sample_col, args.random_effects)

    # ── Auto-detect contrast type ──────────────────────────────────────────────
    is_continuous = pd.api.types.is_numeric_dtype(metadata[args.contrast])
    print(f"  Contrast '{args.contrast}' detected as "
          f"{'continuous' if is_continuous else 'categorical'}.")

    print(f"  Samples: {len(counts)}  |  Cell types: {counts.shape[1]}")
    counts_path  = out_dir / "counts_matrix.csv"
    meta_path    = out_dir / "metadata.csv"
    counts.to_csv(counts_path)
    metadata.to_csv(meta_path)

    # ── Build formula ─────────────────────────────────────────────────────────
    formula = build_formula(args.contrast, args.fixed_effects, args.random_effects)
    print(f"[3/6] Formula: {formula}")

    # ── Run R crumblr ─────────────────────────────────────────────────────────
    print("[4/6] Running crumblr (R) …")
    run_r_analysis(counts, metadata, formula, args.contrast,
                   out_dir, r_script, r_libs)

    # ── Load results ──────────────────────────────────────────────────────────
    results_path = out_dir / "crumblr_results.csv"
    if not results_path.exists():
        raise FileNotFoundError(
            "crumblr_results.csv not produced; check R output above."
        )
    results = pd.read_csv(results_path)
    clr_path = out_dir / "clr_values.csv"
    clr_df   = pd.read_csv(clr_path, index_col=0)

    # ── Generate plots ────────────────────────────────────────────────────────
    print("[5/6] Generating plots …")

    if is_continuous:
        plot_scatter_continuous(counts, metadata, args.contrast, results, out_dir,
                                n_top=args.n_top, color_by=args.color_by)
        plot_clr_heatmap(clr_df, metadata, args.contrast, out_dir,
                         color_by=args.color_by)
    else:
        plot_composition_bar(counts, metadata, args.contrast, out_dir)
        plot_da_class_summary(counts, metadata, args.contrast, out_dir)
        plot_proportion_boxplots(counts, metadata, args.contrast, results, out_dir,
                                 n_top=args.n_top)
        plot_clr_heatmap(clr_df, metadata, args.contrast, out_dir)
    plot_volcano(results, args.contrast, out_dir)
    plot_logfc_bar(results, out_dir)

    # ── Summary table ─────────────────────────────────────────────────────────
    print("[6/6] Summary")
    sig = results[results["adj.P.Val"] < 0.05].sort_values("adj.P.Val")
    print(f"\n  Significant cell types (FDR < 0.05): {len(sig)}")
    if len(sig):
        cols = ["cell_type", "logFC", "adj.P.Val", "coefficient"]
        print(sig[cols].to_string(index=False))

    print(f"\nAll outputs written to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
