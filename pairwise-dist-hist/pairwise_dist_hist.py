#!/usr/bin/env python3
"""
pairwise_dist_hist.py  –  pairwise embedding-distance density plots

For every pair of cell-type labels (including same-type pairs), samples
balanced pairs of cells and plots KDE density histograms of their distances
in an obsm embedding space.  Three overlaid KDE curves per panel show:
  • within-type distances for ct_i  (when applicable)
  • cross-type distances  ct_i ↔ ct_j
  • within-type distances for ct_j  (when applicable)

Environment Setup (one-time)
-----------------------------
    bash /path/to/pairwise-dist-hist/setup_env.sh

Running the Script
------------------
No manual environment activation is needed.  The script automatically
re-execs itself with the adata-venv Python if it is not already running
inside it:

    python /path/to/pairwise-dist-hist/pairwise_dist_hist.py \\
        --h5ad   cells.h5ad   \\
        --obs-col cell_type   \\
        --obsm-key X_pca      \\
        --output  distances.pdf

Usage
-----
python pairwise_dist_hist.py \\
    --h5ad      <path>          (required) \\
    --obs-col   <col>           (required, obs column with cell-type labels) \\
    --obsm-key  <key>           (required, obsm key, e.g. X_pca) \\
    --output    <path.pdf>      (required) \\
    [--n-pairs  10000]          pairs sampled per sub-pairing category \\
    [--metric   cosine]         euclidean | cosine | manhattan \\
    [--seed     42]

Outputs
-------
  <output>.pdf/.png/.svg              grid of all pairs (KDE per panel)
  <output>_joint.pdf/.png/.svg        cross-type pairs, combined + components
  <output>_matrix_joint.pdf/.png/.svg upper-triangle matrix, pooled joint KDE
  <output>_matrix_between.pdf/.png/.svg upper-triangle matrix, cross-only KDE
  <output>_matrix_data.csv            tidy distance data for both matrix figs
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# ── Self-activating: re-exec with adata-venv Python if not already using it ──
_VENV_DIR    = Path(__file__).parent.parent / "adata-venv"
_VENV_PYTHON = _VENV_DIR / "bin" / "python3"
# Compare sys.prefix (set by pyvenv.cfg) rather than the resolved binary path.
# Resolving the symlink would make the venv python appear identical to the
# base conda/system python and prevent the re-exec from firing.
if Path(sys.prefix).resolve() != _VENV_DIR.resolve():
    if not _VENV_PYTHON.is_file():
        sys.exit(
            f"ERROR: adata-venv not found at {_VENV_DIR}.\n"
            f"Run setup first:  bash {Path(__file__).parent / 'setup_env.sh'}"
        )
    os.execv(str(_VENV_PYTHON), [str(_VENV_PYTHON)] + sys.argv)
# ─────────────────────────────────────────────────────────────────────────────

import argparse
import itertools
import math
import warnings

import anndata
import diptest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as scipy_stats
from scipy.spatial.distance import cdist

# ── Distance helpers ──────────────────────────────────────────────────────────

_SUPPORTED_METRICS = ("euclidean", "cosine", "manhattan")


def _compute_distances(
    a: np.ndarray,
    b: np.ndarray,
    metric: str,
) -> np.ndarray:
    """Return per-row pairwise distances between matched rows of a and b."""
    if metric == "euclidean":
        diff = a - b
        return np.sqrt((diff * diff).sum(axis=1))
    if metric == "manhattan":
        return np.abs(a - b).sum(axis=1)
    if metric == "cosine":
        # 1 - cosine_similarity
        norm_a = np.linalg.norm(a, axis=1, keepdims=True)
        norm_b = np.linalg.norm(b, axis=1, keepdims=True)
        # Avoid division by zero
        denom = norm_a * norm_b
        denom = np.where(denom == 0, 1.0, denom)
        sim = (a * b).sum(axis=1) / denom.squeeze()
        return 1.0 - sim
    raise ValueError(
        f"Unsupported metric '{metric}'. Choose from: {_SUPPORTED_METRICS}"
    )


# ── Balanced sampling ─────────────────────────────────────────────────────────

def _balanced_sample(
    idx_a: np.ndarray,
    idx_b: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Return two arrays of equal length by subsampling the larger group.

    If group A has n_a cells and group B has n_b cells with n_a > n_b, then
    group A is subsampled to n_b cells (without replacement).  This makes the
    sampling rate of the majority class 1/n_a and of the minority class
    1/n_b, so the minority class is effectively sampled at rate n_a/n_b times
    higher — exactly balancing the representation as requested.
    """
    n_a, n_b = len(idx_a), len(idx_b)
    if n_a == n_b:
        return idx_a, idx_b
    if n_a > n_b:
        return rng.choice(idx_a, size=n_b, replace=False), idx_b
    return idx_a, rng.choice(idx_b, size=n_a, replace=False)


def _sample_pairs(
    pool: np.ndarray,
    n_pairs: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample n_pairs random pairs (with replacement) from pool."""
    left  = rng.integers(0, len(pool), size=n_pairs)
    right = rng.integers(0, len(pool), size=n_pairs)
    return pool[left], pool[right]


def _cross_pairs(
    pool_a: np.ndarray,
    pool_b: np.ndarray,
    n_pairs: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample n_pairs cross-type pairs (with replacement) from two pools."""
    left  = pool_a[rng.integers(0, len(pool_a), size=n_pairs)]
    right = pool_b[rng.integers(0, len(pool_b), size=n_pairs)]
    return left, right


# ── Core computation ──────────────────────────────────────────────────────────

def compute_pair_distances(
    embedding: np.ndarray,
    obs_labels: np.ndarray,
    ct_i: str,
    ct_j: str,
    n_pairs: int,
    metric: str,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """
    For a pair of cell types (ct_i, ct_j), return distances for three
    sub-pairings:
      "within_i"  – within ct_i
      "cross"     – ct_i ↔ ct_j
      "within_j"  – within ct_j  (omitted when ct_i == ct_j)

    Cells are balanced so that ct_i and ct_j contribute equal numbers to
    each pool before pairs are drawn.
    """
    idx_i = np.where(obs_labels == ct_i)[0]
    idx_j = np.where(obs_labels == ct_j)[0] if ct_j != ct_i else idx_i

    if len(idx_i) == 0:
        raise ValueError(f"No cells found for label '{ct_i}'.")
    if len(idx_j) == 0:
        raise ValueError(f"No cells found for label '{ct_j}'.")

    same_type = ct_i == ct_j

    if same_type:
        # Only one pool; within-type only
        left, right = _sample_pairs(idx_i, n_pairs, rng)
        return {
            "within_i": _compute_distances(
                embedding[left], embedding[right], metric
            )
        }

    # Balanced pools
    pool_i, pool_j = _balanced_sample(idx_i, idx_j, rng)

    left_i, right_i = _sample_pairs(pool_i, n_pairs, rng)
    left_j, right_j = _sample_pairs(pool_j, n_pairs, rng)
    left_c, right_c = _cross_pairs(pool_i, pool_j, n_pairs, rng)

    return {
        "within_i": _compute_distances(
            embedding[left_i], embedding[right_i], metric
        ),
        "cross": _compute_distances(
            embedding[left_c], embedding[right_c], metric
        ),
        "within_j": _compute_distances(
            embedding[left_j], embedding[right_j], metric
        ),
    }


def _build_dist_cache(
    embedding: np.ndarray,
    obs_labels: np.ndarray,
    ct_list: list[str],
    n_pairs: int,
    metric: str,
    rng: np.random.Generator,
) -> dict[tuple[str, str], dict]:
    """Compute and return all within/cross distances for every ct pair.

    Returns a dict keyed by (ct_i, ct_j) tuples.  Diagonal keys (ct, ct)
    contain "within_i".  Off-diagonal keys (ct_i, ct_j) with i < j contain
    "within_i", "cross", "within_j", and "extra_cross" (a second independent
    cross sample used to build the 1:2:1 pooled joint distribution).
    """
    cache: dict[tuple[str, str], dict] = {}

    for ct in ct_list:
        cache[(ct, ct)] = compute_pair_distances(
            embedding, obs_labels, ct, ct, n_pairs, metric, rng
        )

    for ct_i, ct_j in itertools.combinations(ct_list, 2):
        d = compute_pair_distances(
            embedding, obs_labels, ct_i, ct_j, n_pairs, metric, rng
        )
        idx_i = np.where(obs_labels == ct_i)[0]
        idx_j = np.where(obs_labels == ct_j)[0]
        pool_i, pool_j = _balanced_sample(idx_i, idx_j, rng)
        lc2, rc2 = _cross_pairs(pool_i, pool_j, n_pairs, rng)
        d["extra_cross"] = _compute_distances(embedding[lc2], embedding[rc2], metric)
        cache[(ct_i, ct_j)] = d

    return cache


# ── Joint-distribution helpers ────────────────────────────────────────────────

def _bimodality_coefficient(x: np.ndarray) -> float:
    """Sarle's bimodality coefficient.

    BC > 5/9 ≈ 0.556 is the conventional threshold suggesting bimodality.
    """
    n  = len(x)
    sk = float(scipy_stats.skew(x))
    ku = float(scipy_stats.kurtosis(x))          # Fisher's excess kurtosis (normal = 0)
    correction = 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return (sk ** 2 + 1.0) / (ku + correction)


def _fmt_pval(p: float) -> str:
    if p < 1e-15:
        return "< 1×10⁻¹⁵"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def _dip_stars(p: float) -> str:
    """Return significance stars for Hartigan dip test p-value."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def plot_joint(
    embedding: np.ndarray,
    obs_labels: np.ndarray,
    ct_list: list[str],
    n_pairs: int,
    metric: str,
    rng: np.random.Generator,
    ncols: int = 1,
) -> plt.Figure:
    """Grid of combined KDE panels, one per cross-type pair (same-type omitted).

    Each panel pools within_i + cross + within_j into a single distribution,
    overlays the three component KDEs for reference, and annotates with
    Hartigan's dip test and Sarle's bimodality coefficient.
    Legend is placed outside each panel to the upper right.
    """
    cross_pairs = list(itertools.combinations(ct_list, 2))
    n_panels = len(cross_pairs)
    nrows    = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 5.5, nrows * 3.5),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    metric_label = {
        "euclidean": "Euclidean distance",
        "cosine":    "Cosine distance",
        "manhattan": "Manhattan distance",
    }.get(metric, metric)

    # diptest table cap
    _DIP_MAX = 72_000

    for panel_idx, (ct_i, ct_j) in enumerate(cross_pairs):
        ax = axes_flat[panel_idx]

        try:
            dists = compute_pair_distances(
                embedding, obs_labels, ct_i, ct_j,
                n_pairs=n_pairs, metric=metric, rng=rng,
            )
        except ValueError as exc:
            ax.set_visible(False)
            warnings.warn(str(exc))
            continue

        # Build combined with correct proportions: when drawing random pairs from
        # the merged pool, P(i–j) = 2 × P(i–i) = 2 × P(j–j), so cross-type
        # distances should appear twice as often as each within-type category.
        # Sample an extra n_pairs of cross distances (reusing the balanced pools
        # already constructed inside compute_pair_distances) so that
        # combined = within_i (n) + cross (2n) + within_j (n) = 4n total.
        idx_i  = np.where(obs_labels == ct_i)[0]
        idx_j  = np.where(obs_labels == ct_j)[0]
        pool_i, pool_j = _balanced_sample(idx_i, idx_j, rng)
        lc2, rc2 = _cross_pairs(pool_i, pool_j, n_pairs, rng)
        extra_cross = _compute_distances(embedding[lc2], embedding[rc2], metric)

        combined = np.concatenate(
            [dists["within_i"], dists["cross"], extra_cross, dists["within_j"]]
        )

        # ── Bimodality tests ──────────────────────────────────────────────────
        dip_input = combined if len(combined) <= _DIP_MAX else rng.choice(combined, size=_DIP_MAX, replace=False)
        dip_stat, dip_pval = diptest.diptest(dip_input)
        bc = _bimodality_coefficient(combined)
        bc_label = "bimodal" if bc > 5.0 / 9.0 else "unimodal"

        # ── KDEs ──────────────────────────────────────────────────────────────
        color_i     = _ct_color(ct_i, ct_list)
        color_j     = _ct_color(ct_j, ct_list)
        color_cross = tuple((a + b) / 2 for a, b in zip(color_i, color_j))

        # Combined distribution: filled gray background
        sns.kdeplot(
            combined, ax=ax,
            color="0.4", linewidth=2.0,
            fill=True, alpha=0.12,
            label=f"combined (n={len(combined):,})",
            warn_singular=False,
        )
        # Component curves overlaid (same style as plot_grid)
        sns.kdeplot(
            dists["within_i"], ax=ax,
            color=color_i, linewidth=1.8,
            fill=True, alpha=0.18,
            label=f"{ct_i}–{ct_i}", warn_singular=False,
        )
        sns.kdeplot(
            dists["cross"], ax=ax,
            color=color_cross, linewidth=2.2,
            fill=True, alpha=0.18, linestyle="--",
            label=f"{ct_i}–{ct_j}", warn_singular=False,
        )
        sns.kdeplot(
            dists["within_j"], ax=ax,
            color=color_j, linewidth=1.8,
            fill=True, alpha=0.18,
            label=f"{ct_j}–{ct_j}", warn_singular=False,
        )

        ax.set_xlabel(metric_label, fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(f"{ct_i}  ×  {ct_j}", fontsize=9, fontweight="bold")

        # Legend outside the axes, anchored to the upper-right corner
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0,
            fontsize=7,
            frameon=True,
            framealpha=0.85,
            edgecolor="0.75",
        )

        # ── Per-panel bimodality annotation ───────────────────────────────────
        annot = (
            f"Dip D={dip_stat:.3f}  p={_fmt_pval(dip_pval)}\n"
            f"BC={bc:.3f}  ({bc_label})"
        )
        ax.text(
            0.97, 0.97, annot,
            transform=ax.transAxes,
            fontsize=7, verticalalignment="top", horizontalalignment="right",
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="0.75"),
        )

    for ax in axes_flat[n_panels:]:
        ax.set_visible(False)

    # suptitle is patched in main() with the actual obsm key
    fig.suptitle(
        f"Joint pairwise {metric_label} distributions  "
        f"(n_pairs={n_pairs:,})",
        fontsize=11, y=1.01,
    )
    # plt.tight_layout()
    fig.subplots_adjust(hspace=0.55)
    return fig


# ── Plotting ──────────────────────────────────────────────────────────────────

# Seaborn tab10-like palette for cell types
_CT_PALETTE = sns.color_palette("tab10", n_colors=20)


def _ct_color(ct: str, ct_list: list[str]) -> tuple:
    idx = ct_list.index(ct) % len(_CT_PALETTE)
    return _CT_PALETTE[idx]


def plot_grid(
    embedding: np.ndarray,
    obs_labels: np.ndarray,
    ct_list: list[str],
    n_pairs: int,
    metric: str,
    rng: np.random.Generator,
    ncols: int = 3,
) -> plt.Figure:
    """Build the full grid figure.  Returns the Figure."""

    pairs = list(itertools.combinations_with_replacement(ct_list, 2))
    n_panels = len(pairs)
    nrows    = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * 4.5, nrows * 3.5),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    metric_label = {"euclidean": "Euclidean distance",
                    "cosine":    "Cosine distance",
                    "manhattan": "Manhattan distance"}.get(metric, metric)

    for panel_idx, (ct_i, ct_j) in enumerate(pairs):
        ax = axes_flat[panel_idx]

        try:
            dists = compute_pair_distances(
                embedding, obs_labels, ct_i, ct_j,
                n_pairs=n_pairs, metric=metric, rng=rng,
            )
        except ValueError as exc:
            ax.set_visible(False)
            warnings.warn(str(exc))
            continue

        color_i = _ct_color(ct_i, ct_list)
        color_j = _ct_color(ct_j, ct_list)
        # Cross-type: average hue by blending the two colours
        color_cross = tuple((a + b) / 2 for a, b in zip(color_i, color_j))

        same_type = ct_i == ct_j

        if same_type:
            sns.kdeplot(
                dists["within_i"], ax=ax,
                color=color_i, linewidth=2.0,
                label=f"{ct_i}–{ct_i}",
                fill=True, alpha=0.25, warn_singular=False,
            )
        else:
            sns.kdeplot(
                dists["within_i"], ax=ax,
                color=color_i, linewidth=1.8,
                label=f"{ct_i}–{ct_i}",
                fill=True, alpha=0.20, warn_singular=False,
            )
            sns.kdeplot(
                dists["cross"], ax=ax,
                color=color_cross, linewidth=2.2,
                label=f"{ct_i}–{ct_j}",
                fill=True, alpha=0.20, linestyle="--", warn_singular=False,
            )
            sns.kdeplot(
                dists["within_j"], ax=ax,
                color=color_j, linewidth=1.8,
                label=f"{ct_j}–{ct_j}",
                fill=True, alpha=0.20, warn_singular=False,
            )

        ax.set_xlabel(metric_label, fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_title(f"{ct_i}  ×  {ct_j}", fontsize=9, fontweight="bold")
        ax.legend(fontsize=6.5, frameon=False)

    # Hide unused panels
    for ax in axes_flat[n_panels:]:
        ax.set_visible(False)

    fig.suptitle(
        f"Pairwise {metric_label} distributions  "
        f"(obsm key: {{}}, n_pairs={n_pairs:,})",
        fontsize=11, y=1.01,
    )
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.55)
    return fig


def plot_matrix_joint(
    dist_cache: dict[tuple[str, str], dict],
    ct_list: list[str],
    metric: str,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Upper-triangle N×N matrix — diagonal: within-type, upper: pooled joint distribution.

    The pooled joint combines within_i + 2×cross + within_j into a single KDE
    colored with a blend of the two cell-type colors.  Each upper-triangle panel
    is annotated with Hartigan's dip statistic (with significance stars) and
    Sarle's bimodality coefficient.  Lower triangle panels are hidden.

    Returns (figure, data_df) where data_df has columns:
        ct_i, ct_j, pair_type, distance
    """
    N = len(ct_list)
    _DIP_MAX = 72_000
    _rng = np.random.default_rng(0)

    metric_label = {
        "euclidean": "Euclidean distance",
        "cosine":    "Cosine distance",
        "manhattan": "Manhattan distance",
    }.get(metric, metric)

    panel_size = max(2.0, min(3.5, 18.0 / N))
    annot_fs   = max(4.5, 6.5 - 0.3 * N)
    fig, axes  = plt.subplots(N, N, figsize=(N * panel_size, N * panel_size), squeeze=False)

    for row in range(N):
        for col in range(N):
            ax      = axes[row, col]
            ct_row  = ct_list[row]
            ct_col  = ct_list[col]
            color_r = _ct_color(ct_row, ct_list)
            color_c = _ct_color(ct_col, ct_list)
            color_blend = tuple((a + b) / 2 for a, b in zip(color_r, color_c))

            if row == col:
                # ── Diagonal: within-type KDE on gray background ───────────
                ax.set_facecolor("#f5f5f5")
                d = dist_cache[(ct_row, ct_row)]
                sns.kdeplot(
                    d["within_i"], ax=ax,
                    color=color_r, linewidth=2.0,
                    fill=True, alpha=0.35, warn_singular=False,
                )

            elif row < col:
                # ── Upper triangle: single filled KDE of pooled joint ──────
                ct_i = ct_list[row]
                ct_j = ct_list[col]
                d = dist_cache[(ct_i, ct_j)]
                combined = np.concatenate(
                    [d["within_i"], d["cross"], d["extra_cross"], d["within_j"]]
                )

                dip_input = (combined if len(combined) <= _DIP_MAX
                             else _rng.choice(combined, size=_DIP_MAX, replace=False))
                dip_stat, dip_pval = diptest.diptest(dip_input)
                stars = _dip_stars(dip_pval)

                sns.kdeplot(
                    combined, ax=ax,
                    color=color_blend, linewidth=1.8,
                    fill=True, alpha=0.30, warn_singular=False,
                )

                annot = (
                    f"D={dip_stat:.3f}{stars}\n"
                    f"p={_fmt_pval(dip_pval)}"
                )
                ax.text(
                    0.97, 0.97, annot,
                    transform=ax.transAxes,
                    fontsize=annot_fs,
                    verticalalignment="top", horizontalalignment="right",
                    family="monospace",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              alpha=0.8, edgecolor="0.75"),
                )

            else:
                # ── Lower triangle: hidden ─────────────────────────────────
                ax.set_visible(False)
                continue

            ax.tick_params(labelsize=max(5, 7 - N // 4))
            ax.set_xlabel(metric_label, fontsize=max(5, 7 - N // 4))
            ax.set_ylabel("Density", fontsize=max(5, 7 - N // 4))

    # ── Shared X scale, independent Y ─────────────────────────────────────────
    visible = [axes[r, c] for r in range(N) for c in range(r, N)]
    xmin = min(ax.get_xlim()[0] for ax in visible)
    xmax = max(ax.get_xlim()[1] for ax in visible)
    for ax in visible:
        ax.set_xlim(xmin, xmax)

    # ── Column titles (top row) ────────────────────────────────────────────────
    for col in range(N):
        ct    = ct_list[col]
        color = _ct_color(ct, ct_list)
        axes[0, col].set_title(ct, fontsize=max(6, 8 - N // 4),
                               fontweight="bold", color=color)

    # ── Row labels in left margin ──────────────────────────────────────────────
    for row in range(N):
        ct    = ct_list[row]
        color = _ct_color(ct, ct_list)
        bbox  = axes[row, row].get_position()   # use diagonal (leftmost visible)
        y_ctr = (bbox.y0 + bbox.y1) / 2
        fig.text(
            0.01, y_ctr, ct,
            rotation=90, va="center", ha="center",
            fontsize=max(6, 8 - N // 4), fontweight="bold", color=color,
        )

    fig.subplots_adjust(left=0.12, right=0.97, top=0.90, bottom=0.07,
                        hspace=0.35, wspace=0.30)

    # ── Build CSV data ─────────────────────────────────────────────────────────
    rows: list[dict] = []
    for ct in ct_list:
        for val in dist_cache[(ct, ct)]["within_i"]:
            rows.append({"ct_i": ct, "ct_j": ct, "pair_type": "within_i", "distance": float(val)})
    for ct_i, ct_j in itertools.combinations(ct_list, 2):
        d = dist_cache[(ct_i, ct_j)]
        for val in d["within_i"]:
            rows.append({"ct_i": ct_i, "ct_j": ct_j, "pair_type": "within_i", "distance": float(val)})
        for val in d["cross"]:
            rows.append({"ct_i": ct_i, "ct_j": ct_j, "pair_type": "cross", "distance": float(val)})
        for val in d["within_j"]:
            rows.append({"ct_i": ct_i, "ct_j": ct_j, "pair_type": "within_j", "distance": float(val)})
    data_df = pd.DataFrame(rows)

    return fig, data_df


def plot_matrix_between(
    dist_cache: dict[tuple[str, str], dict],
    ct_list: list[str],
    metric: str,
) -> tuple[plt.Figure, pd.DataFrame]:
    """Upper-triangle N×N matrix — diagonal: within-type, upper: cross-type only KDE.

    Each upper-triangle panel shows a single filled KDE of the cross-type
    distances, colored with a blend of the two cell-type colors.  Lower
    triangle panels are hidden.

    Returns (figure, data_df) where data_df has columns:
        ct_i, ct_j, pair_type, distance
    """
    N = len(ct_list)

    metric_label = {
        "euclidean": "Euclidean distance",
        "cosine":    "Cosine distance",
        "manhattan": "Manhattan distance",
    }.get(metric, metric)

    panel_size = max(2.0, min(3.5, 18.0 / N))
    fig, axes  = plt.subplots(N, N, figsize=(N * panel_size, N * panel_size), squeeze=False)

    for row in range(N):
        for col in range(N):
            ax      = axes[row, col]
            ct_row  = ct_list[row]
            ct_col  = ct_list[col]
            color_r = _ct_color(ct_row, ct_list)
            color_c = _ct_color(ct_col, ct_list)
            color_blend = tuple((a + b) / 2 for a, b in zip(color_r, color_c))

            if row == col:
                # ── Diagonal: within-type KDE on gray background ───────────
                ax.set_facecolor("#f5f5f5")
                d = dist_cache[(ct_row, ct_row)]
                sns.kdeplot(
                    d["within_i"], ax=ax,
                    color=color_r, linewidth=2.0,
                    fill=True, alpha=0.35, warn_singular=False,
                )

            elif row < col:
                # ── Upper triangle: cross-type KDE ─────────────────────────
                ct_i = ct_list[row]
                ct_j = ct_list[col]
                d = dist_cache[(ct_i, ct_j)]
                sns.kdeplot(
                    d["cross"], ax=ax,
                    color=color_blend, linewidth=1.8,
                    fill=True, alpha=0.30, warn_singular=False,
                )

            else:
                # ── Lower triangle: hidden ─────────────────────────────────
                ax.set_visible(False)
                continue

            ax.tick_params(labelsize=max(5, 7 - N // 4))
            ax.set_xlabel(metric_label, fontsize=max(5, 7 - N // 4))
            ax.set_ylabel("Density", fontsize=max(5, 7 - N // 4))

    # ── Shared X scale, independent Y ─────────────────────────────────────────
    visible = [axes[r, c] for r in range(N) for c in range(r, N)]
    xmin = min(ax.get_xlim()[0] for ax in visible)
    xmax = max(ax.get_xlim()[1] for ax in visible)
    for ax in visible:
        ax.set_xlim(xmin, xmax)

    # ── Column titles (top row) ────────────────────────────────────────────────
    for col in range(N):
        ct    = ct_list[col]
        color = _ct_color(ct, ct_list)
        axes[0, col].set_title(ct, fontsize=max(6, 8 - N // 4),
                               fontweight="bold", color=color)

    # ── Row labels in left margin ──────────────────────────────────────────────
    for row in range(N):
        ct    = ct_list[row]
        color = _ct_color(ct, ct_list)
        bbox  = axes[row, row].get_position()   # use diagonal (leftmost visible)
        y_ctr = (bbox.y0 + bbox.y1) / 2
        fig.text(
            0.01, y_ctr, ct,
            rotation=90, va="center", ha="center",
            fontsize=max(6, 8 - N // 4), fontweight="bold", color=color,
        )

    fig.subplots_adjust(left=0.12, right=0.97, top=0.90, bottom=0.07,
                        hspace=0.35, wspace=0.30)

    # ── Build CSV data ─────────────────────────────────────────────────────────
    rows: list[dict] = []
    for ct in ct_list:
        for val in dist_cache[(ct, ct)]["within_i"]:
            rows.append({"ct_i": ct, "ct_j": ct, "pair_type": "within_i", "distance": float(val)})
    for ct_i, ct_j in itertools.combinations(ct_list, 2):
        d = dist_cache[(ct_i, ct_j)]
        for val in d["within_i"]:
            rows.append({"ct_i": ct_i, "ct_j": ct_j, "pair_type": "within_i", "distance": float(val)})
        for val in d["cross"]:
            rows.append({"ct_i": ct_i, "ct_j": ct_j, "pair_type": "cross", "distance": float(val)})
        for val in d["within_j"]:
            rows.append({"ct_i": ct_i, "ct_j": ct_j, "pair_type": "within_j", "distance": float(val)})
    data_df = pd.DataFrame(rows)

    return fig, data_df


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--h5ad",      required=True,
                   help="Path to AnnData .h5ad file")
    p.add_argument("--obs-col",   required=True,
                   help="obs column containing cell-type labels")
    p.add_argument("--obsm-key",  required=True,
                   help="obsm key for embedding (e.g. X_pca)")
    p.add_argument("--output",    required=True,
                   help="Output PDF path")
    p.add_argument("--n-pairs",   type=int, default=10_000,
                   help="Pairs sampled per sub-pairing category (default: 10000)")
    p.add_argument("--metric",    default="cosine",
                   choices=list(_SUPPORTED_METRICS),
                   help="Distance metric (default: cosine)")
    p.add_argument("--seed",      type=int, default=42,
                   help="Random seed (default: 42)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    # ── Load data ─────────────────────────────────────────────────────────────
    h5ad_path = Path(args.h5ad)
    if not h5ad_path.is_file():
        raise FileNotFoundError(f"h5ad file not found: {h5ad_path}")

    print(f"[1/8] Loading {h5ad_path} …")
    adata = anndata.read_h5ad(h5ad_path)

    # ── Validate inputs ───────────────────────────────────────────────────────
    if args.obs_col not in adata.obs.columns:
        raise ValueError(
            f"obs column '{args.obs_col}' not found. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    if args.obsm_key not in adata.obsm:
        raise ValueError(
            f"obsm key '{args.obsm_key}' not found. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    labels    = adata.obs[args.obs_col].astype(str).values
    embedding = np.array(adata.obsm[args.obsm_key], dtype=np.float32)

    ct_list = sorted(pd.unique(labels))
    n_types = len(ct_list)
    n_pairs_total = n_types * (n_types + 1) // 2

    print(f"  Cell types ({n_types}): {ct_list}")
    print(f"  Embedding shape: {embedding.shape}")
    print(f"  Panels to generate: {n_pairs_total}")
    print(f"  Metric: {args.metric}  |  n_pairs per category: {args.n_pairs:,}")

    if n_types < 2:
        raise ValueError(
            f"Need at least 2 unique values in '{args.obs_col}', "
            f"found {n_types}."
        )

    # ── Check output path ─────────────────────────────────────────────────────
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    metric_label = {"euclidean": "Euclidean distance",
                    "cosine":    "Cosine distance",
                    "manhattan": "Manhattan distance"}.get(args.metric, args.metric)

    # ── Grid figure ───────────────────────────────────────────────────────────
    print("[2/8] Sampling pairs and computing distances for grid …")
    rng = np.random.default_rng(args.seed)

    fig = plot_grid(
        embedding=embedding,
        obs_labels=labels,
        ct_list=ct_list,
        n_pairs=args.n_pairs,
        metric=args.metric,
        rng=rng,
        ncols=3,
    )

    fig.suptitle(
        f"Pairwise {metric_label} distributions  "
        f"(obsm: {args.obsm_key}, n_pairs={args.n_pairs:,})",
        fontsize=11, y=1.01,
    )

    print(f"[3/8] Saving grid to {out_path} …")
    fig.savefig(out_path, bbox_inches="tight", format="pdf")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    # ── Joint figure ──────────────────────────────────────────────────────────
    joint_base = out_path.parent / (out_path.stem + "_joint")
    joint_path = joint_base.with_suffix(".pdf")
    print(f"[4/8] Computing joint distributions …")
    rng2 = np.random.default_rng(args.seed)

    fig_j = plot_joint(
        embedding=embedding,
        obs_labels=labels,
        ct_list=ct_list,
        n_pairs=args.n_pairs,
        metric=args.metric,
        rng=rng2,
    )
    fig_j.suptitle(
        f"Joint pairwise {metric_label} distributions  "
        f"(obsm: {args.obsm_key}, n_pairs={args.n_pairs:,})",
        fontsize=11, y=1.01,
    )
    fig_j.savefig(joint_path, bbox_inches="tight", format="pdf")
    fig_j.savefig(joint_base.with_suffix(".png"), bbox_inches="tight", dpi=150)
    fig_j.savefig(joint_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig_j)

    # ── Shared distance cache for matrix figures ───────────────────────────────
    matrix_joint_base = out_path.parent / (out_path.stem + "_matrix_joint")
    matrix_btw_base   = out_path.parent / (out_path.stem + "_matrix_between")
    matrix_csv        = out_path.parent / (out_path.stem + "_matrix_data.csv")

    print("[5/8] Computing shared distance cache for matrix figures …")
    rng3 = np.random.default_rng(args.seed)
    dist_cache = _build_dist_cache(
        embedding=embedding,
        obs_labels=labels,
        ct_list=ct_list,
        n_pairs=args.n_pairs,
        metric=args.metric,
        rng=rng3,
    )

    # ── Matrix joint figure ───────────────────────────────────────────────────
    print(f"[6/8] Plotting matrix_joint …")
    fig_mj, data_df = plot_matrix_joint(dist_cache, ct_list, args.metric)
    fig_mj.suptitle(
        f"Joint pairwise {metric_label} matrix — pooled within + cross distributions\n"
        f"(obsm: {args.obsm_key}, n_pairs={args.n_pairs:,})",
        fontsize=10, y=0.995,
    )
    fig_mj.savefig(matrix_joint_base.with_suffix(".pdf"), bbox_inches="tight", format="pdf")
    fig_mj.savefig(matrix_joint_base.with_suffix(".png"), bbox_inches="tight", dpi=150)
    fig_mj.savefig(matrix_joint_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig_mj)

    # ── Matrix between figure ─────────────────────────────────────────────────
    print(f"[7/8] Plotting matrix_between …")
    fig_mb, _ = plot_matrix_between(dist_cache, ct_list, args.metric)
    fig_mb.suptitle(
        f"Between-type {metric_label} matrix — cross-type distances only\n"
        f"(obsm: {args.obsm_key}, n_pairs={args.n_pairs:,})",
        fontsize=10, y=0.995,
    )
    fig_mb.savefig(matrix_btw_base.with_suffix(".pdf"), bbox_inches="tight", format="pdf")
    fig_mb.savefig(matrix_btw_base.with_suffix(".png"), bbox_inches="tight", dpi=150)
    fig_mb.savefig(matrix_btw_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig_mb)

    print(f"[8/8] Saving matrix data CSV …")
    data_df.to_csv(matrix_csv, index=False)

    print("Done.")
    print(f"  Grid:           {out_path.resolve()}")
    print(f"  Joint:          {joint_path.resolve()}")
    print(f"  Matrix joint:   {matrix_joint_base.with_suffix('.pdf').resolve()}")
    print(f"  Matrix between: {matrix_btw_base.with_suffix('.pdf').resolve()}")
    print(f"  Matrix CSV:     {matrix_csv.resolve()}")


if __name__ == "__main__":
    main()
