"""
go_enrichment.py
----------------
GO Biological Process and Molecular Function enrichment for eGenes.
Functions _annotated_genes_from_gmt, _load_gmt_filtered, _run_enrich adapted
from enrichment-analysis/modules/go_enrichment.py.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Dict, List

import pandas as pd
import gseapy as gp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

log = logging.getLogger(__name__)


# ── Three functions copied verbatim from enrichment-analysis/modules/go_enrichment.py ──

def _annotated_genes_from_gmt(gmt_path: Path) -> set:
    """
    Return all gene symbols that appear in at least one term in the GMT file,
    regardless of term size.  Used to define the 'annotated universe' for
    background filtering (equivalent to goseq use_genes_without_cat=FALSE).
    """
    genes: set = set()
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            genes.update(g.strip() for g in parts[2:] if g and g.strip())
    return genes


def _load_gmt_filtered(gmt_path: Path, min_size: int, max_size: int) -> Dict[str, List[str]]:
    """
    Load a GMT file and return only terms whose total gene count falls in
    [min_size, max_size].  This ensures the size filter applies to the
    canonical term size rather than the (background-dependent) overlap count.
    """
    terms: Dict[str, List[str]] = {}
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term_name = parts[0]
            genes = [g.strip() for g in parts[2:] if g and g.strip()]
            if min_size <= len(genes) <= max_size:
                terms[term_name] = genes
    return terms


def _run_enrich(
    gene_list:   List[str],
    background:  List[str],
    gene_sets,
    label:       str,
    out_dir:     Path,
    min_overlap: int = 3,
) -> pd.DataFrame | None:
    import logging as _logging
    _gseapy_log = _logging.getLogger("gseapy")
    _prev_level  = _gseapy_log.level
    _gseapy_log.setLevel(_logging.CRITICAL)
    try:
        enr = gp.enrich(
            gene_list  = gene_list,
            gene_sets  = gene_sets,
            background = background,
            outdir     = None,
            verbose    = False,
            no_plot    = True,
        )
        res = enr.results
        if res is None or not isinstance(res, pd.DataFrame) or res.empty:
            log.info("  [%s] No results.", label)
            return None
        if "Adjusted P-value" in res.columns:
            res = res.sort_values("Adjusted P-value")
        elif "P-value" in res.columns:
            res = res.sort_values("P-value")
        if "Overlap" in res.columns and min_overlap > 1:
            overlap_k = res["Overlap"].str.split("/").str[0].astype(int)
            res = res[overlap_k >= min_overlap].reset_index(drop=True)
        if res.empty:
            return None
        res.to_csv(out_dir / f"go_enrichment_{label}.csv", index=False)
        return res
    except Exception as exc:
        log.error("  [%s] Enrichment failed: %s", label, exc)
        return None
    finally:
        _gseapy_log.setLevel(_prev_level)


# ── Main run function ──────────────────────────────────────────────────────────

def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None,
        go_bp_gmt=None, go_mf_gmt=None,
        padj_thresh: float = 0.05,
        min_size: int = 10, max_size: int = 200, min_overlap: int = 3,
        **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if go_bp_gmt is None or go_mf_gmt is None:
        print("  [go_enrichment] GMT paths not provided — skipping.")
        return

    go_bp_gmt, go_mf_gmt = Path(go_bp_gmt), Path(go_mf_gmt)
    if not go_bp_gmt.exists() or not go_mf_gmt.exists():
        print(f"  [go_enrichment] GMT files not found — skipping.")
        return

    padj_gene = df.get("padj_gene", pd.Series(np.nan, index=df.index))
    sig_df = df[padj_gene < padj_thresh].drop_duplicates(subset=["ensg_id"])
    all_df = df.drop_duplicates(subset=["ensg_id"])

    # Use gene_symbol column for GMT lookup (display-only); ensg_id identifies the gene
    def _to_sym(row_df: "pd.DataFrame") -> "List[str]":
        """Return symbol list from gene_symbol col (fallback: gene_loc_v2, then skip)."""
        if "gene_symbol" in row_df.columns:
            syms = row_df["gene_symbol"].astype(str).tolist()
        elif gene_loc_df is not None:
            sym_map = dict(zip(gene_loc_df["ensg_id"].astype(str),
                               gene_loc_df["gene_symbol"].astype(str)))
            syms = [sym_map.get(str(g), "nan") for g in row_df["ensg_id"].astype(str)]
        else:
            syms = ["nan"] * len(row_df)
        return [s for s in syms if s not in ("nan", "None", "")]

    sig_syms = _to_sym(sig_df)
    bg_syms  = _to_sym(all_df)

    if len(sig_syms) < 5:
        print(f"  [go_enrichment] Only {len(sig_syms)} sig gene symbols — skipping.")
        return

    print(f"  [go_enrichment] {len(sig_syms)} sig genes, {len(bg_syms)} background genes")

    results = {}
    for label, gmt_path in [("BP", go_bp_gmt), ("MF", go_mf_gmt)]:
        ann_genes = _annotated_genes_from_gmt(gmt_path)
        bg_filtered = [g for g in bg_syms if g in ann_genes]
        gene_sets = _load_gmt_filtered(gmt_path, min_size, max_size)
        res = _run_enrich(sig_syms, bg_filtered, gene_sets, label, out_dir, min_overlap)
        results[label] = res

    # Lollipop plot: top 20 terms combined
    frames = []
    for label, res in results.items():
        if res is not None and not res.empty:
            sub = res.head(20).copy()
            sub["namespace"] = label
            frames.append(sub)

    if not frames:
        print("  [go_enrichment] No enrichment results to plot.")
        return

    plot_df = pd.concat(frames, ignore_index=True)
    padj_col = "Adjusted P-value" if "Adjusted P-value" in plot_df.columns else "P-value"
    plot_df["-log10_padj"] = -np.log10(plot_df[padj_col].clip(1e-300))
    plot_df = plot_df.sort_values("-log10_padj", ascending=True).tail(30)

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.35)))
    colors = {"BP": "steelblue", "MF": "darkorange"}
    for _, row in plot_df.iterrows():
        ns = row.get("namespace", "BP")
        ax.plot([0, row["-log10_padj"]], [row["Term"], row["Term"]], color=colors.get(ns, "gray"), lw=2)
        ax.scatter([row["-log10_padj"]], [row["Term"]], color=colors.get(ns, "gray"), s=60, zorder=3)
    ax.axvline(-np.log10(0.05), color="red", linestyle="--", lw=1, label="adj-p=0.05")
    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor=c, label=ns) for ns, c in colors.items()]
    ax.legend(handles=legend_els + [plt.Line2D([0], [0], color="red", linestyle="--", label="FDR=0.05")])
    ax.set_xlabel("−log₁₀(adjusted p-value)")
    ax.set_title("GO enrichment — top terms (BP + MF)")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"go_enrichment.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [go_enrichment] Outputs → {out_dir}")
