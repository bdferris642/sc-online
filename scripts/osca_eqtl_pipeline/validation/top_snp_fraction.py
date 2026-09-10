"""
top_snp_fraction.py
-------------------
For each eGene, rank all SNPs by p-value and examine the rank of the
top significant SNP. Well-calibrated FDR → most sig SNPs should be rank 1.
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
    sig_mask = padj_snp < padj_thresh

    # Build ensg_id → gene_symbol map for output
    sym_map: dict = {}
    if "gene_symbol" in df.columns:
        sym_map = df.drop_duplicates("ensg_id").set_index("ensg_id")["gene_symbol"].to_dict()

    records = []
    for gene, grp in df.groupby("ensg_id"):
        n_sig = int(sig_mask.reindex(grp.index).fillna(False).sum())
        if n_sig == 0:
            continue
        # Rank all SNPs in this gene by p (ascending → rank 1 = best p)
        ranked = grp["p"].rank(method="min", ascending=True)
        sig_ranks = ranked[sig_mask.reindex(grp.index).fillna(False)].values.astype(int)
        records.append({
            "ensg_id": gene,
            "gene_symbol": sym_map.get(gene),
            "n_sig_snps": n_sig,
            "rank_of_top_sig_snp": int(sig_ranks.min()),
            "all_sig_ranks": ",".join(map(str, sorted(sig_ranks))),
        })

    csv_df = pd.DataFrame(records)
    csv_df.to_csv(out_dir / "top_snp_fraction.csv", index=False)

    if csv_df.empty:
        print("  [top_snp_fraction] No significant SNPs found.")
        return

    ranks = csv_df["rank_of_top_sig_snp"].values
    frac_rank1 = float(np.mean(ranks == 1))

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))
    max_rank = min(ranks.max(), 50)
    bins = np.arange(0.5, max_rank + 1.5, 1)
    ax.hist(ranks[ranks <= max_rank], bins=bins, color="steelblue", edgecolor="black", alpha=0.8)
    ax.set_xlabel("Rank of top significant SNP within gene (lower = closer to top)")
    ax.set_ylabel("Count (eGenes)")
    ax.set_title(f"Top significant SNP rank per eGene\n(rank-1 fraction = {frac_rank1:.1%}; showing ≤50)")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"top_snp_fraction.{ext}", dpi=150)
    plt.close(fig)
    print(f"  [top_snp_fraction] {len(csv_df)} eGenes, rank-1 fraction={frac_rank1:.1%} → {out_dir}")
