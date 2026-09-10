"""
pi1_replication.py
------------------
Estimate π₁ (proportion of true positives) among sn-vta eQTLs that match
GTEx v8 Brain_Substantia_nigra significant pairs using Storey's method.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ENSG_TO_SYMBOL = Path("/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv")


def _storey_pi1(pvals: np.ndarray, n_bootstrap: int = 1000) -> tuple[float, float, float]:
    """Return (pi1, ci_lo, ci_hi) using lambda=0.9 conservative estimate + bootstrap CI."""
    pvals = np.asarray(pvals)
    lambdas = np.arange(0.05, 0.95, 0.05)
    pi0s = [np.mean(pvals > lam) / (1 - lam) for lam in lambdas]
    pi0_est = min(pi0s[-1], 1.0)
    pi1 = 1 - pi0_est

    # Bootstrap CI
    rng = np.random.default_rng(42)
    boot_pi1s = []
    for _ in range(n_bootstrap):
        boot_p = rng.choice(pvals, size=len(pvals), replace=True)
        b_pi0s = [np.mean(boot_p > lam) / (1 - lam) for lam in lambdas]
        boot_pi1s.append(1 - min(b_pi0s[-1], 1.0))
    ci_lo, ci_hi = float(np.percentile(boot_pi1s, 2.5)), float(np.percentile(boot_pi1s, 97.5))
    return float(pi1), ci_lo, ci_hi


def run(df: pd.DataFrame, out_dir: Path, gene_loc_df=None, gtex_sn=None, **kwargs) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if gtex_sn is None or not Path(gtex_sn).exists():
        print(f"  [pi1_replication] GTEx SN file not found — skipping.")
        return

    # ── Load GTEx SN significant pairs ────────────────────────────────────────
    print(f"  [pi1_replication] Loading GTEx SN pairs from {gtex_sn} …")
    gtex = pd.read_csv(gtex_sn, sep="\t")
    # variant_id format: chr_pos_ref_alt_b38
    parts = gtex["variant_id"].str.split("_", expand=True)
    gtex["chrom"] = parts[0].str.replace("chr", "", regex=False)
    gtex["pos"] = pd.to_numeric(parts[1], errors="coerce")

    # ── Ensembl → symbol → Entrez bridge ──────────────────────────────────────
    ensg_sym = {}
    if ENSG_TO_SYMBOL.exists():
        es = pd.read_csv(ENSG_TO_SYMBOL)
        ensg_sym = dict(zip(es["ensembl_gene_id"], es["hgnc_symbol"]))

    sym_entrez = {}
    if gene_loc_df is not None:
        # gene_loc named cols: entrez, chr, TSS, symbol, strand
        for _, row in gene_loc_df.iterrows():
            sym_entrez[str(row["symbol"])] = str(row["entrez"])

    # Map GTEx gene_id (ENSG.version) → Entrez
    gtex["ensg_base"] = gtex["gene_id"].str.split(".").str[0]
    gtex["symbol"] = gtex["ensg_base"].map(ensg_sym)
    gtex["entrez"] = gtex["symbol"].map(sym_entrez)
    gtex_matched = gtex.dropna(subset=["entrez", "pos"])

    # ── Match to sn-vta eQTLs ─────────────────────────────────────────────────
    df_match = df.copy()
    df_match["Gene"] = df_match["Gene"].astype(str)
    df_match["Chr"] = df_match["Chr"].astype(str)
    df_match["BP"] = pd.to_numeric(df_match["BP"], errors="coerce")

    matched_rows = []
    for _, grow in gtex_matched.iterrows():
        ent = str(grow["entrez"])
        chrom = str(grow["chrom"])
        pos = grow["pos"]
        hits = df_match[
            (df_match["Gene"] == ent) &
            (df_match["Chr"] == chrom) &
            ((df_match["BP"] - pos).abs() <= 10)
        ]
        for _, hrow in hits.iterrows():
            matched_rows.append({
                "entrez": ent,
                "chr": chrom,
                "pos_gtex": pos,
                "BP_sn": hrow["BP"],
                "p_gtex": grow["pval_nominal"],
                "p_sn": hrow["p"],
                "SNP": hrow["SNP"],
            })

    if not matched_rows:
        print("  [pi1_replication] No matched SNP-gene pairs found — skipping π₁.")
        return

    match_df = pd.DataFrame(matched_rows).drop_duplicates(subset=["entrez", "chr", "pos_gtex"])
    match_df.to_csv(out_dir / "pi1_gtex_sn.csv", index=False)
    n_matched = len(match_df)
    print(f"  [pi1_replication] {n_matched} matched pairs for π₁ estimation.")

    pvals = match_df["p_sn"].values
    pvals = np.clip(pvals, 1e-300, 1.0)
    pi1, ci_lo, ci_hi = _storey_pi1(pvals)
    print(f"  [pi1_replication] π₁ = {pi1:.3f} [{ci_lo:.3f}, {ci_hi:.3f}], n={n_matched}")

    # ── Plot ───────────────────────────────────────────────────────────────────
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar([0], [pi1], color="steelblue", width=0.5, zorder=3)
    ax.errorbar([0], [pi1], yerr=[[pi1 - ci_lo], [ci_hi - pi1]],
                fmt="none", color="black", capsize=6, linewidth=2)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, min(1.05, ci_hi + 0.1))
    ax.set_xticks([0])
    ax.set_xticklabels(["sn-vta vs\nGTEx SN"])
    ax.set_ylabel("π₁ (replication rate)")
    ax.set_title(f"π₁ = {pi1:.3f}  (n={n_matched} matched pairs)\n95% CI [{ci_lo:.3f}, {ci_hi:.3f}]")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"pi1_gtex_sn.{ext}", dpi=150)
    plt.close(fig)
