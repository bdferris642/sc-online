#!/usr/bin/env python
"""Partial-residual rank purity per gene, from NEBULA results + pseudobulk_h5ad.py outputs.

Port of compute_rank_purity() in ~/sc-online/add_rank_purity.R (groupby = donor, as called by
combine-nebula-results.R), with these deliberate differences:

1. Sign check: the purity is zeroed when it disagrees in sign with the contrast logFC. Upstream
   looked up the coefficient column named exactly `contrast_col` (e.g. "case_control"), which does
   not exist for a categorical contrast (NEBULA names it "case_controlpd"), so the check never ran.
   Here the column is contrast_col + case_val for a categorical contrast (contrast_col for a
   numeric one), and its absence is an error.
2. Categorical covariates that vary within a donor (e.g. chemistry) enter the design as the
   donor's fraction of cells at that level instead of the first cell's level. Pseudobulk = cell
   mean of log-expression and the fitted covariate effect is linear, so this is the exact cell
   mean of the fitted covariate contribution; donor-constant categoricals are unchanged.
3. Donors with a missing covariate or contrast value are dropped (logged) instead of being
   ranked as NA; genes with a missing coefficient get rank_purity = NaN. Donors with fewer than
   --min-cells-per-donor cells are dropped, matching the donors NEBULA was fit on.
4. Covariates are matched to coefficient columns by type (numeric: exact name; categorical:
   name + observed level) instead of by prefix, which upstream rejects as ambiguous once
   covariates share a prefix (control_PC1 vs control_PC10).

Steps: residual = pseudobulk - B X^T (B = genes x covariate-coefficients, X = donors x covariate
values), minus the per-gene median; then purity = Mann-Whitney U scaled to [-1, 1] (categorical
contrast; U counts case > control strictly, average ranks for ties) or Spearman rho (numeric
contrast); sign check; abs().
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from h5ad_io import ensure_env_python  # noqa: E402


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_coefs(combined_csv):
    """gene + logFC_* columns of the combined NEBULA summary, with the logFC_ prefix stripped."""
    summ = pd.read_csv(combined_csv, index_col=0)
    if "gene" not in summ.columns:
        raise KeyError(f"'gene' column missing in {combined_csv}")
    lfc = [c for c in summ.columns if "logFC" in c]
    coef = summ[["gene", *lfc]].copy()
    coef.columns = ["gene"] + [c[len("logFC_"):] if c.startswith("logFC_") else c for c in lfc]
    if coef["gene"].duplicated().any():
        raise ValueError("duplicated genes in NEBULA summary")
    return coef


def match_covariates(covariates, coef_cols, contrast_coef_col, types, levels):
    """Type-aware matching: {cov: [coef columns]}.

    Numeric covariate -> the column named exactly `cov`. Categorical -> columns `cov + level` for
    levels observed in the data (the reference level has none). Upstream matched by prefix, which
    breaks once covariates share a prefix (control_PC1 vs control_PC10..PC19). Every coefficient
    column other than the intercept and the contrast must be claimed by exactly one covariate.
    """
    out, claimed = {}, {}
    for cov in covariates:
        if cov not in types:
            raise KeyError(f"covariate '{cov}' missing from pseudobulk outputs (rerun pseudobulk_h5ad.py)")
        if types[cov] == "numeric":
            matched = [cov] if cov in coef_cols else []
        else:
            lv = set(levels.loc[levels.covariate == cov, "level"])
            matched = [c for c in coef_cols if c.startswith(cov) and c[len(cov):] in lv]
        if not matched:
            raise KeyError(f"No coefficient column found for {types[cov]} covariate '{cov}'. "
                           f"Columns: {list(coef_cols)}")
        for c in matched:
            if c in claimed:
                raise ValueError(f"Coefficient column '{c}' matches both '{claimed[c]}' and '{cov}'")
            claimed[c] = cov
        out[cov] = matched
    unclaimed = [c for c in coef_cols if c not in claimed and c != contrast_coef_col and c != "(Intercept)"]
    if unclaimed:
        raise ValueError(f"Coefficient columns not assigned to any covariate: {unclaimed}")
    return out


def build_design(matches, types, donors, numeric, levels):
    """Return (column names, X donors x k) following .expand_covariate_design with fractional dummies."""
    names, cols = [], []
    for cov, matched in matches.items():
        if cov not in types:
            raise KeyError(f"covariate '{cov}' missing from pseudobulk outputs (rerun pseudobulk_h5ad.py)")
        for col in matched:
            suffix = col[len(cov):]
            if suffix == "":
                if types[cov] != "numeric":
                    raise ValueError(f"'{cov}' has an exact-match coefficient column but is categorical; "
                                     f"categoricals need dummy columns like '{cov}<level>'")
                x = numeric.set_index("donor").loc[donors, cov].to_numpy(dtype=float)
                log(f"  [numeric]  {col:<28} <- donor mean of {cov}")
            else:
                # match_covariates only assigns cov+level columns for observed levels
                lv = levels[(levels.covariate == cov)]
                piv = lv.pivot(index="donor", columns="level", values="fraction").reindex(donors)
                x = piv[suffix].to_numpy(dtype=float)
                log(f"  [dummy]    {col:<28} <- fraction of cells with {cov} == '{suffix}'")
            names.append(col)
            cols.append(x)
    if not cols:
        raise ValueError("No covariate columns could be matched")
    return names, np.column_stack(cols)


def purity_categorical(res, is_case):
    """Mann-Whitney U (case > control, strict) scaled to [-1, 1], vectorised over rows."""
    n_genes, _ = res.shape
    ranks = rankdata(res, axis=1)  # average ties, as R rank()
    case, ctrl = ranks[:, is_case], ranks[:, ~is_case]
    big = res.shape[1] + 1.0
    offs = np.arange(n_genes)[:, None] * big
    sorted_ctrl = np.sort(ctrl + offs, axis=1).ravel()
    pos = np.searchsorted(sorted_ctrl, (case + offs).ravel(), side="left").reshape(case.shape)
    u = (pos - np.arange(n_genes)[:, None] * ctrl.shape[1]).sum(axis=1)
    return 2.0 * u / (case.shape[1] * ctrl.shape[1]) - 1.0


def purity_spearman(res, y):
    """Spearman correlation of each row of res with y."""
    rx = rankdata(res, axis=1)
    ry = rankdata(y)
    rx = rx - rx.mean(axis=1, keepdims=True)
    ry = ry - ry.mean()
    return (rx @ ry) / (np.sqrt((rx ** 2).sum(axis=1)) * np.sqrt((ry ** 2).sum()))


def main():
    ensure_env_python()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combined-csv", required=True, help="combined NEBULA summary CSV")
    ap.add_argument("--pb-dir", required=True, help="pseudobulk_h5ad.py output dir")
    ap.add_argument("--covs", default="", help="same --covs string given to NEBULA")
    ap.add_argument("--contrast-col", required=True)
    ap.add_argument("--case-val", default=None, help="case level for a categorical contrast")
    ap.add_argument("--min-cells-per-donor", type=int, default=5,
                    help="drop donors with fewer cells, as run-nebula.R's --min-num-cells-per-id (default 5)")
    ap.add_argument("--out", required=True, help="output TSV (gene, rank_purity)")
    a = ap.parse_args()

    pb_dir = Path(a.pb_dir)
    pb = np.load(pb_dir / "pb_mean.npy")
    donors = np.loadtxt(pb_dir / "donors.txt", dtype=str, delimiter="\t", ndmin=1)
    pb_genes = np.loadtxt(pb_dir / "genes.txt", dtype=str, delimiter="\t", ndmin=1)
    types = json.load(open(pb_dir / "covariate_types.json"))
    numeric = pd.read_csv(pb_dir / "donor_numeric.tsv", sep="\t", dtype={"donor": str})
    levels = pd.read_csv(pb_dir / "donor_levels.tsv", sep="\t", dtype={"donor": str, "level": str, "covariate": str})

    coef = load_coefs(a.combined_csv)
    covariates = [c for c in a.covs.split(",") if c and c != a.contrast_col]
    coef_cols = [c for c in coef.columns if c != "gene"]

    # ---- contrast values per donor
    if a.contrast_col not in types:
        raise KeyError(f"contrast '{a.contrast_col}' missing from pseudobulk outputs")
    if types[a.contrast_col] == "categorical":
        if a.case_val is None:
            raise ValueError(f"contrast '{a.contrast_col}' is categorical; --case-val is required")
        lv = levels[levels.covariate == a.contrast_col]
        if a.case_val not in set(lv.level):
            raise ValueError(f"case value '{a.case_val}' not in '{a.contrast_col}': {sorted(set(lv.level))}")
        y = lv.pivot(index="donor", columns="level", values="fraction").reindex(donors)[a.case_val].to_numpy()
        contrast_coef_col = f"{a.contrast_col}{a.case_val}"
    else:
        y = numeric.set_index("donor").loc[donors, a.contrast_col].to_numpy(dtype=float)
        contrast_coef_col = a.contrast_col
    if contrast_coef_col not in coef.columns:
        raise KeyError(f"contrast coefficient column '{contrast_coef_col}' not in NEBULA summary: {coef_cols}")
    matches = match_covariates(covariates, coef_cols, contrast_coef_col, types, levels)

    # ---- genes
    gidx = pd.Index(pb_genes).get_indexer(coef["gene"])
    keep = gidx >= 0
    log(f"{int(keep.sum())} / {len(coef)} genes found in pseudobulk")
    if not keep.any():
        raise ValueError("No genes in common between NEBULA results and the h5ad")
    coef = coef[keep].reset_index(drop=True)
    pb_sub = pb[:, gidx[keep]].T  # genes x donors

    # ---- design and donors with complete data
    log("Design:")
    names, x = build_design(matches, types, donors, numeric, levels)
    log(f"  formula: ~ {' + '.join(names)}")
    n_cells = pd.read_csv(pb_dir / "donor_n_cells.tsv", sep="\t", dtype={"donor": str}).set_index("donor")
    few = n_cells.loc[donors, "n_cells"].to_numpy() < a.min_cells_per_donor
    if few.any():
        log(f"dropping {int(few.sum())} donors with < {a.min_cells_per_donor} cells (as run-nebula.R): "
            f"{list(donors[few][:10])}")
    ok = np.isfinite(x).all(axis=1) & np.isfinite(y)
    if not ok.all():
        log(f"dropping {int((~ok).sum())} donors with missing covariate/contrast values: {list(donors[~ok][:10])}")
    ok &= ~few
    x, y, pb_sub = x[ok], y[ok], pb_sub[:, ok]
    b = coef[names].to_numpy(dtype=float)

    res = pb_sub - b @ x.T
    res = res - np.nanmedian(res, axis=1, keepdims=True)
    bad_gene = ~np.isfinite(res).all(axis=1)
    if bad_gene.any():
        log(f"{int(bad_gene.sum())} genes have missing coefficients/expression -> rank_purity NaN")

    purity = np.full(len(coef), np.nan)
    good = ~bad_gene
    if types[a.contrast_col] == "categorical":
        is_case = y == 1.0
        if is_case.all() or not is_case.any():
            raise ValueError("contrast has only one group after filtering donors")
        log(f"categorical contrast: {int(is_case.sum())} case / {int((~is_case).sum())} control donors")
        purity[good] = purity_categorical(res[good], is_case)
    else:
        log(f"numeric contrast: Spearman over {len(y)} donors")
        purity[good] = purity_spearman(res[good], y)

    lfc = coef[contrast_coef_col].to_numpy(dtype=float)
    purity[purity * lfc < 0] = 0.0
    purity = np.abs(purity)

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": coef["gene"], "rank_purity": purity}).to_csv(out, sep="\t", index=False, na_rep="NA")
    log(f"wrote {out} (median rank_purity {np.nanmedian(purity):.3f})")


if __name__ == "__main__":
    main()
