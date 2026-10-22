"""Simulate a small h5ad (CSR counts + norm_log) for testing the run_nebula_parallel workflow.

Donors: 24 (12 pd / 12 ctr); covariates sex, age (numeric), brain_bank (numeric-looking
categorical levels "1"/"2"), chemistry (varies within some donors). ~150 cells/donor, 400 genes,
a subset of genes shifted in pd and by chemistry.

Usage: simulate_h5ad.py OUT.h5ad [--pc-test]
  --pc-test (for the control-PC step): 600 genes, 40 donors, three latent cell-state programs
  (program 1 shifted in pd), obs frac_mito / frac_intronic / log10_nUMI, plus OLD-format
  control_PC1..30 obs columns, an old-format uns/control_pca and an unrelated uns/keep_me entry,
  so the in-place PC replacement can be checked.
"""
import sys
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

out = sys.argv[1]
pc_test = "--pc-test" in sys.argv[2:]
rng = np.random.default_rng(7)
n_donors, n_genes = (40, 600) if pc_test else (24, 400)
rows = []
for d in range(n_donors):
    n = rng.integers(100, 200)
    cc = "pd" if d % 2 else "ctr"
    sex = "Male" if d % 3 else "Female"
    age = float(rng.integers(55, 90))
    bank = str(1 + d % 2 if d < 20 else 2)
    chem_p = 0.0 if d < 8 else (1.0 if d < 16 else 0.4)  # donors 16-23 are mixed
    chem = np.where(rng.random(n) < chem_p, "v4", "v3")
    for k in range(n):
        rows.append((f"donor{d:02d}", cc, sex, age, bank, chem[k]))
obs = pd.DataFrame(rows, columns=["participant_id", "case_control", "sex", "age", "brain_bank", "chemistry"])
for c in ["participant_id", "case_control", "sex", "brain_bank", "chemistry"]:
    obs[c] = obs[c].astype("category")
obs.index = [str(i) for i in range(len(obs))]
n_cells = len(obs)

base = rng.lognormal(-1.0, 1.5, n_genes)
lib = rng.lognormal(7.5, 0.4, n_cells)
donor_eff = rng.normal(0, 0.3, (n_donors, n_genes))
pd_eff = np.zeros(n_genes); pd_eff[:40] = rng.choice([-0.8, 0.8], 40)
chem_eff = np.zeros(n_genes); chem_eff[40:80] = 0.7
d_idx = obs.participant_id.cat.codes.to_numpy()
mu = (lib[:, None] / lib.mean()) * base[None, :] * np.exp(
    donor_eff[d_idx] + (obs.case_control == "pd").to_numpy()[:, None] * pd_eff
    + (obs.chemistry == "v4").to_numpy()[:, None] * chem_eff)
if pc_test:
    # three latent cell-state programs on gene blocks; program 1 is shifted in pd donors
    z = rng.normal(0, 1, (n_cells, 3))
    z[:, 0] += 0.8 * (obs.case_control == "pd").to_numpy()
    load = np.zeros((3, n_genes))
    for k, (a, b) in enumerate([(100, 180), (180, 260), (260, 340)]):
        load[k, a:b] = 0.6
    mu = mu * np.exp(z @ load)
counts = rng.negative_binomial(5, 5 / (5 + mu)).astype(np.float32)
counts = sp.csr_matrix(counts)
cpm = sp.diags(1e4 / np.asarray(counts.sum(1)).ravel()) @ counts
norm_log = cpm.copy(); norm_log.data = np.log1p(norm_log.data).astype(np.float32)
norm_log = sp.csr_matrix(norm_log, dtype=np.float32)
var = pd.DataFrame(index=[f"ENSG{i:011d}" for i in range(n_genes)])
if pc_test:
    tot = np.asarray(counts.sum(1)).ravel()
    obs["log10_nUMI"] = np.log10(tot + 1)
    obs["frac_mito"] = rng.beta(2, 60, n_cells)
    obs["frac_intronic"] = rng.beta(30, 20, n_cells)
    for i in range(30):                      # stale, old-format PCs that must be replaced
        obs[f"control_PC{i + 1}"] = rng.normal(0, 1, n_cells).astype(np.float32)
adata = ad.AnnData(X=norm_log.copy(), obs=obs, var=var, layers={"counts": counts, "norm_log": norm_log})
if pc_test:
    adata.uns["control_pca"] = {"n_elbow_pcs": 3, "elbow_pcs": np.array(["control_PC1", "control_PC2", "control_PC3"]),
                                "significant_pcs": np.array(["control_PC5"])}
    adata.uns["keep_me"] = {"value": np.arange(5)}
adata.write_h5ad(out)
print(adata, "nnz", counts.nnz)
