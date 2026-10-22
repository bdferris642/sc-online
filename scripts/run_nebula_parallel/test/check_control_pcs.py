"""Verify control_pcs.py against a scanpy + sklearn reference and check the in-place h5ad update.

Run with a python that has scanpy + scikit-learn (the reference implementation; e.g. the
r-4.3.0 env), from run_nebula_parallel/test:
    python check_control_pcs.py ORIGINAL.h5ad UPDATED.h5ad MODEL_DIR N_HVGS
ORIGINAL = sim before control_pcs.py compute; UPDATED = the same file after it (written in place).

Reference = what clean_h5ad/compute_pcs.py does, but normalising by the full-row total (as
control_pcs.py does in both fit and projection): normalize_total(1e4) on control cells, log1p,
genes detected in > 1% of controls, scanpy HVG flavor='seurat', dense, scale (ddof=0, sd+1e-8),
TruncatedSVD(arpack, exact). Then donor-mean OLS checked with R lm() via the r-4.3.0 Rscript.
"""
import glob
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

orig, upd, model_dir, n_hvgs = sys.argv[1], sys.argv[2], Path(sys.argv[3]), int(sys.argv[4])
slogan = Path(upd).name[:-5]
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")


A = ad.read_h5ad(orig)
B = ad.read_h5ad(upd)
m = np.load(model_dir / f"{slogan}_control_pca.npz", allow_pickle=False)

# ---- reference PCA
ctrl = (A.obs["case_control"].astype(str) == "ctr").to_numpy()
C = ad.AnnData(sp.csr_matrix(A.layers["counts"][ctrl], dtype=np.float64), var=A.var.copy())
sc.pp.normalize_total(C, target_sum=1e4)
sc.pp.log1p(C)
det = np.asarray((C.X > 0).sum(axis=0)).ravel() / C.n_obs
C = C[:, det > 0.01].copy()
sc.pp.highly_variable_genes(C, flavor="seurat", n_top_genes=min(n_hvgs, C.n_vars))
ref_hvgs = C.var_names[C.var["highly_variable"].to_numpy()]
check("HVG set identical", set(ref_hvgs) == set(m["hvgs"]), f"(ref {len(ref_hvgs)}, ours {len(m['hvgs'])}, "
      f"sym diff {len(set(ref_hvgs) ^ set(m['hvgs']))})")
Xh = np.asarray(C[:, list(m["hvgs"])].X.todense())
mu, sd = Xh.mean(0), Xh.std(0) + 1e-8
check("scaling mean/sd", np.allclose(mu, m["mean"], atol=1e-10) and np.allclose(sd, m["sd"], atol=1e-10))
Z = (Xh - mu) / sd
k = m["loadings"].shape[1]
svd = TruncatedSVD(n_components=k, algorithm="arpack", random_state=0).fit(Z)
ref_load = svd.components_.T
corr = np.abs(np.sum(ref_load * m["loadings"], axis=0))   # unit vectors -> |cos|
check("loadings match up to sign", corr.min() > 0.9999, f"(min |cos| {corr.min():.6f})")
ref_ve = svd.singular_values_ ** 2 / (svd.singular_values_ ** 2).sum()
check("variance explained", np.allclose(ref_ve, m["variance_explained"], rtol=1e-6))

# ---- projection of ALL cells with the reference model (full-row normalisation)
Aall = ad.AnnData(sp.csr_matrix(A.layers["counts"], dtype=np.float64), var=A.var.copy())
sc.pp.normalize_total(Aall, target_sum=1e4)
sc.pp.log1p(Aall)
P_ref = ((np.asarray(Aall[:, list(m["hvgs"])].X.todense()) - mu) / sd) @ ref_load
names = [str(x) for x in m["pcs"]]
P = B.obs[names].to_numpy()
pc_corr = np.array([abs(np.corrcoef(P_ref[:, j], P[:, j])[0, 1]) for j in range(k)])
check("projected PCs match up to sign", pc_corr.min() > 0.9999, f"(min |r| {pc_corr.min():.6f})")
check("PC names zero-padded", names[:3] == ["control_PC_01", "control_PC_02", "control_PC_03"])

# ---- donor OLS vs R lm
res = pd.read_csv(model_dir / f"{slogan}_control_pc_test.csv")
obs = A.obs.copy()
obs["log10_nUMI"] = np.log10(np.asarray(A.layers["counts"].astype(np.float64).sum(1)).ravel() + 1)
d = obs.groupby("participant_id", observed=True)
tab = pd.DataFrame({"pd": d["case_control"].agg(lambda s: float((s == "pd").mean())),
                    "age": d["age"].mean(), "log10_nUMI": d["log10_nUMI"].mean(),
                    "frac_mito": d["frac_mito"].mean(), "frac_intronic": d["frac_intronic"].mean(),
                    "sexMale": d["sex"].agg(lambda s: float((s == "Male").mean())),
                    "chemv4": d["chemistry"].agg(lambda s: float((s == "v4").mean()))})
pcm = pd.DataFrame(P, index=obs.index, columns=names).groupby(obs["participant_id"].astype(str)).mean()
tab.index = tab.index.astype(str)
tab = tab.join(pcm)
with tempfile.TemporaryDirectory() as td:
    tab.to_csv(f"{td}/t.csv")
    r = ("d <- read.csv('%s/t.csv', row.names=1); pcs <- grep('^control_PC_', names(d), value=TRUE);"
         "o <- t(sapply(pcs, function(p) summary(lm(as.formula(paste(p, '~ pd + age + log10_nUMI + frac_mito + frac_intronic + sexMale + chemv4')), d))$coef['pd', c(1,2,4)]));"
         "write.csv(o, '%s/o.csv')" % (td, td))
    subprocess.run(["/home/ferris/anaconda3/envs/r-4.3.0/bin/Rscript", "-e", r], check=True, capture_output=True)
    o = pd.read_csv(f"{td}/o.csv", index_col=0)
o.columns = ["beta", "se", "p"]
for c in ["beta", "se", "p"]:
    check(f"OLS {c} = R lm (sign may flip for beta)", np.allclose(np.abs(o[c].to_numpy()) if c == "beta" else o[c].to_numpy(),
          np.abs(res[c].to_numpy()) if c == "beta" else res[c].to_numpy(), rtol=1e-6, atol=1e-12))

# ---- in-place write
other = [c for c in A.obs.columns if not c.startswith("control_PC")]
check("non-PC obs columns identical", all(A.obs[c].equals(B.obs[c]) for c in other) and A.obs.index.equals(B.obs.index))
check("old control_PC1..30 removed", not any(c.startswith("control_PC") and not c.startswith("control_PC_") for c in B.obs.columns))
check("layers identical", all((A.layers[l] != B.layers[l]).nnz == 0 for l in A.layers) and (A.X != B.X).nnz == 0)
check("other uns kept", np.array_equal(A.uns["keep_me"]["value"], B.uns["keep_me"]["value"]))
cp = B.uns["control_pca"]
check("uns/control_pca new keys", {"pcs", "elbow_pcs", "significant_pcs", "created", "case_control_test", "model_path"} <= set(cp))
bk = sorted(glob.glob(str(Path(upd).parent / "backup_*" / f"{slogan}_obs_uns.h5")))[-1]
with h5py.File(bk) as h:
    old_obs = ad.io.read_elem(h["obs"])
    old_cp = ad.io.read_elem(h["uns_control_pca"])
check("backup holds original obs", old_obs.equals(A.obs))
check("backup holds original uns/control_pca", sorted(old_cp) == sorted(A.uns["control_pca"]))
print("ALL PASS" if ok else "SOME CHECKS FAILED")
sys.exit(0 if ok else 1)
