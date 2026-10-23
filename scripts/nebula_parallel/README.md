# run_nebula_parallel — gene-chunked NEBULA from an h5ad, with control PCs and rank purity

A portable, self-activating workflow. Every script it runs lives in this directory and executes
under its own micromamba env (`run-nebula-parallel-env`). Nothing is called from `~/sc-online` or
from other sandbox workflows at run time. It replaces `nebula_h5ad/` and the manual
"compute control PCs → test them → type PCs into `--covs`" routine.

## Setup (once)
```bash
bash /home/ferris/cc-sandbox/run_nebula_parallel/setup.sh
```
- Creates `<sandbox>/micromamba_root/envs/run-nebula-parallel-env`:
  - R 4.5: Seurat 5.5, qs 0.27, Matrix, dplyr, getopt, glue, future, data.table, jsonlite;
  - Python 3.12: numpy, scipy, pandas, h5py, anndata;
  - GNU parallel.
- **nebula is pinned to 1.5.5**, the version the existing sn-vta DE results were fit with. That
  version was never released on CRAN (1.5.3 → 1.5.6), so it is built from the GitHub tag `v1.5.5`,
  and setup checks its DESCRIPTION version.
- Ends with R and Python verification blocks.
- **R isolation:** `setup.sh` and both runners export `R_PROFILE_USER=/dev/null` and
  `R_ENVIRON_USER=/dev/null` and unset `R_LIBS`/`R_LIBS_USER`. A personal `~/.Rprofile` that
  loads packages missing from this env (e.g. `library(repr)` for Jupyter) would otherwise abort
  every R start. If you call an R script directly, prefix it with `R_PROFILE_USER=/dev/null`.

## Usage
```bash
/home/ferris/cc-sandbox/run_nebula_parallel/run_nebula_parallel.sh \
  --path /mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/sn-vta-immune-latest.h5ad \
  --outdir /path/to/outdir \
  --covs case_control,age,sex,frac_mito,frac_intronic,region_SN,region_VTA,brain_bank,chemistry \
  --contrast-col case_control --case-val pd \
  --compute-control-pcs 1 --pcs-to-add elbow significant --pc-max 20 \
  --get-rank-purity 1 \
  --n-folds 60 --n-cores 60 --suffix bb_chemistry__sig_and_elbow_pc
```

| arg | default | meaning |
|---|---|---|
| `--path` | – | h5ad with a CSR raw-count layer |
| `--outdir` | – | results, logs, tmp |
| `--layer` | `counts` | CSR raw-count layer (`X` = adata.X) |
| `--id-col` | `participant_id` | donor column (NEBULA random effect) |
| `--covs` | – | NEBULA fixed effects (interactions `a:b`, `a*b` allowed) |
| `--offset-col` | `nCount_RNA` | `nCount_RNA` = per-cell total of `--layer` over **all** genes; anything else = obs column |
| `--contrast-col`, `--case-val` | `case_control`, `pd` | the tested contrast. **Must be in the final covariates** |
| `--run-without-contrast` | 0 | 1 = allow a model without the contrast term (rank purity is then refused) |
| `--n-folds`, `--n-cores` | – | gene chunks (nnz-balanced, each < 2^31 nnz) / parallel NEBULA jobs |
| `--chunk-workers`, `--io-workers` | 8, 24 | parallel chunk writers / h5ad readers |
| `--compute-control-pcs` | 0 | 1 = run step 0 (below) and **rewrite the h5ad's control PCs in place** |
| `--control-val` | `ctr` | control level of `--contrast-col` |
| `--n-hvgs`, `--n-pcs`, `--pc-fdr` | 2000, 30, 0.05 | PCA size; FDR cut for "significant" |
| `--pc-test-covs` | `age,sex,log10_nUMI,frac_mito,frac_intronic,chemistry` | donor-level covariates of the PC case/control test |
| `--pca-model-dir` | `<h5ad dir>/pca_models` | where the PCA model and test table go |
| `--pcs-to-add` | – | one or more of `elbow`, `significant`, `all`; their **union** (index ≤ `--pc-max`) is appended to `--covs`. `all` = every computed PC, i.e. `control_PC_01..--pc-max` |
| `--pc-max` | 20 | highest PC index that may be added |
| `--get-rank-purity` | 0 | 1 = add a `rank_purity` column |
| `--purity-scale` | `log-rate` | `log-rate` (NEBULA's scale, recommended) or `mean` of `--purity-layer` |
| `--purity-layer` | `norm_log` | only used with `--purity-scale mean` |
| `--derive-region-dummies` | 0 | rebuild `region_SN`/`region_VTA` from `region` (contains rule: SN_VTA → 1,1) |
| `--save-tmp`, `--suffix` | 0, – | keep the tmp dir / output-name suffix |

**Fail-fast rules:**
- The contrast must be in the covariates, unless `--run-without-contrast 1` is passed. Without
  this rule, NEBULA would silently fit no contrast coefficient.
- `--pcs-to-add` together with `control_PC*` already in `--covs` is an error.
- `--pcs-to-add` on an h5ad whose `uns/control_pca` is old-format (unpadded `control_PC1`, no
  `created` stamp) is an error; run once with `--compute-control-pcs 1`.
- Any covariate with missing values stops chunking.
- Re-using chunks built from different parameters, **including PCs recomputed since**
  (`uns/control_pca/created` is part of the chunk parameters), is an error.

## Steps
0. **Control PCs** (`python/control_pcs.py compute`, if `--compute-control-pcs 1`). Counts are
   streamed in parallel CSR row blocks. Nothing larger than one block × n_hvgs is ever dense.
   1. Pass 1 (control cells): each row is normalised to 1e4 by its total over **all** genes. Per
      gene, it accumulates detection, Σx and Σx². Genes detected in > 1% of controls are
      eligible. HVGs = top 2000 by scanpy `flavor="seurat"` normalised dispersion, computed
      exactly from those moments.
   2. Pass 2 (control cells, HVGs only, sparse): accumulates Σ, Σ² and XᵀX of log1p values.
      - Scaling: mean/sd (ddof=0, sd + 1e-8). The top 30 eigenvectors of the scaled Gram matrix
        are the loadings, which is exactly PCA/TruncatedSVD of the scaled dense matrix.
      - Sign: each loading vector's largest |entry| is positive.
      - The elbow (`pick_pcs`) is taken over the 30 PCs.
      - The model is saved to `<slogan>_control_pca.npz` + `.json`.
   3. Pass 3: projects **all** cells with the saved model → `control_PC_01 … control_PC_30`.
   4. Donor-level case/control test per PC (see "Choosing PCs").
      - Output: `<slogan>_control_pc_test.csv` (beta, se, t, df, p, fdr, n_donors).
   5. In-place h5ad update (details below).
1. **PC selection** (`control_pcs.py select`, if `--pcs-to-add`): the union of the requested sets
   (elbow / significant / all), index ≤ `--pc-max`, read from `uns/control_pca`. The final covariates are logged and saved to
   `logs/run_args_*.json` (`covs_final`, `pcs_added`, `control_pca_created`).
2. **Chunking** (`python/h5ad_to_gene_chunks.py` + `R/write_chunk_qs.R`): gene-range chunks are
   written straight from the CSR layer as Seurat `.qs` files (all cells × a gene range). The
   offset `nUMI_nebula` is the full-gene `nCount_RNA`.
3. **NEBULA** (`R/run_nebula.R`, one GNU-parallel job per chunk, LN method).
4. **Combine** (`R/combine_nebula_results.R`): summary + padj + −log10 columns.
5. **Rank purity** (`python/pseudobulk_h5ad.py` → `python/rank_purity.py` →
   `R/add_rank_purity_to_combined.R`). The pseudobulk runs in the background during step 3.
6. The tmp dir is removed unless `--save-tmp 1`.

### In-place h5ad update and backup
Only the `obs` `control_PC*` datasets and `uns/control_pca` change; the obs group is not rewritten.
- **Backup (always):** the current `obs` group and `uns/control_pca` are copied with h5py into
  `<h5ad dir>/backup_<timestamp>/<slogan>_obs_uns.h5`. This is a few MB to a few hundred MB, not a
  copy of the whole h5ad. `<h5ad dir>/metadata/<slogan>_obs.csv`, if present, is backed up and
  then rewritten with the same column change.
- **Swap:**
  - New PCs and the new `uns/control_pca` are written to temporary root groups and read back.
  - Then every existing `obs/control_PC*` (old unpadded or padded) is deleted, the new datasets
    are moved in, and obs `column-order` is updated.
  - `uns/control_pca` is replaced the same way.
- **Verification:**
  - The on-disk addresses of all other obs columns, the layers/X shapes and nnz, and the other
    uns keys must be unchanged.
  - Obs is re-read with anndata.
  - A leftover temporary group from an aborted run is an error.
- **Restore** from a backup:
  ```python
  import h5py
  with h5py.File(h5ad, "r+") as f, h5py.File(backup, "r") as b:
      del f["obs"]; f.copy(b["obs"], f, "obs")
      del f["uns/control_pca"]; f.copy(b["uns_control_pca"], f["uns"], "control_pca")
  ```
- HDF5 does not shrink files, so the file grows by the size of the deleted PC datasets.
- New `uns/control_pca` keys: `pcs`, `variance_explained`, `n_elbow_pcs`, `elbow_pcs`,
  `significant_pcs`, `fdr_threshold`, `test_method` (`pseudobulk_ols`), `model`,
  `covariates_dropped`, `case_control_test{pc,beta,se,t,p,fdr}`, `n_hvgs`, `contrast_col`,
  `control_val`, `case_val`, `model_path`, `created`.

### Difference from the old `clean_h5ad/compute_pcs.py`
The old code normalised control cells by their total over > 1%-expressed genes when fitting, but
by the full total when projecting. Here both use the full total. Dry runs (`--no-write`) on the real
objects:
- **Oligodendrocyte** (1.2M cells, 4.26B nnz): 46 s, peak 19.5 GB with 24 workers (memory scales
  with `--io-workers`). New vs old PCs: |r| ≥ 0.998 for PC1–22. Same elbow (5) and the same
  variance explained.
- **Immune**: 5 s. |r| 0.96–1.00 for PC1–16. A few low-variance PCs with nearly tied eigenvalues
  (PC22/23, PC30) rotate. Same elbow (8).

## Choosing PCs: donor-level OLS instead of a cell-level mixed model
The previous selection test was `PC ~ case_control + covs + (1|participant_id)`, fitted on every
cell with lmerTest. The new test is OLS on per-donor PC means:
`mean PC ~ case_control + donor means / level fractions of --pc-test-covs`. Donors with < 5 cells
are dropped, as in `run_nebula.R`.

- **Mixed model, pros:** within-donor covariates are adjusted at the cell level, and donors are
  weighted by precision.
- **Mixed model, cons:**
  - It assumes one shared residual variance across donors.
  - When between-donor variance is small, donors with many cells dominate.
  - With ~10⁵–10⁶ cells, small mean shifts on low-variance PCs reach significance. Immune flagged
    15/30 PCs, including PC23–30 at about 1% variance each.
  - It is slow (17–29 s per fit) and can give singular fits.
- **Donor-level OLS, pros:** the donor is the replicate, each donor counts once, it is immune to
  cell-count imbalance, and it runs in about a second.
- **Donor-level OLS, cons:** donors with few cells are noisy (hence the minimum-cell filter), and
  within-donor covariates become donor means/fractions.
- Significant PCs, OLS vs lmer: immune 9 vs 15 (all 9 are in the lmer set); oligodendrocyte 13 vs 19
  (12 of the 13 are in the lmer set).

**Caveat on "significant" PCs in general.** The PCs are learned on controls and the cases are
projected onto them. A PC that differs between PD and control may *be* part of the disease effect,
for example a shift in cell state or composition. Adjusting for it removes real DE signal (adjusting
for a mediator). Selecting covariates by their association with the tested contrast, in the same
data, also biases the DE results. Elbow PCs are chosen without looking at the contrast, so
`--pcs-to-add elbow` is the conservative choice. `--pc-max` keeps low-variance PCs out either way.
`--pcs-to-add all --pc-max K` takes PCs 01..K regardless of elbow or test, which is also chosen
without looking at the contrast and covers cell-type / composition axes beyond the elbow (useful
when PCs are shared across objects, e.g. one pan-neuron PCA applied to per-cell-type subsets).

## Rank purity
`python/rank_purity.py` ports `compute_rank_purity()` (`~/sc-online/add_rank_purity.R`,
groupby = donor) with these fixes:
- **Scale (default `--purity-scale log-rate`):** the donor pseudobulk is `log(r + eps)`.
  - `r` = the donor mean of counts ÷ offset.
  - `eps` = 1% of the gene's mean `r`, the same for every donor.
  - NEBULA's coefficients are on this scale, so subtracting the fitted `B·X` is unit-consistent.
  - On oligodendrocyte, a `log1p`-normalised pseudobulk disagreed in sign with logFC for 57% of
    FDR<0.05 genes, against 1.3% with `log-rate`.
- **Sign check:** uses the coefficient `contrast_col + case_val` (e.g. `case_controlpd`).
  Upstream looked up `case_control` and silently skipped the check.
- **Categoricals that vary within a donor** (chemistry, `dapi_nurr`) enter as per-donor level
  fractions, not the first cell's level.
- **Covariate matching by type:** numeric → exact name; categorical → name + observed level.
  Upstream's prefix match fails on `control_PC_1` vs `control_PC_10`, and every model coefficient
  must be claimed by exactly one covariate.
- **Donor filtering:** donors with < 5 cells or missing values are dropped (logged).

`add_rank_purity_existing.sh` adds rank purity to an existing combined result (e.g. from upstream
`run-nebula-parallel.sh`), writing copies and never modifying the source.

## Outputs
```
<outdir>/
|-- de_results/nebula_ln/
|   |-- <slogan>__nebula_ln_combined[__<suffix>].qs / .csv
|   |-- <slogan>__nebula_ln_combined[__<suffix>]__rank_purity.tsv
|-- logs/run_<slogan>[__<suffix>]_<ts>.log, logs/run_args_<ts>.json
|-- <slogan>__nebula_ln/          # tmp, removed unless --save-tmp 1
<h5ad dir>/pca_models/<slogan>_control_pca.{npz,json}, <slogan>_control_pc_test.csv   (step 0)
<h5ad dir>/backup_<ts>/<slogan>_obs_uns.h5                                            (step 0)
```
Because the h5ad and `pca_models/` live under `/mnt/accessory`, step 0 must be run by a user with
write access there.

## Vendored code (copied, not referenced)
| file | origin | changes |
|---|---|---|
| `R/run_nebula.R` | `~/sc-online/scripts/run-nebula.R` | self-activation; `utils.R` source removed (unused, pulls biomaRt/sva); `nebula::scToNeb` → local `scToNeb_seurat` (same code, `GetAssayData(layer=)`, since `slot=` is defunct in SeuratObject ≥ 5.4) |
| `R/combine_nebula_results.R` | `~/sc-online/scripts/combine-nebula-results.R` | self-activation; Seurat-based rank-purity branch removed (Python here) |
| `python/control_pcs.py::pick_pcs` | `clean_h5ad/compute_pcs.py::_pick_pcs` | none |
| `python/*`, `R/write_chunk_qs.R`, `R/add_rank_purity_to_combined.R` | `nebula_h5ad/` | env paths |

## Tests (`test/`)
- `simulate_h5ad.py OUT [--pc-test]`: the default sim is identical to the `nebula_h5ad` one.
  `--pc-test` adds latent cell-state programs (one PD-shifted), QC covariates, stale old-format
  `control_PC1..30`, an old `uns/control_pca` and an unrelated `uns/keep_me`.
- **Parity:** this workflow (`--purity-scale mean`) reproduces the `nebula_h5ad` / upstream result
  on the default sim. All 49 statistics agree within 1e-14, including rank purity, and the
  combined `.qs` loads in r-4.3.0.
- `check_control_pcs.py` (run with a python that has scanpy + sklearn, e.g. r-4.3.0):
  - Against a scanpy + `TruncatedSVD(arpack)` reference: identical HVG set, loadings and
    projections |cos| / |r| = 1.000000, and identical variance explained.
  - The donor OLS equals R `lm()`.
  - The in-place write leaves all other obs columns, layers and uns entries identical, and the
    backup holds the original obs/uns.
- **Orchestrator checks:**
  - a full run with `--compute-control-pcs 1 --pcs-to-add elbow significant`;
  - the contrast missing from covs → error, and `--run-without-contrast 1` runs (with purity → error);
  - `--pcs-to-add` + `control_PC*` in covs → error;
  - recomputed PCs + stale chunks → error;
  - `select` on old-format uns → error;
  - a leftover temp group → error;
  - rerunning step 0 is idempotent.
