# pairwise-dist-hist

Computes and visualises pairwise distance distributions between cell types in single-cell data. Given an AnnData `.h5ad` file and an embedding (e.g. `X_pca`), it samples pairs of cells within and across each cell-type combination, computes distances, and produces KDE density histograms — useful for assessing how well-separated cell populations are in embedding space.

For cross-type pairs, the joint plot annotates each panel with **Hartigan's dip statistic** (and Sarle's bimodality coefficient) to quantify whether the combined within/cross distance distribution is unimodal or bimodal. The matrix figures use dip test p-values with significance stars only.

---

## Outputs

Running the main script produces five sets of files:

| File | Contents |
|------|----------|
| `<output>.pdf/.png/.svg` | Grid of all cell-type pairs (including same-type within). 3 columns, KDE curves for within-i, cross, within-j. |
| `<output>_joint.pdf/.png/.svg` | Cross-type pairs only. Combined distribution (gray) overlaid with three component KDEs. Bimodality statistics annotated per panel. |
| `<output>_matrix_joint.pdf/.png/.svg` | Upper-triangle N×N matrix. Diagonal: within-type KDE. Upper triangle: single filled KDE of pooled joint distribution (within_i + 2×cross + within_j) in blended color, annotated with Hartigan's dip statistic and p-value (significance stars). Shared X axis, independent Y. |
| `<output>_matrix_between.pdf/.png/.svg` | Upper-triangle N×N matrix. Diagonal: within-type KDE. Upper triangle: cross-type-only KDE in blended color. Shared X axis, independent Y. |
| `<output>_matrix_data.csv` | Tidy distance data for both matrix figures. Columns: `ct_i`, `ct_j`, `pair_type`, `distance`. |

---

## Setup

Dependencies are managed in a Python virtual environment. Run the setup script once from any location:

```bash
bash /path/to/pairwise-dist-hist/setup_env.sh
```

This creates an `adata-venv` environment in the **parent directory** of `pairwise-dist-hist/` and installs all required packages. The analysis scripts activate this environment automatically when run directly with `python3`.

**Requirements:**
- Python 3.6+ with the `venv` module
- `curl` (used to bootstrap `pip`)
- Internet access (first-time setup only)

---

## Usage

### Generate test data (optional)

```bash
python3 simulate_data.py --out test/my_test.h5ad --n-cells 3000 --seed 0
```

This creates a synthetic `.h5ad` with 5 cell types (TypeA/B/C/D/E), 50 features, and a 20-dimensional PCA embedding in `obsm['X_pca']`. TypeE is drawn from the same distribution as TypeA (same centre, same noise), making A↔E cross-type distances indistinguishable from within-A — a unimodal negative control.

### Run the analysis

```bash
python3 pairwise_dist_hist.py \
    --h5ad     <path/to/data.h5ad> \
    --obs-col  <cell_type_column> \
    --obsm-key <embedding_key> \
    --output   <output_prefix.pdf> \
    [--n-pairs  10000] \
    [--metric   cosine] \
    [--seed     42]
```

**Arguments:**

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--h5ad` | yes | — | Path to AnnData `.h5ad` file |
| `--obs-col` | yes | — | `obs` column containing cell-type labels |
| `--obsm-key` | yes | — | `obsm` key for the embedding (e.g. `X_pca`) |
| `--output` | yes | — | Output PDF path (stem is reused for all output files) |
| `--n-pairs` | no | `10000` | Number of pairs sampled per within/cross category |
| `--metric` | no | `cosine` | Distance metric: `euclidean`, `cosine`, or `manhattan` |
| `--seed` | no | `42` | Random seed for reproducibility |

**Example using test data:**

```bash
python3 pairwise_dist_hist.py \
    --h5ad     test/test_data.h5ad \
    --obs-col  cell_type \
    --obsm-key X_pca \
    --output   test/my_output.pdf
```

Expected outputs (10 files, cosine metric by default):
```
test/my_output.pdf/.png/.svg
test/my_output_joint.pdf/.png/.svg
test/my_output_matrix_joint.pdf/.png/.svg
test/my_output_matrix_between.pdf/.png/.svg
test/my_output_matrix_data.csv
```

---

## Directory structure

```
pairwise-dist-hist/
├── pairwise_dist_hist.py   # Main analysis script
├── simulate_data.py        # Synthetic data generator for testing
├── setup_env.sh            # One-time environment setup
└── test/
    ├── test_data.h5ad      # Pre-generated test dataset (3000 cells; TypeE drawn from same distribution as TypeA)
    ├── test_output.pdf     # Example grid output (cosine)
    ├── test_output_joint.pdf
    ├── test_output_matrix_joint.pdf
    └── test_output_matrix_between.pdf
```

The virtual environment is created **one level up** from this folder (`../adata-venv/`) and is excluded from version control.

---

## Dependencies

| Package | Purpose | Version |
|---------|---------|---------|
| `anndata` | Read/write `.h5ad` files | ≥0.10 |
| `numpy` | Array ops, sampling, distance computation | ≥1.24 |
| `pandas` | Label handling | ≥2.0 |
| `scipy` | Skewness/kurtosis statistics | ≥1.10 |
| `scikit-learn` | PCA (simulation only) | ≥1.3 |
| `matplotlib` | Figure/PDF rendering | ≥3.7 |
| `seaborn` | KDE plots, colour palettes | ≥0.13 |
| `h5py` | HDF5 backend for AnnData | ≥3.9 |
| `diptest` | Hartigan's dip test for multimodality | ≥0.10 |

---

## Methods notes

- **Balanced sampling:** When two cell types differ in size, the larger is subsampled to match the smaller before cross-type pairing, ensuring equal representation.
- **Dip test:** Hartigan's dip statistic tests departure from unimodality; small p-values suggest multimodality. Significance is indicated by stars: `*` p<0.05, `**` p<0.01, `***` p<0.001. The `_joint` flat-list plot additionally reports Sarle's bimodality coefficient (BC > 5/9 ≈ 0.556 suggests bimodality); the matrix figures show dip p-values only.
- **Pooled joint distribution:** For `matrix_joint`, each upper-triangle panel combines within_i + 2×cross + within_j to reflect the 1:2:1 frequency ratio expected when drawing random pairs from the merged cell pool.
- **Sampling:** Distances are estimated from `n-pairs` randomly sampled pairs per category rather than exhaustive pairwise computation, making the tool practical for large datasets.
- **Default metric:** Cosine distance is the default, appropriate for normalized single-cell embeddings.
