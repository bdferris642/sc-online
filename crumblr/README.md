# crumblr workflow

A Python/R pipeline for **cell-type composition analysis** of single-nucleus RNA-seq data using the [crumblr](https://github.com/GabrielHoffman/crumblr) package. The workflow aggregates cell-level data to sample-level counts, applies a CLR (centered log-ratio) transformation, fits a mixed-effects model via `dream()`, and generates a suite of publication-ready figures.

---

## Requirements

- **Linux x86-64** (micromamba binary is architecture-specific)
- `curl` and `bash` available on `PATH`
- Internet access for the one-time setup

---

## Installation

Clone or copy the `crumblr/` directory to your machine, then run the setup script once:

```bash
bash /path/to/crumblr/setup_env.sh
```

This will:
1. Download `micromamba` into `../bin/` (relative to `crumblr/`)
2. Create a conda environment named `crumblr` with R 4.4, Bioconductor packages (`variancePartition`, `edgeR`, `limma`, `SingleCellExperiment`, `ggtree`, …), and the `crumblr` R package from GitHub
3. Create a Python virtual environment (`../crumblr-venv/`) with `numpy`, `pandas`, `matplotlib`, `seaborn`, and `scipy`

All files are confined to the directory containing `crumblr/` — nothing is installed system-wide.

**Expected layout after setup:**

```
parent-dir/
├── bin/                        # micromamba binary
├── micromamba_root/            # conda environments
│   └── envs/crumblr/          # R 4.4 + Bioconductor
├── crumblr-venv/               # Python virtual environment
└── crumblr/                    # this repo
    ├── README.md
    ├── setup_env.sh
    ├── run_workflow.py
    ├── crumblr_analysis.R
    ├── simulate_data.py
    ├── environment.yml         # pinned conda packages
    └── requirements.txt        # pinned Python packages
```

> **Note:** `environment.yml` and `requirements.txt` pin the exact package versions used during development. The setup script uses them automatically to produce a reproducible environment.

---

## Quick start with simulated data

Generate a small synthetic snRNA-seq dataset (40 donors, 33 cell types, Parkinson's disease vs. control) and run the full analysis:

```bash
# 1. Generate simulated data
python /path/to/crumblr/simulate_data.py \
    --output /path/to/data/simulated_cells.csv

# 2. Run the analysis
python /path/to/crumblr/run_workflow.py \
    --csv            /path/to/data/simulated_cells.csv \
    --proportion-col cell_type \
    --sample-col     sample_id \
    --donor-col      donor_id \
    --contrast       case_control \
    --ref-level      ctr \
    --fixed-effects  age sex region_SN dapi_nurr brain_bank study \
    --random-effects donor_id \
    --output-dir     /path/to/results/
```

No manual environment activation is needed — `run_workflow.py` detects and uses the `crumblr-venv` automatically.

---

## Running with your own data

```bash
python /path/to/crumblr/run_workflow.py \
    --csv            cells.csv \
    --proportion-col cell_type \
    --sample-col     sample_id \
    --donor-col      donor_id \
    --contrast       diagnosis \
    --ref-level      control \
    --fixed-effects  age sex batch \
    --random-effects donor_id \
    --output-dir     results/
```

---

## Input format

The input is a **cell-level CSV** — one row per cell:

| sample_id | donor_id | cell_type     | diagnosis | age | sex    | batch |
|-----------|----------|---------------|-----------|-----|--------|-------|
| S001_SN   | D001     | Da_SOX6_EYA4  | pd        | 72  | Male   | B1    |
| S001_SN   | D001     | Inh_PENK      | pd        | 72  | Male   | B1    |
| S002_SN   | D002     | Da_CALB1_CALCR| ctr       | 65  | Female | B2    |

- Rows are individual cells
- `--sample-col` identifies the observation unit passed to `dream()` (e.g. `sample_id` when donors contribute multiple samples)
- `--donor-col` tracks the subject for random-effects bookkeeping; if each donor is a single sample, `--sample-col` can be omitted and `--donor-col` will be used for both

---

## Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--csv` | yes | — | Path to cell-level CSV |
| `--proportion-col` | yes | — | Column whose values define cell type labels (proportions computed over this) |
| `--contrast` | yes | — | Main variable of interest (categorical or continuous) |
| `--output-dir` | yes | — | Directory for all outputs (created if absent) |
| `--sample-col` | no | `--donor-col` | Sample identifier; rows are aggregated within each unique value |
| `--donor-col` | no | — | Subject/donor identifier (used for random effects) |
| `--fixed-effects` | no | `[]` | Space-separated covariate column names added to the model as fixed effects |
| `--random-effects` | no | `[]` | Column names modelled as random effects `(1\|col)` via `dream()` |
| `--ref-level` | no | — | Reference level for the contrast variable (sets the baseline group) |
| `--n-top` | no | `60` | Number of cell types shown in boxplot/scatter panels (split evenly enriched/depleted) |
| `--min-n-cells` | no | `100` | Drop samples with fewer than this many cells before analysis |
| `--min-n-cells-category` | no | `150` | Drop cell types with fewer than this many total cells before analysis |
| `--color-by` | no | — | Column used to colour heatmap rows / scatter points when `--contrast` is continuous (e.g. `case_control` when contrasting by `age`) |
| `--r-script` | no | `crumblr_analysis.R` | Path to the R analysis script |
| `--r-libs` | no | auto | Path to a custom R library directory |

---

## Outputs

All files are written to `--output-dir`:

| File | Description |
|---|---|
| `crumblr_results.csv` | Full `topTable` output for all cell types and contrast coefficients |
| `crumblr_results_sorted.csv` | Same, sorted by p-value |
| `clr_values.csv` | CLR-transformed proportions (samples × cell types) |
| `precision_weights.csv` | Per-observation precision weights from `crumblr()` |
| `counts_matrix.csv` | Aggregated sample × cell-type count matrix passed to R |
| `metadata.csv` | Sample-level metadata passed to R |
| `run_args.json` | Record of all CLI arguments used |
| `volcano_<coeff>.pdf/png` | Volcano plot: log2FC vs −log10(FDR) per contrast coefficient |
| `forest_<coeff>.pdf` | Forest plot of effect sizes with 95% CIs (generated in R) |
| `logfc_bar_<coeff>.pdf/png` | Horizontal bar chart of log2FC for all cell types |
| `proportion_boxplots.pdf/png` | Box + strip plots for top N cell types *(categorical contrast)* |
| `proportion_scatter.pdf/png` | Scatter with OLS regression line per cell type *(continuous contrast)* |
| `composition_bar.pdf/png` | Stacked bar of mean proportions per group *(categorical contrast)* |
| `cell_class_summary.pdf/png` | Broad cell-class boxplots *(categorical contrast)* |
| `clr_heatmap.pdf/png` | Hierarchically-clustered CLR heatmap (samples × cell types) |

---

## How it works

1. **Aggregate** — cell-level rows are grouped by `--sample-col` and `--proportion-col` to produce a sample × cell-type count matrix
2. **Filter** — samples and cell types below minimum cell thresholds are dropped
3. **CLR transform** — `crumblr()` applies a centered log-ratio transformation and computes precision weights that account for library-size variation
4. **Model fit** — `dream()` (lme4-based mixed-effects) or `lmFit` (fixed-effects only) is fit per cell type; collinear covariates are detected and dropped automatically with a warning
5. **Results** — `topTable()` extracts log2FC, t-statistics, p-values, and BH-adjusted FDR per cell type per contrast coefficient
6. **Plots** — Python-side figures (volcano, boxplots, heatmap, etc.) and an R-side forest plot are written to `--output-dir`
