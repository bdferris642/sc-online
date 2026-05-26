# cnmf-analysis

Pipeline for characterizing cNMF gene expression programs (GEPs) from consensus spectra files. Given a manifest describing one or more cNMF runs, it runs three complementary analyses per run:

1. **Regression** — beta regression of GEP scores against covariates, with collinearity detection, Spearman correlation matrices, and a full plot suite
2. **Intersection** — Fisher's exact test gene-set enrichment using a top-gene threshold, with redundancy filtering and cross-GEP summary heatmaps
3. **GSEA** — ranked GSEA (fgsea) on the full loading distribution per GEP, with redundancy filtering

All three steps are managed by the orchestrator `run_gep_analysis.py`, which can run individual steps or the full pipeline sequentially or in parallel across manifest rows.

---

## Directory layout

```
cnmf-analysis/
  run_gep_analysis.py              ← Orchestrator (entry point)
  get_cnmf_spectra_gsea.py         ← GSEA per spectra file
  get_cnmf_spectra_intersection.py ← Fisher intersection per spectra file
  run_gep_regression.R             ← Regression per manifest row
  run_intersection.R               ← Per-GEP Fisher test (called by intersection wrapper)

../scripts/
  run-gsea.R                       ← Per-GEP GSEA (called by GSEA wrapper)
```

The Python wrappers and R scripts resolve their dependencies (Rscript path, sibling scripts) relative to their own file locations — no PATH manipulation or manual env activation is needed.

---

## Prerequisites

### Conda environment

One environment must exist at `cc-sandbox/micromamba_root/envs/cnmf-analysis-env/`. All R scripts self-activate it — you do not need to activate anything before running.

**Packages included:** R 4.4, fgsea, ggplot2, ggrepel, glue, getopt, dplyr, stringr, scales, igraph, betareg, nnls, svglite, tidyr, SeuratObject, qs, qs2, optparse

To create (if missing):
```bash
bash /path/to/cc-sandbox/cnmf-analysis/setup.sh
```

### External dependencies (absolute paths, not in this repo)

| Resource | Path | Used by |
|----------|------|---------|
| Shared R library (`gsea.R`, `plot.R`) | `/home/ferris/sc-online/` | `run-gsea.R`, `run_intersection.R` |
| Gene set GMT files | `/home/ferris/sc-online/genesets/latest/` | `run-gsea.R`, `run_intersection.R` |
| Ensembl→HGNC symbol map | `/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv` | Python wrappers |

---

## Manifest format

Tab-separated file (`.tsv`). One row per cNMF run.

### Required columns

| Column | Description |
|--------|-------------|
| `name` | Short identifier used as output subdirectory name (e.g. `sox6-agtr1-k12`) |
| `spectra_path` | Absolute path to the GEPs×genes consensus spectra file (e.g. `*.spectra.k_12.dt_0_05.consensus.txt`) |
| `qs_path` | Absolute path to the Seurat object (`.qs` or `.qs2`) |
| `feature_loadings_path` | Absolute path to the cNMF feature loadings file |
| `qc_geps` | Comma-separated 1-based indices of QC-passing GEPs (e.g. `1,3,5,7`) |
| `id_col` | Column in Seurat metadata holding the participant/sample identifier |
| `covs` | Comma-separated covariate column names for regression (e.g. `age,sex,batch`) |

### Optional columns (empty or `NA` → use defaults)

| Column | Default | Description |
|--------|---------|-------------|
| `ncores` | `detectCores()-1` | Parallel cores for regression |
| `intersection_top_n` | `150` | Max genes selected per GEP after z-score filter |
| `intersection_zscore` | `3.0` | Z-score threshold: select genes with loading > mean + k·sd |
| `intersection_min_size` | `15` | Min gene-set size after background filtering |
| `intersection_max_size` | `250` | Max gene-set size after background filtering |
| `intersection_jobs` | n_geps | Parallel R jobs for intersection step |
| `gsea_jobs` | n_geps | Parallel R jobs for GSEA step |

### Example manifest

```tsv
name	spectra_path	qs_path	feature_loadings_path	qc_geps	id_col	covs	ncores	intersection_top_n
sox6-agtr1-k12	/mnt/.../my-run.spectra.k_12.dt_0_05.consensus.txt	/mnt/.../sox6.qs	/mnt/.../loadings.txt	1,3,5,7,9,11	donor_id	age,sex,pmi	48	200
sox6-k15	/mnt/.../my-run.spectra.k_15.dt_0_05.consensus.txt	/mnt/.../sox6.qs	/mnt/.../loadings.txt	2,4,6,8,10		age,sex,pmi
```

---

## Usage

### Full pipeline (all three steps)

```bash
python /home/ferris/cc-sandbox/cnmf-analysis/run_gep_analysis.py \
  --manifest /mnt/accessory/analysis/pd-freeze/sox6-gep-manifest.tsv \
  --outdir   /mnt/accessory/analysis/pd-freeze/gep-analysis
```

### Specific steps only

```bash
python /home/ferris/cc-sandbox/cnmf-analysis/run_gep_analysis.py \
  --manifest manifest.tsv \
  --outdir   /path/to/output \
  --steps    intersection,gsea
```

### Parallel rows

```bash
python /home/ferris/cc-sandbox/cnmf-analysis/run_gep_analysis.py \
  --manifest manifest.tsv \
  --outdir   /path/to/output \
  --jobs     4
```

Note: each step within a row is always sequential (regression → intersection → GSEA). `--jobs` parallelises across manifest rows only.

---

## Output layout

```
<outdir>/
  <name>/
    gep_regression/
      1_<qs_basename>/           ← created by run_gep_regression.R
        gep_regression_results.csv
        spearman_corr_*.csv
        spearman_pval_*.csv
        volcano_*.svg
        coef_heatmap_*.svg
        r2_barplot.svg
        gep_distribution_*.svg
        pca_biplot_*.svg
        sample_heatmap_*.svg

    intersection/
      <basename>-transpose.txt                  ← cached GEP×gene transpose (shared with gsea)
      <basename>-intersection-input.tsv
      intersection/
        <GEP_01>/
          <basename>_GEP_01_kegg_2026_human_intersection.csv
          <basename>_GEP_01_kegg_2026_human_intersection_redundancy_filtered.csv
          figures/svg/
            <basename>_GEP_01_kegg_2026_human_intersection.svg
          figures/svg/redundancy_filtered/
            ...
        cross_gep/
          kegg_2026_human_cross_gep_neg_log10_padj.svg
          kegg_2026_human_cross_gep_neg_log10_padj.csv
          kegg_2026_human_cross_gep_log2_OR.svg
          kegg_2026_human_cross_gep_log2_OR.csv
          ... (one pair per collection)

    gsea/
      <basename>-transpose.txt                  ← same cache, shared if already present
      <basename>-gsea-input.csv
      gsea/signed/
        <basename>-gsea-input_GEP_01_gsea.csv
        <basename>-gsea-input_GEP_01_gsea_redundancy_filtered.csv
        figures/svg/
          <basename>-gsea-input_GEP_01_kegg_2026_human_gsea.svg
        figures/svg/redundancy_filtered/
          ...
```

---

## Gene set collections

Both intersection and GSEA use the same eight collections:

| Key | File |
|-----|------|
| `kegg_2026_human` | `KEGG_2026.gmt` |
| `reactome_2024` | `Reactome_Pathways_2024.gmt` |
| `gwas_2025` | `GWAS_Catalog_2025.gmt` |
| `jensen_disease_2025` | `Jensen_DISEASES_Curated_2025.gmt` |
| `trrust_2019` | `TRRUST_Transcription_Factors_2019.gmt` |
| `jaspar_tf_2025` | `JASPAR_PWM_Human_2025.gmt` |
| `go_bp_2025` | `GO_Biological_Process_2025.gmt` |
| `go_mf_2025` | `GO_Molecular_Function_2025.gmt` |

---

## Running individual components

Each wrapper can be called directly, bypassing the orchestrator.

### GSEA only

```bash
python /home/ferris/cc-sandbox/cnmf-analysis/get_cnmf_spectra_gsea.py \
  --path   /path/to/my-run.spectra.k_12.dt_0_05.consensus.txt \
  --outdir /path/to/output \
  --jobs   12
```

### Intersection only

```bash
python /home/ferris/cc-sandbox/cnmf-analysis/get_cnmf_spectra_intersection.py \
  --path     /path/to/spectra.txt \
  --outdir   /path/to/output \
  --top-n    200 \
  --zscore   3.0 \
  --jobs     12
```

### Regression only

```bash
# run_gep_regression.R self-activates cnmf-analysis-env; any Rscript will do
Rscript --no-save --no-restore --no-site-file --no-init-file \
  /home/ferris/cc-sandbox/cnmf-analysis/run_gep_regression.R \
  --input      /path/to/manifest.tsv \
  --output-dir /path/to/output \
  --ncores     48
```

---

## Caching behaviour

The Python wrappers write cache files on first run:

| File | Contents | Shared between |
|------|----------|----------------|
| `<basename>-transpose.txt` | Transposed spectra (genes×GEPs) + HGNC symbols | intersection + GSEA |
| `<basename>-gsea-input.csv` | gene + GEP score columns (no gene_name) | GSEA R script |
| `<basename>-intersection-input.tsv` | gene_name + gene_id + GEP columns | Intersection R script |

If the cache files already exist they are reused. Delete them to force a rebuild (e.g. after the gene symbol map is updated).

---

## Portability notes

**This pipeline is not fully portable as-is.** The following are hardcoded to the current machine and must be updated when moving to a new environment:

### 1. Gene symbol map — Python wrappers

Both `get_cnmf_spectra_gsea.py` and `get_cnmf_spectra_intersection.py` contain:
```python
ENSG_TO_SYMBOL = "/mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/ensg_to_symbol.csv"
```
Update this constant at the top of each file.

### 2. Shared R library paths — R scripts

`run_intersection.R` and `scripts/run-gsea.R` both contain:
```r
SCONLINE_DIR = "/home/ferris/sc-online"
```
This directory holds `gsea.R`, `plot.R`, and `genesets/latest/*.gmt`. Either copy it alongside the pipeline or update `SCONLINE_DIR` in both R files.

### 3. Conda environment — must be recreated

The `cnmf-analysis-env` environment has absolute paths baked into its conda metadata and cannot be moved; it must be recreated at the new `micromamba_root` location. The R self-activation blocks find it relative to the script's location, so the directory tree just needs to be structurally preserved.

### 4. What is portable without changes

- All scripts in `cnmf-analysis/` (paths computed relative to `__file__`)
- The self-activation logic in both R scripts (relative to script location)
- Any output already written (CSVs, SVGs)

### Summary of changes needed for a new machine

| Item | File(s) to edit | Change |
|------|----------------|--------|
| Gene symbol map path | `get_cnmf_spectra_gsea.py`, `get_cnmf_spectra_intersection.py` | Update `ENSG_TO_SYMBOL` |
| Shared R library | `run_intersection.R`, `scripts/run-gsea.R` | Update `SCONLINE_DIR` |
| Conda env | — | Recreate with micromamba (see Prerequisites above) |
| micromamba binary | — | Download to `<sandbox>/bin/micromamba` |
