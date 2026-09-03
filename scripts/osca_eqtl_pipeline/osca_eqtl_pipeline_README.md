# OSCA eQTL Pipeline

Single-nucleus RNA-seq pseudobulk eQTL analysis pipeline: pseudobulk aggregation → OSCA linear-model association → mashr effect-size sharing.

---

## Quick Start

```bash
./00-run-osca-eqtl-pipeline.sh \
  --ct-id cell_class \
  --gene-log-expr-threshold 0.01 \
  --pb-output-dir /mnt/accessory/seq_data/pd_all/240514/dapi_nurr_merged_seurat_clean_subsets/h5ad_pseudobulk \
  --h5ad-cat-covar "sex case_control study brain_bank dapi_nurr" \
  --h5ad-quant-covar "age pmi" \
  --input-files /mnt/accessory/seq_data/pd-freeze/sn-vta/subsets/latest/latest-h5ad-paths.txt \
  --mash-eps 1e-6 \
  --mash-num-random 100000 \
  --mash-padj-thresh 0.05 \
  --min-num-cells 10 \
  --osca-input-dir /mnt/accessory/seq_data/gtex \
  --participants /mnt/accessory/seq_data/gtex/gtex_participants__10_cells_per_class.txt \
  --pipeline-slogan gtex_osca_0.01 \
  --sample-id participant_id \
  --start-at-step 3 \
  --stop-after-step 6 \
  --vcf-slogan gtex_merged_v8_v9_normed_dbsnp_annot_filtered_dedupe
```

Note: `--h5ad-cat-covar` and `--h5ad-quant-covar` match the current defaults; shown explicitly for clarity. `example-command.txt` in the original pipeline directory is outdated (predates these arguments).

---

## Overview

**eQTLs** (expression quantitative trait loci) are genetic variants that explain variation in gene expression. This pipeline identifies cis-eQTLs (variants within ~1 kb of a gene) in each snRNA-seq cell class independently, then applies **mashr** (multivariate adaptive shrinkage) to estimate the degree to which effect sizes are shared or specific across cell classes.

**OSCA** (OMIC data analysis using a Summary-statistics based Causal inference Approach) is used for the per-cell-class association testing. It encodes expression phenotypes in a binary BOD format and fits a linear model with categorical and quantitative covariates.

**When to use this pipeline:**
- You have snRNA-seq data with matched whole-genome sequencing (WGS/genotypes).
- You want cell-class-resolved cis-eQTL maps.
- You want to compare effect-size sharing across cell classes (mashr).
- You have ≥ ~50 participants with both expression and genotype data.

---

## Prerequisites

### System Tools
| Tool | Purpose |
|------|---------|
| `osca` | Core eQTL association engine |
| `plink` (v1.9) | Generates PLINK bfiles from VCF |
| `bcftools` | VCF preprocessing (not called directly by pipeline, needed for input prep) |
| `GNU parallel` | Parallelizes step 3 across cell classes |
| `gcloud` CLI | Steps 4 and 8: copy outputs to Google Cloud Storage |

### Conda Environment
The orchestrator calls `conda activate mashr` at startup. All R and Python packages below must be installed in this environment.

### R Packages
`getopt`, `Matrix`, `sva`, `tidyverse`, `dplyr`, `ggplot2`, `ggrepel`, `mashr`, `ashr`, `tidyr`

### Python Packages
`decoupler`, `scanpy`, `sccoda`, `pandas`, `numpy`

---

## Input Files

| Input | Format | Notes |
|-------|--------|-------|
| h5ad expression files (`--input-files`) | AnnData `.h5ad` | One file per cell type; must contain a `counts` layer with raw integer counts |
| h5ad file list (`--input-files`) | Plain text | Newline-delimited absolute paths to the h5ad files to process |
| Participants list (`--participants`) | Plain text | One participant ID per line; must match `--sample-id` column in `.obs` |
| ID map CSV (`--id-map`) | CSV | Optional. Columns: `nucseq_participant_id`, `vcf_sample_id`. Remaps h5ad participant IDs to VCF IDs and restricts to matched participants. |
| PLINK bfiles (`--vcf-slogan`) | `.bed/.bim/.fam` | Placed in `--osca-input-dir`; VCF-derived genotype data |
| Genotype PCs | `{vcf-slogan}_pca.eigenvec` | Placed in `--osca-input-dir`; first 5 PCs used as covariates |
| Gene annotation (`--gene-anot`) | Tab-delimited | Columns: `NAME`, `probe`, `chr`, `TSS`, `strand`; default `/mnt/accessory/analysis/eqtl/gene_loc.txt` |
| Metadata (`--metadata`) | CSV | Must contain: `participant_id`, `age`, `sex`, `pmi`; default path is project-specific (see Known Issues) |

---

## Pipeline Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--ct-id` | `cell_class` | Column in `.obs` holding cell type labels |
| `--pb-output-dir` | *(required)* | Directory for pseudobulk output CSVs (written by step 1, read by step 2) |
| `--input-files` | *(required)* | Path to newline-delimited text file listing h5ad files to process |
| `--gene-log-expr-threshold` | `0.01` | Min mean log-expression to retain a gene |
| `--h5ad-cat-covar` | `"sex case_control study brain_bank dapi_nurr"` | Space-delimited list of categorical covariate column names from `.obs`. Passed to step 2 as `cov1` (converted to factors). Also used to auto-build the SVA formula unless `--sva-formula` is provided. |
| `--h5ad-quant-covar` | `"age pmi"` | Space-delimited list of quantitative covariate column names from `.obs`. Passed to step 2 as `cov2` (min-max scaled). Also used to auto-build the SVA formula unless `--sva-formula` is provided. |
| `--id-map` | *(none)* | Optional path to participant ID mapping CSV (see above) |
| `--mash-eps` | `1e-6` | Floor for SE and beta values passed to mashr |
| `--mash-num-random` | `1000000` | Number of random SNP-probe pairs used to estimate null correlation in mashr |
| `--mash-padj-thresh` | `0.01` | BH-adjusted p-value threshold for "strong" tests in mashr |
| `--min-num-cells` | `10` | Minimum cells per sample for pseudobulk inclusion |
| `--osca-input-dir` | *(required)* | Directory for OSCA inputs/outputs; must contain PLINK bfiles |
| `--participants` | *(required unless `--id-map` given)* | Path to newline-delimited VCF participant ID list. If omitted, IDs are derived from the `vcf_sample_id` column of `--id-map`. |
| `--pipeline-slogan` | *(required)* | Short label for this run; used in output paths and GCS bucket |
| `--sample-id` | `participant_id` | Column in `.obs` holding sample/participant IDs |
| `--start-at-step` | `1` | Resume from this step (inclusive) |
| `--stop-after-step` | `1000` | Stop after this step (inclusive) |
| `--strs-to-skip` | *(none)* | Comma-separated strings; `.h5ad` files containing any of these are skipped |
| `--vcf-slogan` | *(required)* | Prefix of `.bed/.bim/.fam` and `_pca.eigenvec` files in `--osca-input-dir` |

---

## Pipeline Steps and Data Flow

```
h5ad files (one per cell type)  +  cell_class list
        │
        ▼  Step 1: 01-make-eqtl-pseudobulk.py
        │  ID remap + restrict (if --id-map) → Pseudobulk sum aggregation → CPM normalize → log1p → scale
        │
   {pb-output-dir}/
     {cell_class}_expression_matrix_ds.csv   (samples × genes, log-normalized, scaled; cell_class from adata.obs[CT_ID])
     {cell_class}_composition_matrix_ds.csv  (samples × cell types, raw counts)
        │
        ▼  Step 2: 02-run-osca-formatting-scanpy.R
        │  SVA → PCA (30 PCs) → merge genotype PCs + SVs + metadata
        │  → remove SVs correlated with covariates (Spearman r > 0.9)
        │  → write OSCA-format files
        │
   {osca_input_dir}/
     Phenotype_{cc}_osca.txt      (FID/IID × genes, wide format)
     Upprobe_{cc}.opi             (probe annotation: chr, NAME, TSS, probe, strand)
     cov1_{cc}.txt                (categorical covariates: sex)
     cov2_{cc}_reduced.txt        (quantitative covariates: age, pmi, PCs 1-20, G-PCs 1-5, SVs)
        │
        ▼  Step 3: 03-build-eqtl.sh  [parallelized via GNU parallel, one job per cell class]
        │  osca --make-bod → --update-opi → --eqtl --cis → --query 1
        │
   {osca_output_dir}/
     eqtl_{cc}.tsv                (large; all cis SNP-probe pairs)
     befile_{cc}.*                (OSCA binary expression files)
     eqtl_{cc}_osca_*.log         (per-step OSCA logs)
        │
        ▼  Step 4 (optional): gcloud storage cp
        │  Backup raw .tsv files to GCS before deletion
        │
        ▼  Step 5: 05-process-and-plot-osca-tsv.R  [parallelized per cell class]
        │  Two-stage FDR: within-gene Bonferroni → BH across genes (padj_gene)
        │              + BH within eGenes (padj_snp)
        │  Saves RDS; generates Manhattan + histogram plots
        │  Deletes .tsv files
        │
   {osca_output_dir}/
     eqtl_{cc}.rds                (all SNP-probe pairs + FDR columns)
     eqtl_{cc}_sig.rds            (significant rows: padj_gene < 0.05 & padj_snp < 0.05)
     plots/eqtl_{cc}_manhattan.png
     plots/eqtl_{cc}_min_p_gene_hist.png
     plots/eqtl_{cc}_padj_gene_hist.png
     plots/eqtl_{cc}_pval_hist.png
     plots/eqtl_{cc}_padj_snp_hist.png
        │
        ▼  Step 6: 06-get-common-snp-probes-osca-rds.R
        │  Intersect SNP-probe pairs across ALL cell classes
        │  Subset to pairs significant in at least one class
        │
   {osca_output_dir}/
     eqtl_present_in_all.rds          (all common SNP-probe pairs, long format)
     eqtl_present_in_all_sig_in_one.rds  (pairs significant [padj_snp < 0.01] in ≥1 class)
        │
        ▼  Step 7: 07-run-mashr.R
        │  mashr: PCA + ED covariances + canonical + custom (neurons/glia)
        │  Fit on random subset → apply to strong tests
        │
   {osca_output_dir}/
     eqtl_present_in_all__mash_results_sig.rds    (mashr posteriors, strong tests)
     eqtl_present_in_all__mash_results_random.rds (mashr fit object)
     eqtl_present_in_all__mash_data_sig.rds       (mash_set_data object, strong)
     eqtl_present_in_all__mash_data_random.rds    (mash_set_data object, random)
        │
        ▼  Step 8 (optional): gcloud storage cp -r
           Archive entire osca_input_dir (inputs + outputs) to GCS
```

Output root: `{osca_input_dir}/eqtl_final_outs/{pipeline_slogan}/`

---

## Per-Script Reference

### Step 1: `01-make-eqtl-pseudobulk.py`

**Purpose:** Aggregate raw counts per sample × cell type, normalize, and filter genes for each `.h5ad` file in the input directory.

**Usage** (called by orchestrator):
```bash
python 01-make-eqtl-pseudobulk.py \
  --input-files /path/to/h5ad_files.txt \
  --output-dir DIR/pseudobulk \
  --min-num-cells 10 \
  --gene-log-expr-threshold 0.01 \
  --ct-id cell_class \
  --sample-id participant_id \
  [--strs-to-skip opc] \
  [--id-map /path/to/id_map.csv]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--input-files` | *(required)* | Path to newline-delimited text file listing h5ad files to process |
| `--output-dir` | `pseudobulk` | Output directory |
| `--min-num-cells` | `10` | Min cells per sample |
| `--gene-log-expr-threshold` | `1.5` | Min mean log-expression to retain gene |
| `--ct-id` | `cell_class` | Cell type column in `.obs` |
| `--sample-id` | `participant_id` | Sample column in `.obs` |
| `--strs-to-skip` | *(none)* | Comma-separated substrings to skip |
| `--id-map` | *(none)* | CSV with `nucseq_participant_id` and `vcf_sample_id` columns; remaps and restricts participants |

**Processing steps per h5ad:**
1. Loads `counts` layer as `X` (raises an error if layer is absent)
2. Validates `--sample-id` and `--ct-id` columns exist in `.obs` (raises an error otherwise)
3. Removes genes/cells with zero counts; removes NaN sample/cell-type labels
4. Casts `--sample-id` column to `str` and strips whitespace
5. If `--id-map` provided: restricts to cells whose participant ID appears in `nucseq_participant_id`; remaps IDs to `vcf_sample_id`; raises an error if no cells remain
6. Validates that each h5ad contains exactly one cell type in `--ct-id` column (raises an error otherwise)
7. Pseudobulk sum aggregation with `decoupler.get_pseudobulk` (min cells = `--min-num-cells`, min counts = 100)
8. CPM normalization (`target_sum=1e6`) + log1p
9. Gene filter: mean log-expression ≥ `--gene-log-expr-threshold`
10. Scale to max value 10

**Outputs per input file — named by `adata.obs[--ct-id].unique()[0]`, not the h5ad filename stem:**
- `{output-dir}/{cell_class}_expression_matrix_ds.csv` — pseudobulk expression (samples × genes)
- `{output-dir}/{cell_class}_composition_matrix_ds.csv` — cell composition (samples × cell types)

**Requirement:** Each input h5ad must contain exactly one cell type (one unique value in `--ct-id` column). Mixed-cell-type files will cause an error.

**Note:** `sc.settings.n_jobs = 30` is hardcoded.

---

### Step 2: `02-run-osca-formatting-scanpy.R`

**Purpose:** Format pseudobulk CSVs into OSCA input files, running SVA and PCA to build covariate matrices.

**Usage** (called by orchestrator):
```bash
Rscript 02-run-osca-formatting-scanpy.R \
  --expression-dir DIR/pseudobulk \
  --output-dir DIR/osca_inputs \
  --vcf-slogan vcf_prefix \
  --participants /path/to/participants.txt \
  [--gene-anot /path/to/gene_loc.txt] \
  [--metadata /path/to/metadata.csv]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--expression-dir` | *(required)* | Directory with pseudobulk CSVs |
| `--output-dir` | *(required)* | OSCA input/output directory (also contains bfiles) |
| `--vcf-slogan` | *(required)* | Prefix of `.fam`, `_pca.eigenvec` files |
| `--participants` | *(required)* | Path to newline-delimited VCF participant ID list. Passed by orchestrator; auto-derived from `--id-map` if not supplied directly. |
| `--h5ad-cat-covars` | *(required unless `--sva-formula` given)* | Space-delimited categorical covariate names (go to `cov1` as factors; used to build SVA formula) |
| `--h5ad-quant-covars` | *(required unless `--sva-formula` given)* | Space-delimited quantitative covariate names (go to `cov2` min-max scaled; used to build SVA formula) |
| `--sva-formula` | *(auto-built from covariate lists)* | Override SVA formula string (e.g. `"~ age + sex + pmi"`); does not affect `cov1`/`cov2` contents |
| `--gene-anot` | `/mnt/accessory/analysis/eqtl/gene_loc.txt` | Gene annotation table |
| `--metadata` | project-specific default (see Known Issues) | CSV with `participant_id` and all covariate columns |

**Covariate flow — what goes where:**

| Source | File | OSCA flag |
|--------|------|-----------|
| `--h5ad-cat-covars` | `cov1_{cc}.txt` (as factors) | `--covar` |
| Expression PCs 1–20, genotype PCs 1–5, SVs 1–5 | `cov2_{cc}_reduced.txt` | `--qcovar` |
| `--h5ad-quant-covars` (min-max scaled) | `cov2_{cc}_reduced.txt` | `--qcovar` |
| Cell composition counts | `cov2_{cc}_reduced.txt` | `--qcovar` |

SVA model is `mod = model.matrix(~ cat_covar1 + ... + quant_covar1 + ..., data = pheno)`, built automatically from the two covariate lists. Use `--sva-formula` only if you need to specify a formula that differs from the covariate lists (e.g. to include interaction terms).

**Processing per cell class:**
1. Loads expression matrix; subsets to participants in `--participants` file
2. Removes low-variance genes (SD < 1e-6)
3. Builds SVA formula from `--h5ad-cat-covars` + `--h5ad-quant-covars`; runs SVA (5 SVs)
4. Runs PCA (30 PCs) on WGS-intersected expression data
5. Merges genotype PCs, SVs, scaled quant covars, cat covars, cell composition
6. Removes SVs with Spearman |r| > 0.9 to any other covariate; drops PCs 21–30
7. Annotates genes via `gene_annotation` table; drops unannotated genes

**Outputs per cell class `{cc}`:**
- `Phenotype_{cc}_osca.txt` — FID/IID + gene expression (OSCA phenotype format)
- `Upprobe_{cc}.opi` — probe annotation (chr, NAME, TSS, probe, strand)
- `cov1_{cc}.txt` — categorical covariates (all `--h5ad-cat-covars`, as factors)
- `cov2_{cc}_reduced.txt` — quantitative covariates (expression PCs 1–20, G-PCs 1–5, SVs, min-max scaled `--h5ad-quant-covars`, cell composition)

---

### Step 3: `03-build-eqtl.sh`

**Purpose:** Run OSCA eQTL analysis for a single cell class (called in parallel for all cell classes).

**Usage** (called by orchestrator via GNU parallel):
```bash
03-build-eqtl.sh \
  Phenotype_{cc}_osca.txt \
  befile_{cc} \
  {vcf_slogan} \
  Upprobe_{cc}.opi \
  cov1_{cc}.txt \
  cov2_{cc}_reduced.txt \
  5 \
  {osca_output_dir}/eqtl_{cc}.tsv
```

**OSCA sub-steps:**
1. `--make-bod` — encode expression phenotype as binary BOD
2. `--update-opi` — update probe info with gene locations
3. `--eqtl --cis --cis-wind 1000` — cis-eQTL test within 1 Mb window (`--cis-wind` unit is Kb; hardcoded; see Known Issues #6)
4. `--query 1 --to-smr` — export all cis results as TSV

**Outputs:**
- `tempeqtl_{cc}.tsv` → renamed to `eqtl_{cc}.tsv` after query
- `befile_{cc}.*` — OSCA BOD files
- `eqtl_{cc}_osca_{output,error,progress}.log`

---

### Step 4: `gcloud storage cp` (optional)

Copies all `.tsv` files from the output directory to:
`gs://macosko_data/ferris/eqtl_outs/{pipeline_slogan}/`

Skip with `--start-at-step 5` if GCS backup is not needed.

---

### Step 5: `05-process-and-plot-osca-tsv.R`

**Purpose:** Apply two-stage FDR to OSCA output; save as RDS; generate Manhattan and histogram plots.

**Usage** (called via GNU parallel):
```bash
Rscript 05-process-and-plot-osca-tsv.R --path eqtl_{cc}.tsv [--out-dir DIR]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--path` | *(required)* | Path to `eqtl_{cc}.tsv` |
| `--out-dir` | `dirname(path)` | Output directory for RDS and plots |

**FDR procedure (two-stage):**
1. Within each gene (Probe): Bonferroni correction across all its cis-SNPs → `p_bonf`
2. Take minimum `p_bonf` per gene → BH correction across genes → `padj_gene`
3. eGenes: those with `padj_gene < 0.05`
4. Within eGenes: BH correction across all SNPs → `padj_snp`
5. Significant eSNPs: `padj_gene < 0.05 AND padj_snp < 0.05`

**Constants:** `FDR_THRESH = 0.05`, `HARD_CAP = 40` (capped -log10 p for plotting)

**Outputs:**
- `eqtl_{cc}.rds` — full data frame with FDR columns
- `eqtl_{cc}_sig.rds` — significant rows only
- `plots/eqtl_{cc}_manhattan.png`
- `plots/eqtl_{cc}_min_p_gene_hist.png`
- `plots/eqtl_{cc}_padj_gene_hist.png`
- `plots/eqtl_{cc}_pval_hist.png`
- `plots/eqtl_{cc}_padj_snp_hist.png`

**cell_class extraction:** `cell_class = sub("^eqtl_", "", tools::file_path_sans_ext(basename(path)))`. Handles underscores within cell class names (e.g. `eqtl_da_neuron.tsv` → `da_neuron`). Requires input file named `eqtl_{cell_class}.tsv`.

**Note:** Chromosomes 1–22 assumed in Manhattan plot; others emit a warning (see Known Issues #5).

---

### Step 6: `06-get-common-snp-probes-osca-rds.R`

**Purpose:** Find SNP-probe pairs tested in all cell classes; identify those significant in at least one.

**Usage:**
```bash
Rscript 06-get-common-snp-probes-osca-rds.R --base DIR
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--base` | *(required)* | Directory containing per-cell-class `.rds` files (output of step 5) |

**Processing:**
1. Reads all `*.rds` files in `--base` (excluding `_sig.rds`)
2. Adds `cell_class` column and `SNP_probe` composite key
3. Intersects `SNP_probe` sets across all cell classes
4. Subsets all data to common pairs
5. Saves combined long-format RDS
6. Saves subset significant in at least one class (threshold: `padj_snp < 0.01` — **see Known Issues #3**)

**Outputs:**
- `eqtl_present_in_all.rds`
- `eqtl_present_in_all_sig_in_one.rds`

---

### Step 7: `07-run-mashr.R`

**Purpose:** Estimate sharing of eQTL effect sizes across cell classes using mashr.

**Usage:**
```bash
Rscript 07-run-mashr.R \
  --path eqtl_present_in_all.rds \
  --padj-thresh 0.01 \
  --num-random 1000000 \
  --eps 1e-6
```

| Argument | Default | Description |
|----------|---------|-------------|
| `--path` | *(required)* | Path to `eqtl_present_in_all.rds` |
| `--padj-thresh` | `0.01` | BH adjusted p-value threshold for "strong" tests |
| `--num-random` | `1000000` | Random tests for null correlation estimation |
| `--eps` | `1e-6` | Floor for SE and near-zero beta values |

**mashr workflow:**
1. Pivot long → wide (SNP-probe × cell-class matrices for b, SE, padj_snp)
2. Clamp SE and beta to `±eps`
3. Identify "strong" tests: `min(padj_snp across classes) < padj-thresh`
4. Draw random subset of size `num-random`
5. Estimate null correlation (`estimate_null_correlation_simple`) from random subset
6. Compute PCA + empirical covariances (`cov_pca`, `cov_ed`) from strong tests
7. Add canonical covariances + custom neuron/glia covariances — **hard-coded; see Known Issues #1**
8. Fit mash on random data; apply to strong tests with fixed mixture weights

**Outputs:**
- `eqtl_present_in_all__mash_results_sig.rds` — mashr posteriors for strong tests
- `eqtl_present_in_all__mash_results_random.rds` — fitted mash model
- `eqtl_present_in_all__mash_data_sig.rds` — mash data object (strong)
- `eqtl_present_in_all__mash_data_random.rds` — mash data object (random)

---

### Step 8: `gcloud storage cp -r` (optional)

Archives the entire `--osca-input-dir` (inputs + all outputs) to:
`gs://macosko_data/ferris/eqtl_outs/{pipeline_slogan}/`

---

## Outputs Reference

| File | Step | Format | Description |
|------|------|--------|-------------|
| `{stem}_expression_matrix_ds.csv` | 1 | CSV | Pseudobulk expression (samples × genes) |
| `{stem}_composition_matrix_ds.csv` | 1 | CSV | Cell composition per sample |
| `Phenotype_{cc}_osca.txt` | 2 | TSV | OSCA phenotype file (genes × participants) |
| `Upprobe_{cc}.opi` | 2 | TSV | OSCA probe info |
| `cov1_{cc}.txt` | 2 | TSV | Categorical covariates |
| `cov2_{cc}_reduced.txt` | 2 | TSV | Quantitative covariates (SVA, PCA, genotype PCs) |
| `befile_{cc}.*` | 3 | binary | OSCA BOD files |
| `eqtl_{cc}.tsv` | 3 | TSV | Raw OSCA output (deleted after step 5) |
| `eqtl_{cc}.rds` | 5 | RDS | Full results with FDR columns |
| `eqtl_{cc}_sig.rds` | 5 | RDS | Significant results (padj_gene & padj_snp < 0.05) |
| `plots/eqtl_{cc}_manhattan.png` | 5 | PNG | Manhattan plot |
| `plots/eqtl_{cc}_{min_p,padj_gene,pval,padj_snp}_hist.png` | 5 | PNG | Diagnostic histograms |
| `eqtl_present_in_all.rds` | 6 | RDS | Long format, common SNP-probe pairs |
| `eqtl_present_in_all_sig_in_one.rds` | 6 | RDS | Pairs significant in ≥ 1 class |
| `eqtl_present_in_all__mash_results_sig.rds` | 7 | RDS | mashr posteriors (strong tests) |
| `eqtl_present_in_all__mash_results_random.rds` | 7 | RDS | mashr fitted model |
| `eqtl_present_in_all__mash_data_{sig,random}.rds` | 7 | RDS | mash_set_data objects |

RDS column schema (steps 5–7):

| Column | Description |
|--------|-------------|
| `SNP` | SNP identifier |
| `Probe` | Gene/probe identifier |
| `Chr` | Chromosome |
| `BP` | Base pair position |
| `b` | Effect size (beta) |
| `SE` | Standard error |
| `p` | Nominal p-value |
| `p_bonf` | Within-gene Bonferroni p |
| `padj_gene` | BH-adjusted gene-level p |
| `padj_snp` | BH-adjusted SNP-level p (within eGenes) |
| `is_significant_gene` | `padj_gene < 0.05` |
| `is_significant_snp` | `padj_snp < 0.05` |
| `cell_class` | Cell class label (added in step 6) |
| `SNP_probe` | `{SNP}_{Probe}` composite key (added in step 6) |

---

## Partial Execution

Use `--start-at-step` and `--stop-after-step` to run subsets of the pipeline:

| Use case | Command flags |
|----------|--------------|
| Pseudobulk only | `--stop-after-step 1` |
| Format OSCA inputs only | `--start-at-step 2 --stop-after-step 2` |
| Run OSCA from scratch, skip pseudobulk | `--start-at-step 2` |
| Just run OSCA (inputs already formatted) | `--start-at-step 3 --stop-after-step 3` |
| Process outputs only (OSCA done) | `--start-at-step 5 --stop-after-step 7` |
| Re-run mashr only | `--start-at-step 7` |
| Skip GCS backups (steps 4 and 8) | `--start-at-step 1 --stop-after-step 7` |

Note: step 1 uses `--start-at-step` without `--stop-after-step` check (always runs if `START_AT_STEP <= 1`). All other steps check both flags.

---

## Known Issues & Limitations

### 1. mashr hard-coded cell classes (`07-run-mashr.R`, lines 156–182)

Custom covariance structures for "neurons only" and "glia only" groups are built by hard-coding column names `astro`, `da`, `mg`, `nonda`, `oligo`:

```r
astro_col_ind = which(colnames(eqtl_wide_b) == "astro")
...
neurons_only[da_col_ind, da_col_ind] = 1
```

**Risk:** If any expected cell class is absent, `which()` returns `integer(0)`, causing silent incorrect matrix indexing. If cell classes differ from these five, the covariances will be wrong.

**Patched version:** `/home/ferris/cc-sandbox/osca_eqtl_pipeline_patches/07-run-mashr.R`

---

### 2. Hard-coded SVA formula (`02-run-osca-formatting-scanpy.R`, line 94)

SVA always uses:
```r
mod = model.matrix(~ age + sex + pmi, data = pheno)
```

This is not configurable; different studies with different covariates must edit the script.

**Patched version:** `/home/ferris/cc-sandbox/osca_eqtl_pipeline_patches/02-run-osca-formatting-scanpy.R`
Adds `--sva-formula` argument (default `"~ case_control + age + sex + study + brain_bank"`).

---

### 3. Inconsistent FDR thresholds between steps 5 and 6

Step 5 (`05-process-and-plot-osca-tsv.R`) uses `FDR_THRESH = 0.05` for `padj_snp` significance.
Step 6 (`06-get-common-snp-probes-osca-rds.R`, line 63) uses `padj_snp < 0.01` for the `sig_in_one` subset.

This means the "significant in one" file uses a stricter threshold than step 5's significance calls, which can be confusing.

**Patched version:** `/home/ferris/cc-sandbox/osca_eqtl_pipeline_patches/06-get-common-snp-probes-osca-rds.R`
Adds `--padj-thresh` argument (default `0.05` to match step 5).

---

### 4. Hard-coded metadata path (`02-run-osca-formatting-scanpy.R`)

The default `--metadata` path is:
```
/mnt/accessory/seq_data/pd_all/240514/dapi_nurr_merged_seurat_clean_subsets/dapi_nurr_metadata_new.csv
```
This is project-specific. Always pass `--metadata` explicitly, or the script will silently use stale data.

The project-specific prefix strip (`gsub("dapi_nurr_merged_seurat_clean__", ...)`) has been removed from the patched version. Output CSVs from step 1 are now named `{cell_class}_expression_matrix_ds.csv` where `cell_class` comes from `adata.obs[CT_ID]`, so no stripping is needed.

---

### 5. Chromosomes 1–22 only (`05-process-and-plot-osca-tsv.R`, line 248)

The Manhattan plot orders chromosomes as `1:22`. Sex chromosomes (X, Y) and mitochondrial DNA are excluded from the plot layout.

---

### 6. `CIS_WIND=1000` not parameterized (`03-build-eqtl.sh`)

`CIS_WIND=1000` is hardcoded. The `--cis-wind` flag unit is **Kb**, so 1000 = **1 Mb** (a standard cis-eQTL window; OSCA's default is 2000 Kb = 2 Mb). This is documented in the patched script with a comment. Verify this is intentional for your use case, or edit the variable.

---

## Logging

The orchestrator redirects stdout and stderr to:
```
/mnt/accessory/analysis/eqtl/logs/osca_eqtl_pipeline_{pipeline_slogan}_{datetime}.log
```

Each cell class in step 3 writes its own logs:
```
{osca_output_dir}/eqtl_{cc}_osca_{output,error,progress}.log
```

The git branch and commit hash are logged at startup for reproducibility.
