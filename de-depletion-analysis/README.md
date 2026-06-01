# de-depletion-analysis

Correlates per-gene DE statistics (and per-pathway GSEA NES values) with
cell-type depletion estimates from crumblr, across multiple cell types.
Produces scatter plots, volcano plots, and heatmaps — then calls `run-gsea.R`
on the gene r-values to run pathway enrichment on the correlation results.

---

## Setup

Run once from any directory:

```bash
bash /path/to/cc-sandbox/de-depletion-analysis/setup.sh
```

Creates the `de-depletion-env` conda environment inside
`cc-sandbox/micromamba_root/`. Prints "All R packages OK" on success.

---

## Usage

The script self-activates its environment — no manual `conda activate` needed:

```bash
Rscript cc-sandbox/de-depletion-analysis/run_depletion_analysis.R \
  --crumblr   /path/to/crumblr_results.csv \
  --manifest  /path/to/manifest.csv \
  --output-dir /path/to/output/
```

---

## Manifest format

CSV with a header row. Required columns:

| Column | Description |
|--------|-------------|
| `ct_name` | Cell type name — must match `cell_type` in the crumblr CSV |
| `de_csv_path` | Absolute path to the DE results CSV for this cell type |
| `gsea_csv_path` | Absolute path to the GSEA results CSV for this cell type |
| `de_rank_col` | Column in the DE CSV to use as the ranking/correlation statistic (e.g. `z_stat`, `logFC`) |
| `gsea_rank_col` | Column in the GSEA CSV to use as the NES value (default: `NES` if empty/absent) |

The DE CSV must contain a gene column named `gene_short_name`, `gene`, or
`gene_id` (tried in that order).

The GSEA CSV must contain a pathway column named `pathway`, `Pathway`, or
`name`; and optionally a `gene_set` column (used to create unique pathway IDs
across databases).

### Example manifest

```csv
ct_name,de_csv_path,gsea_csv_path,de_rank_col,gsea_rank_col
DAN,/data/dan_de.csv,/data/dan_gsea.csv,z_stat,NES
OPC,/data/opc_de.csv,/data/opc_gsea.csv,z_stat,NES
Astrocyte,/data/ast_de.csv,/data/ast_gsea.csv,z_stat,
```

---

## CLI reference

| Flag | Default | Description |
|------|---------|-------------|
| `--crumblr` | required | crumblr results CSV |
| `--manifest` | required | Manifest CSV |
| `--output-dir` | required | Output directory (created if missing) |
| `--min-de-sets` | 5 | Gene/pathway must appear in ≥N files to enter universe |
| `--top-n-genes` | 50 | Top genes by \|r\| shown in heatmap |
| `--top-n-pathways` | 50 | Top pathways by \|r\| shown in heatmap |
| `--focus-ct` | NULL | Partial string match for focus CT; default = most depleted CT |
| `--n-perm` | 10000 | Permutations for empirical p-value |
| `--n-procs` | 60 | Parallel cores for gene permutation loop |
| `--force` | FALSE | Use intersection on manifest/crumblr mismatch instead of failing |
| `--hsp-genes` | `/mnt/accessory/seq_data/pd-freeze/HSP-genes.txt` | One gene per line |
| `--pd-gwas-genes` | `/mnt/accessory/seq_data/pd-freeze/GP2-pd-genes.txt` | One gene per line |
| `--n-label` | 15 | Genes/pathways labelled per direction on scatter plots |
| `--padj-threshold` | 0.05 | FDR threshold for heatmap border significance marks |

---

## Output file reference

All outputs are written to `--output-dir`. Each plot has a paired `.csv` with
the plotted data.

### Union set (all genes / pathways)

| File pattern | Description |
|-------------|-------------|
| `gene_depletion_correlation_scatter.{svg,csv}` | All genes; X=r, Y=mean\|stat\|; used as input for run-gsea.R |
| `gene_depletion_volcano.{svg,csv}` | All genes; X=r, Y=-log10(padj) |
| `gene_depletion_heatmap.{svg,csv}` | Top genes by \|r\|; rows=genes, cols=CTs |
| `gene_depletion_correlation_scatter_hsp.{svg,csv}` | HSP gene subset |
| `gene_depletion_volcano_hsp.{svg,csv}` | HSP gene subset |
| `gene_depletion_heatmap_hsp.{svg,csv}` | HSP gene subset |
| `gene_depletion_correlation_scatter_pdgwas.{svg,csv}` | PD GWAS gene subset |
| `gene_depletion_volcano_pdgwas.{svg,csv}` | PD GWAS gene subset |
| `gene_depletion_heatmap_pdgwas.{svg,csv}` | PD GWAS gene subset |
| `pathway_depletion_correlation_scatter.{svg,csv}` | All pathways; faceted by gene_set |
| `pathway_depletion_heatmap.{svg,csv}` | Top pathways by \|r\|; faceted by gene_set |

### Focus set (genes/pathways sig in focus CT or by global padj)

Same files as above with `_focus_{ct_name}` suffix, plus `_hsp` / `_pdgwas`
variants. The focus CT defaults to the most depleted cell type (lowest crumblr
logFC); override with `--focus-ct`.

### Other outputs

| File | Description |
|------|-------------|
| `run.log` | Timestamped log of all analysis steps |
| `gsea/` | GSEA results from run-gsea.R (placed by that script in dirname of scatter CSV) |

---

## Color scheme

- **Blue–gray–red** (`FILL_R`): used for all DE stat / NES / r values.
  Red = positive r = gene more DE in depleted cell types.
  Blue = negative r = gene more DE in resilient (non-depleted) cell types.
  Gray = center (neutral). Palette: `#2166ac → #92c5de → #d4d4d4 → #f4a582 → #b2182b`.

- **Purple–white–green** (`FILL_LFC`): used for the crumblr logFC annotation
  strip at the top of each heatmap. Green = depleted (negative logFC in crumblr
  convention); purple = enriched.

---

## Example command

```bash
Rscript /home/ferris/cc-sandbox/de-depletion-analysis/run_depletion_analysis.R \
  --crumblr    /home/ferris/cc-sandbox/dan-analysis/crumblr_results.csv \
  --manifest   /home/ferris/cc-sandbox/de-depletion-analysis/test_manifest.csv \
  --output-dir /home/ferris/cc-sandbox/de-depletion-analysis/test_output/ \
  --min-de-sets 3 \
  --n-perm 1000 \
  --n-procs 8
```
