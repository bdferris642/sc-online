library(qs)
library(clusterProfiler)
library(org.Hs.eg.db)

# Function to write commands to a file
write_command <- function(gene_set, background, output_file, ontology, gene_set_file, background_file) {
  cat(sprintf(
    "Rscript -e \"library(qs); library(clusterProfiler); library(org.Hs.eg.db); genes <- readRDS('%s'); background <- readRDS('%s'); enrich_result <- enrichGO(gene=genes, universe=background, OrgDb=org.Hs.eg.db, keyType='SYMBOL', ont='%s', pAdjustMethod='BH', qvalueCutoff=0.05, readable=TRUE); qsave(enrich_result, '%s')\"\n",
    gene_set_file,
    background_file,
    ontology,
    output_file
  ), file="~/go-commands.txt", append=TRUE)
}
# Generate commands for parallel execution
file.remove("~/go-commands.txt")


base="/mnt/accessory/seq_data/pd_all/240514"
go_subdir = "markers/go/ctr"

astro_markers = read.csv(file.path(base, "markers/dapi_astro_clean_ctr_markers_cell_type.csv"))
endo_markers = read.csv(file.path(base, "markers/dapi_endo_clean_ctr_markers_cell_type.csv"))
immune_markers = read.csv(file.path(base, "markers/dapi_immune_clean_ctr_markers_cell_type.csv"))
mg_markers = read.csv(file.path(base, "markers/dapi_mg_clean_ctr_markers_cell_type.csv"))
myeloid_markers = read.csv(file.path(base, "markers/dapi_myeloid_ctr_markers_cell_type.csv"))
opc_markers = read.csv(file.path(base, "markers/dapi_opc_seurat_ctr_markers_cell_type.csv"))
oligo_markers = read.csv(file.path(base, "markers/dapi_oligo_seurat_ctr_markers_cell_type.csv"))

inh_markers = read.csv(file.path(base, "markers/dapi_nurr_merged_seurat_clean_inh_ctr_markers_cell_type.csv"))
da_markers = read.csv(file.path(base, "markers/dapi_nurr_merged_seurat_clean_da_ctr_markers_cell_type.csv"))
ex_markers = read.csv(file.path(base, "markers/dapi_nurr_merged_seurat_clean_ex_ctr_markers_cell_type.csv"))

astro_markers$cell_class = "astro"
endo_markers$cell_class = "endo"
immune_markers$cell_class = "immune"
mg_markers$cell_class = "mg"
myeloid_markers$cell_class = "myeloid"
opc_markers$cell_class = "opc"
oligo_markers$cell_class = "oligo"
inh_markers$cell_class = "inh"
da_markers$cell_class = "da"
ex_markers$cell_class = "ex_n"
inh_markers$cell_class = "inh_n"

bkg_gene_filenames = list.files(file.path(base, "genes_expressed"), full.names = TRUE)
names(bkg_gene_filenames) = sapply(bkg_gene_filenames, function(x) {
    b = basename(x)
    cc = strsplit(b, "__")[[1]][[2]]
    return(cc)
})

bkg_genes= list()
for (cc in names(bkg_gene_filenames)) {
    bkg_genes[[cc]] = qread(bkg_gene_filenames[[cc]])
}

marker_df_list = list(
    astro = astro_markers,
    endo = endo_markers,
    immune = immune_markers,
    mg = mg_markers,
    myeloid = myeloid_markers,
    opc = opc_markers,
    oligo = oligo_markers,
    inh = inh_markers,
    da = da_markers,
    ex_n = ex_markers
)

for (name in names(marker_df_list)){
    this_marker_df = marker_df_list[[name]]
    this_marker_df = this_marker_df[this_marker_df$p_val_adj < 0.05,]
    marker_df_list[[name]] = this_marker_df
}
marker_df = do.call(rbind, marker_df_list)
cts = sort(unique(marker_df$cluster))

for (i in 1:length(cts)) {
    ct = cts[[i]]
    markers_ct = marker_df[marker_df$cluster == ct,] 
    cell_class = unique(markers_ct$cell_class)[[1]]
    background_genes = bkg_genes[[cell_class]]
    print(paste("Number of markers for", ct, "is", nrow(markers_ct)))
    if (nrow(markers_ct) == 0) {
        next
    }
    print(paste("Number of background genes for", ct, "is", length(background_genes)))
    if (length(background_genes) == 0) {
        next
    }

    # Save the gene set and background to files
    gene_set = markers_ct$gene
    gene_set_file <- file.path(base, paste0(go_subdir, "/inputs/", ct, "_gene_set.rds"))
    if (!file.exists(gene_set_file)){
        saveRDS(gene_set, gene_set_file)
    }
        
    background_file <- file.path(base, paste0(go_subdir, "/inputs/", cell_class, "_background.rds"))
    if (!file.exists(background_file)){
        saveRDS(background_genes, background_file)
    }

    # Write commands for BP, MF, and CC ontologies
    write_command(gene_set, background_genes, 
        file.path(base, paste0(go_subdir, "/outs/", ct, "__go_bp_enrichment.qs")), "BP", gene_set_file, background_file)
    write_command(gene_set, background_genes, 
        file.path(base, paste0(go_subdir, "/outs/", ct, "__go_mf_enrichment.qs")), "MF", gene_set_file, background_file)
    write_command(gene_set, background_genes, 
        file.path(base, paste0(go_subdir, "/outs/", ct, "__go_cc_enrichment.qs")), "CC", gene_set_file, background_file)
}


# Execute the commands in parallel using GNU parallel
#system("parallel -j 30 < ~/go-commands.txt")
