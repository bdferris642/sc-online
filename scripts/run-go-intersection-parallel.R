library(qs)
library(clusterProfiler)
library(org.Hs.eg.db)

# change depending on number of cores on your VM
ncores=30

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

for (i in 1:(length(cts) - 1)) {
    for (j in (i+1):length(cts)) {
        ct1 = cts[[i]]
        ct2 = cts[[j]]
        if (ct1 == ct2) {
            next
        }

        markers_ct1 = marker_df[marker_df$cluster == ct1,] 
        markers_ct2 = marker_df[marker_df$cluster == ct2,]

        cell_class_1 = unique(markers_ct1$cell_class)[[1]]
        cell_class_2 = unique(markers_ct2$cell_class)[[1]]

        background_genes_1 = bkg_genes[[cell_class_1]]
        background_genes_2 = bkg_genes[[cell_class_2]]

        background_genes = intersect(background_genes_1, background_genes_2)
        print(paste("Number of background genes for", paste0(ct1, "_x_", ct2), "is", length(background_genes)))
        if (length(background_genes) == 0) {
            next
        }

        gene_set = intersect(markers_ct1$gene, markers_ct2$gene)
        print(paste("Number of markers for", paste0(ct1, "_x_", ct2), "is", length(gene_set)))
        if (length(gene_set) == 0) {
            next
        }

        comparison_str = paste0(ct1, "_x_", ct2)
        cell_class_str = paste0(cell_class_1, "_x_", cell_class_2)

        # Save the gene set and background to files
        gene_set_file <- file.path(base, paste0("markers/go/inputs/", comparison_str, "_gene_set.rds"))
        if (!file.exists(gene_set_file)){
            saveRDS(gene_set, gene_set_file)
        }
        
        background_file <- file.path(base, paste0("markers/go/inputs/", cell_class_str, "_background.rds"))
        if (!file.exists(background_file)){
            saveRDS(background_genes, background_file)
        }

        # Write commands for BP, MF, and CC ontologies
        write_command(gene_set, background_genes, 
            file.path(base, paste0("markers/go/outs", comparison_str, "_go_bp_enrichment.qs")), "BP", gene_set_file, background_file)
        write_command(gene_set, background_genes, 
            file.path(base, paste0("markers/go/outs", comparison_str, "_go_mf_enrichment.qs")), "MF", gene_set_file, background_file)
        write_command(gene_set, background_genes, 
            file.path(base, paste0("markers/go/outs", comparison_str, "_go_cc_enrichment.qs")), "CC", gene_set_file, background_file)
    }
}

# Execute the commands in parallel using GNU parallel
system(paste("parallel -j", ncores, "< ~/go-commands.txt"))
