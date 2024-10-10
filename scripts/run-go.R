library(qs)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)

base="/mnt/accessory/seq_data/pd_all/240514"

run_fisher_exact_test = function(list_A, list_B, total_genes, alternative = "greater") {
    # Calculate the overlap and non-overlap
    a <- length(intersect(list_A, list_B))
    b <- length(setdiff(list_A, list_B))
    c <- length(setdiff(list_B, list_A))
    d <- total_genes - (a + b + c)

    # Create the contingency table
    contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    rownames(contingency_table) <- c("In_List_B", "Not_in_List_B")
    colnames(contingency_table) <- c("In_List_A", "Not_in_List_A")

    # Perform Fisher's exact test
    fisher_test_result <- fisher.test(contingency_table, alternative = alternative)

    # Display the results
    return(fisher_test_result)
}


nalls = read.csv("~/nalls-2019-s2.csv", header=TRUE, skip=1)
nalls = nalls[nalls$Failed.final.filtering.and.QC == 0,]
nalls$gene_coalesced = ifelse(
    nalls$QTL.Nominated.Gene..nearest.QTL %in% c("", NA, NULL, "NA"),
    nalls$Nearest.Gene, 
    nalls$QTL.Nominated.Gene..nearest.QTL)

nalls = nalls[!duplicated(nalls$gene_coalesced),]
nalls = nalls[!is.na(nalls$gene_coalesced) & !nalls$gene_coalesced %in% c("", "NA", "NULL"),]



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
marker_df = do.call(rbind, marker_df_list)

cts = sort(unique(marker_df$cluster))

for (i in 1:length(cts)) {
    for (j in (i+1):length(cts)) {
        ct1 = cts[[i]]
        ct2 = cts[[j]]
        if (ct1 == ct2) {
            next
        }

        markers_ct1 = marker_df[marker_df$cluster == ct1,] 
        markers_ct2 = marker_df[marker_df$cluster == ct2,]

        cc1 = unique(markers_ct1$cell_class)[[1]]
        
        # OJO! TODO! SKIPPING ASTRO FOR NOW ONLY
        if (cc1 == "astro") {
            next         
        }
        
        cc2 = unique(markers_ct2$cell_class)[[1]]

        bkg_genes_cc1 = bkg_genes[[cc1]]
        bkg_genes_cc2 = bkg_genes[[cc2]]

        bkg_genes_ct1_ct2 = intersect(bkg_genes_cc1, bkg_genes_cc2)
        print(paste("Number of background genes for", paste0(ct1, "_x_", ct2), "is", length(bkg_genes_ct1_ct2)))
        if (length(bkg_genes_ct1_ct2) == 0) {
            next
        }

        ct1_x_ct2 = intersect(markers_ct1$gene, markers_ct2$gene)
        print(paste("Number of markers for", paste0(ct1, "_x_", ct2), "is", length(ct1_x_ct2)))
        if (length(ct1_x_ct2) == 0) {
            next
        }

        ct1_x_ct2_str = paste0(ct1, "_x_", ct2)

        print(paste("Runnng Go Enrichment for ", ct1_x_ct2_str))
        # Perform GO enrichment analysis
        bp = enrichGO(
            gene          = ct1_x_ct2,
            universe      = bkg_genes_ct1_ct2,
            OrgDb         = org.Hs.eg.db,      # Replace with appropriate organism database
            keyType       = "SYMBOL",          # Or "ENTREZID" depending on your gene IDs
            ont           = "BP",              # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
            pAdjustMethod = "BH",              # Adjust for multiple testing using Benjamini-Hochberg
            qvalueCutoff  = 0.05,              # Adjust significance level if needed
            readable      = TRUE               # Convert results to readable gene symbols
        )

        mf = enrichGO(
            gene          = ct1_x_ct2,
            universe      = bkg_genes_ct1_ct2,
            OrgDb         = org.Hs.eg.db,      # Replace with appropriate organism database
            keyType       = "SYMBOL",          # Or "ENTREZID" depending on your gene IDs
            ont           = "MF",              # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
            pAdjustMethod = "BH",              # Adjust for multiple testing using Benjamini-Hochberg
            qvalueCutoff  = 0.05,              # Adjust significance level if needed
            readable      = TRUE               # Convert results to readable gene symbols
        )

        cc = enrichGO(
            gene          = ct1_x_ct2,
            universe      = bkg_genes_ct1_ct2,
            OrgDb         = org.Hs.eg.db,      # Replace with appropriate organism database
            keyType       = "SYMBOL",          # Or "ENTREZID" depending on your gene IDs
            ont           = "CC",              # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
            pAdjustMethod = "BH",              # Adjust for multiple testing using Benjamini-Hochberg
            qvalueCutoff  = 0.05,              # Adjust significance level if needed
            readable      = TRUE               # Convert results to readable gene symbols
        )

        print(paste("Saving results for", ct1_x_ct2_str))
        qsave(bp, file.path(base, paste0("markers/go/", ct1_x_ct2_str, "_go_bp_enrichment.qs")))
        qsave(mf, file.path(base, paste0("markers/go/", ct1_x_ct2_str, "_go_mf_enrichment.qs")))
        qsave(cc, file.path(base, paste0("markers/go/", ct1_x_ct2_str, "_go_cc_enrichment.qs")))
    }
}
