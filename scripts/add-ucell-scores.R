# path to a tsv file with columns "score_name" and "genes". 
# The genes column should be in the format of a comma-delimited string.

print("Loading Libraries")
suppressWarnings(suppressMessages(library(getopt)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(qs)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(UCell)))

spec <- matrix(c(
    'sobj-path', 'p', 1, "character",
    'score-path', 's', 1, "character",
    'assay', 'a', 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)

SOBJ_PATH = opt[["sobj-path"]]
SCORE_PATH = opt[["score-path"]]
if(is.null(opt[['assay']])) {
    ASSAY = "RNA"
} else {
    ASSAY = opt[["assay"]]
}

print( paste0("Reading gene scores from ", SCORE_PATH))
score_list_df = read.table(SCORE_PATH, header=TRUE, sep="\t", stringsAsFactors=FALSE)
print(head(score_list_df))
gene_sets = lapply(score_list_df$genes, function(genes) {
    strsplit(genes, ",")[[1]]
})
names(gene_sets) = score_list_df$score_name

print( paste("Reading Seurat object from", SOBJ_PATH))
sobj = qread(SOBJ_PATH)

for (name in names(gene_sets)){
    gene_sets[[name]] = gene_sets[[name]][gene_sets[[name]] %in% rownames(sobj@assays[[ASSAY]]$counts)]
    print(paste0("Calculating score for ", name))
    print("Genes:")
    print(paste(gene_sets[[name]][1:5], collapse=', '))
    print("...")
}
print(sapply(gene_sets, length))

missing_genes <- setdiff(unlist(gene_sets), rownames(sobj@assays[[ASSAY]]$counts))
if (length(missing_genes) > 0) {
    print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
}

print("Calculating U Scores")
u_scores = ScoreSignatures_UCell(
    as.matrix(
        GetAssayData(sobj, slot="counts", assay = ASSAY),
        gene_sets, 
        assay = ASSAY)

u_scores = u_scores[match(rownames(sobj@meta.data), rownames(u_scores)),]
u_scores = as.data.frame(u_scores)

if (all(rownames(u_scores) == rownames(sobj@meta.data))){
    for (score_name in colnames(u_scores)){
        sobj@meta.data[[paste0(score_name, "__", ASSAY)]] = u_scores[[score_name]]
    }
}

print("Saving Seurat object")
qsave(sobj, SOBJ_PATH)