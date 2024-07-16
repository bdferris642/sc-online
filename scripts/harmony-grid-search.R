library(Seurat)
library(dplyr)
library(qs)
library(harmony)

BASE_PATH="/mnt/accessory/seq_data/pd_all/20240410"
OUT_PATH = file.path(BASE_PATH, "harmony")
dir.create(OUT_PATH)

THETAS = c(0.5, 1, 2, 4, 6, 8)
SIGMAS = c(0.05, 0.1, 0.2, 0.4, 0.8)


nurr_da_merged = qread(file.path(BASE_PATH, "seurat_nurr_da_merged_20240424.qs"))

for (this_theta in THETAS){
    for (this_sigma in SIGMAS){

        suffix = paste0("theta_", this_theta, "_sigma_", this_sigma)
        nurr_da_merged = (nurr_da_merged
            %>% RunHarmony(
                group.by.vars="batch",
                dims.use = 1:30,
                early_stop = T,
                max_iter=25,
                reference_values = "gte100_ctr",
                plot_convergence = F,
                theta=this_theta,
                sigma=this_sigma)
            %>% FindNeighbors(dims = 1:30, reduction="harmony")
            %>% FindClusters(resolution = 0.2)
            %>% FindClusters(resolution = 0.4)
            %>% FindClusters(resolution = 0.5)
            %>% FindClusters(resolution = 0.6)
            %>% FindClusters(resolution = 0.8)
            %>% FindClusters(resolution = 1)
            %>% RunUMAP(dims = 1:30, verbose=F, reduction="harmony")
        )

        qsave(nurr_da_merged, file.path(OUT_PATH, paste0("seurat_nurr_da_merged_20240424_harmony_", suffix, ".qs")))
    }
}
