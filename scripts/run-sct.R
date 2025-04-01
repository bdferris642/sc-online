library(dplyr)
library(getopt)
library(qs)
library(Seurat)
library(sctransform)

spec <- matrix(c(
    'path', 'p', 1, "character",
    'workers', 'w', 1, "numeric",
    'max-size-mb', 'm', 1, "numeric"
), byrow = TRUE, ncol = 4)

opt = getopt(spec)
PATH = opt[['path']]
if (is.null(opt[['workers']])) {
    WORKERS = 4
} else {
    WORKERS = opt[['workers']]
}

if (is.null(opt[['max-size-mb']])) {
    MAX_SIZE = 80000 * 1024^2
} else {
    MAX_SIZE = opt[['max-size-mb']] * 1024^2
}

# Load the Seurat object
print(paste("Reading Seurat Object", basename(PATH), "from", dirname(PATH)))
sobj = qread(PATH)

# embiggen resources
plan(strategy = "multicore", workers = WORKERS)
options(future.globals.maxSize = MAX_SIZE)

# Run SCTransform
print("Running SC Transform")
sobj = SCTransform(sobj, verbose=TRUE)

print(paste("Saving back to", PATH))
qsave(sobj, PATH)


