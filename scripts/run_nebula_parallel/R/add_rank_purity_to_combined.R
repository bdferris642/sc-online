#!/usr/bin/env Rscript
# Left-join rank_purity.py output (gene, rank_purity) into a combined NEBULA result and rewrite
# both the .qs (summary element) and the companion .csv, in the format combine-nebula-results.R uses.

# ---- self-activation: re-exec under the workflow env Rscript (env path derived from this file)
local({
    args_all = commandArgs(trailingOnly = FALSE)
    script = normalizePath(sub("^--file=", "", args_all[grep("^--file=", args_all)]))
    env_prefix = file.path(dirname(dirname(dirname(script))), "micromamba_root", "envs", "run-nebula-parallel-env")
    if (!dir.exists(env_prefix)) stop("workflow env not found: ", env_prefix, " (run run_nebula_parallel/setup.sh)")
    if (normalizePath(R.home()) != normalizePath(file.path(env_prefix, "lib", "R"))) {
        status = system2(file.path(env_prefix, "bin", "Rscript"), c(shQuote(script), shQuote(commandArgs(trailingOnly = TRUE))))
        quit(status = status, save = "no")
    }
})

suppressMessages(suppressWarnings({
    library(getopt)
    library(qs)
    library(dplyr)
}))

spec = matrix(c(
    "combined-qs", "q", 1, "character",
    "purity-tsv",  "p", 1, "character"
), byrow = TRUE, ncol = 4)
opt = getopt(spec)
if (is.null(opt[["combined-qs"]]) || is.null(opt[["purity-tsv"]])) stop("Required: --combined-qs --purity-tsv")

combined_qs = opt[["combined-qs"]]
combined_csv = sub("\\.qs$", ".csv", combined_qs)

combined = qread(combined_qs)
purity = read.delim(opt[["purity-tsv"]], stringsAsFactors = FALSE)
if (!all(c("gene", "rank_purity") %in% colnames(purity))) stop("purity TSV needs columns gene, rank_purity")
if ("rank_purity" %in% colnames(combined$summary)) stop("combined summary already has a rank_purity column")

n_before = nrow(combined$summary)
combined$summary = combined$summary %>% left_join(purity, by = "gene")
if (nrow(combined$summary) != n_before) stop("join changed the number of summary rows")
message(sprintf("rank_purity joined for %d / %d genes", sum(!is.na(combined$summary$rank_purity)), n_before))

write.csv(combined$summary, combined_csv)
qsave(combined, combined_qs)
message("Updated ", combined_qs, " and ", combined_csv)
