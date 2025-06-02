setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(foreach)
library(Signac)
library(ClustAssess)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
qc_folder <- file.path(project_folder, "preprocessing", "qc", "ATAC")
ca_folder <- file.path(project_folder, "objects", "R", "clustassess")
app_folder <- file.path(project_folder, "shiny_apps", "clustassess")
mt_folder <- file.path(project_folder, "objects", "expr_matrix")

so <- qread(file.path(objects_folder, "aggregate_norm.qs"), nthreads = ncores)
ca <- qread(file.path(ca_folder, "aggregate.qs"), nthreads = ncores)

nup <- 50000
ndown <- 50000

for (base_inf in c("fragment", "peak")[2]) {
    expr_matrix <- qread(file.path(mt_folder, paste0(base_inf, "_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")), nthreads = ncores)

    write_shiny_app(
        object = expr_matrix,
        metadata = so@meta.data,
        clustassess_object = ca,
        shiny_app_title = paste0("ATAC aggregate - RNA expression based on ", base_inf, " up ", nup, " down ", ndown),
        project_folder = file.path(app_folder, paste0("aggregate_rna_", base_inf, "_up_", nup, "_down_", ndown)),
        prompt_feature_choice = FALSE
    )
}
