# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(".")
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(foreach)
library(Signac)
library(ClustAssess)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
qc_folder <- file.path(project_folder, "preprocessing", "qc", "ATAC")
ca_folder <- file.path(project_folder, "objects", "R", "clustassess")
if (!dir.exists(ca_folder)) {
    dir.create(ca_folder, recursive = TRUE)
}
shiny_folder <- file.path(project_folder, "shiny_apps", "clustassess")
if (!dir.exists(shiny_folder)) {
    dir.create(shiny_folder, recursive = TRUE)
}

so <- qs::qread(file.path(objects_folder, "aggregate_norm.qs"), nthreads = ncores)

my_cluster <- parallel::makeCluster(
    16,
    type = "PSOCK"
)
doParallel::registerDoParallel(cl = my_cluster)

var_peaks <- lapply(paste0("q", c(80, 85, 90, 95, 96, 97, 98, 99)), function(x) {
    temp_so <- FindTopFeatures(so, min.cutoff = x)
    VariableFeatures(temp_so)
})
ordered_var_peaks <- c()
for (i in seq(from = length(var_peaks), to = 1, by = -1)) {
    ordered_var_peaks <- c(ordered_var_peaks, setdiff(var_peaks[[i]], ordered_var_peaks))
}

gene_list <- list(
    "Highly_Variable_peaks" = ordered_var_peaks
)
steps_list <- list(
    "Highly_Variable_peaks" = rev(sapply(var_peaks, length))
)

stype <- so@meta.data[, "sample_type", drop = FALSE]
# convert stype to vector
print(stype)
matrix_proc_function <- function(dt_matrix, actual_npcs = 30) {
    dt_matrix <- t(dt_matrix)
    actual_npcs <- min(actual_npcs, ncol(dt_matrix) %/% 2)

    RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
    embedding <- Signac::RunSVD(dt_matrix, n = actual_npcs, verbose = FALSE)@cell.embeddings

    rownames(embedding) <- colnames(dt_matrix)
    if (is.null(colnames(embedding))) {
        colnames(embedding) <- paste0("LSI_", seq_len(actual_npcs))
    }
    print(dim(embedding))

    embedding <- harmony::RunHarmony(embedding, stype, vars_use = "sample_type", verbose = FALSE)

    return(embedding)
}

expr_matrix <- as.matrix(GetAssayData(so, slot = "data")[ordered_var_peaks, ])
rm(so)
gc()
test_automm <- automatic_stability_assessment(
      expression_matrix = expr_matrix,
      n_repetitions = 100,
      temp_file = "clustassess_temp.rds",
      save_temp = FALSE,
      n_neigh_sequence = seq(from = 5, to = 50, by = 5),
      resolution_sequence = seq(from = 0.1, to = 2, by = 0.1),
      features_sets = gene_list,
      steps = steps_list,
      n_top_configs = 4,
      umap_arguments = list(
          min_dist = 0.3,
          n_neighbors = 30,
          metric = "cosine"
      ),
      matrix_processing = matrix_proc_function,
      verbose = TRUE
)

qs::qsave(test_automm, file.path(ca_folder, "aggregate.qs"), nthreads = ncores)
parallel::stopCluster(cl = my_cluster)
foreach::registerDoSEQ()
gc()
so <- qread(file.path(objects_folder, "aggregate_norm.qs"), nthreads = ncores)
write_shiny_app(
    object = expr_matrix,
    metadata = so@meta.data,
    clustassess_object = test_automm,
    shiny_app_title = "ATAC aggregate - Harmony integration",
    project_folder = file.path(shiny_folder, "aggregate"),
    prompt_feature_choice = FALSE
)



