setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)
library(qs)

objects_folder <- file.path(project_folder, "objects", "R")
ca_folder <- file.path(objects_folder, "clustassess")
mtd_folder <- file.path(project_folder, "metadata")
if (!dir.exists(ca_folder)) {
    dir.create(ca_folder, recursive = TRUE)
}

mtd_configs <- read.csv(file.path(mtd_folder, "configurations.csv"), header = TRUE, comment.char = "#")
nreps <- 100
neigh_seq <- seq(from = 5, to = 50, by = 5)
assay_name <- "RNA"
res_seq <- seq(from = 0.1, to = 2, by = 0.1)

for (i in seq_len(nrow(mtd_configs))) {
    id <- mtd_configs$id[i]
    so_path <- file.path(objects_folder, "seurat", paste0(id, "_norm.qs"))
    ca_path <- file.path(ca_folder, paste0(id, ".qs"))

    if (file.exists(ca_path) && file.size(ca_path) > 0) {
        next
    }

    print(id)

    so <- qread(file.path(objects_folder, "seurat", paste0(id, "_norm.qs")), nthreads = ncores)
    DefaultAssay(so) <- assay_name


    # extract the matrix and the features
    expr_matrix <- GetAssayData(
        object = so,
        layer = "scale.data"
    )
    features <- dimnames(so@assays[[assay_name]])[[1]]
    var_features <- VariableFeatures(so)
    max_ngenes <- min(length(var_features), 4500)

    most_abundant_genes <- rownames(expr_matrix)[order(Matrix::rowSums(expr_matrix), decreasing = TRUE)][seq_len(max_ngenes)]

    gene_list <- list(
        "Most_Abundant" = most_abundant_genes,
        "Highly_Variable" = var_features
    )

    steps_list <- list(
        "Most_Abundant" = seq(from = 500, by = 500, to = max_ngenes),
        "Highly_Variable" = seq(from = 500, by = 500, to = max_ngenes)
    )

    harmony_mtd <- mtd_configs$harmony_mtd[i]
    mtd_vals <- c()
    if (is.na(harmony_mtd) || harmony_mtd == "") {
        matrix_processing_function <- function(dt_mtx, actual_npcs = 30, ...) {
            actual_npcs <-
            min(actual_npcs, ncol(dt_mtx)%/%2)
            
            RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
            embedding <-
        stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x
            
            RhpcBLASctl::blas_set_num_threads(1)
        rownames(embedding) <- rownames(dt_mtx)
        
        colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))
            
            return(embedding)
        }
    } else {
        mtd_vals <- so@meta.data[[harmony_mtd]]
        matrix_processing_function <- function(dt_mtx, actual_npcs = 30) {
            actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)

            RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
            embedding <- stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x

            RhpcBLASctl::blas_set_num_threads(1)
            rownames(embedding) <- rownames(dt_mtx)
            colnames(embedding) <- paste0("PC_", seq_len(actual_npcs))

            embedding <- harmony::RunHarmony(embedding, mtd_vals, verbose = FALSE)

            return(embedding)
        }
    }

    expr_matrix <- expr_matrix[union(most_abundant_genes, var_features), ]
    rm(so)
    gc()

    RhpcBLASctl::blas_set_num_threads(1)
    my_cluster <- parallel::makeCluster(
        ncores,
        type = "PSOCK"
    )
    doParallel::registerDoParallel(cl = my_cluster)

    test_automm <- automatic_stability_assessment(
        expression_matrix = expr_matrix,
        n_repetitions = nreps,
        temp_file = "clustassess_temp.rds",
        save_temp = FALSE,
        n_neigh_sequence = neigh_seq,
        resolution_sequence = res_seq,
        features_sets = gene_list,
        steps = steps_list,
        n_top_configs = 3,
        matrix_processing = matrix_processing_function,
        umap_arguments = list(
            min_dist = 0.3,
            n_neighbors = 30,
            metric = "cosine"
        ),
        verbose = TRUE
    )

    parallel::stopCluster(cl = my_cluster)
    qsave(test_automm, file.path(ca_folder, paste0(id, ".qs")), nthreads = ncores)
}
