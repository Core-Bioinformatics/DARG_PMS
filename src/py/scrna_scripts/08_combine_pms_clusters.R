setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)
library(qs)

objects_folder <- file.path(project_folder, "objects", "R")
ca_folder <- file.path(objects_folder, "clustassess")
seurat_folder <- file.path(objects_folder, "seurat")
mtd_folder <- file.path(project_folder, "metadata")

dmso <- qread(file.path(ca_folder, "105_2g_DMSO_sub_0_4.qs"), nthreads = ncores)
abt <- qread(file.path(ca_folder, "105_2g_ABT_263.qs"), nthreads = ncores)

pms <- qread(file.path(seurat_folder, "pms_norm.qs"), nthreads = ncores)


dmso_clusters <- paste0("dmso_cluster_", get_clusters_from_clustassess_object(
    dmso,
    feature_type = "Highly_Variable",
    feature_size = 2500,
    clustering_method = "SLM",
    nclusters = 9
)[[1]]$partitions[[1]]$mb)
names(dmso_clusters) <- rownames(dmso[[2]][[1]]$pca)

abt_clusters <- paste0("abt_cluster_", get_clusters_from_clustassess_object(
    abt,
    feature_type = "Highly_Variable",
    feature_size = 2000,
    clustering_method = "SLM",
    nclusters = 11
)[[1]]$partitions[[1]]$mb)
names(abt_clusters) <- rownames(abt[[2]][[1]]$pca)

combined_clusters <- c(dmso_clusters, abt_clusters)
combined_clusters <- combined_clusters[colnames(pms)]

pms$"clusters_from_samples" <- combined_clusters

qsave(pms, file = file.path(seurat_folder, "pms_norm.qs"), nthreads = ncores)


healthy <- qread(file.path(seurat_folder, "healthy_norm.qs"), nthreads = ncores)
healthy_ca <- qread(file.path(ca_folder, "healthy.qs"), nthreads = ncores)

healthy_clusters <-  get_clusters_from_clustassess_object(
    healthy_ca,
    feature_type = "Highly_Variable",
    feature_size = 3000,
    clustering_method = "SLM",
    nclusters = 13
)[[1]]$partitions[[1]]$mb
healthy$stable_clusters <- healthy_clusters - 1

qsave(healthy, file = file.path(seurat_folder, "healthy_norm.qs"), nthreads = ncores)
