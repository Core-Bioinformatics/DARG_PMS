setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(Signac)
library(ClustAssess)
library(qs)
library(rhdf5)
library(dplyr)
library(future)
library(doParallel)
library(ComplexHeatmap)
library(ggplot2)

seurat_path <- file.path(project_folder, "objects", "R", "seurat")
so <- qread(file.path(seurat_path, "aggregate_norm.qs"), nthreads = ncores)

so$stable_clusters_9
DimPlot(so, group.by = "stable_clusters_9")

ctg_table <- table(so$sample_type, so$stable_clusters_9)
ctg_table <- ctg_table / rowSums(ctg_table) * 100

ctg_table
# combine 1,2 ,4
ctg_table <- cbind(ctg_table, rowSums(ctg_table[, c("1", "2", "4")]))
colnames(ctg_table)[ncol(ctg_table)] <- "1_2_4"
ctg_table <- ctg_table[, setdiff(colnames(ctg_table), c("1", "2", "4"))]

ctg_table

chisq.test(ctg_table[, c("1_2_4")])
chisq.test(ctg_table[, "5"])


ctg_table <- table(so$sample_type, so$stable_clusters_9)
ctg_table <- cbind(ctg_table, rowSums(ctg_table[, c("1", "2", "4")]))
colnames(ctg_table)[ncol(ctg_table)] <- "1_2_4"
ctg_table <- ctg_table[, setdiff(colnames(ctg_table), c("1", "2", "4"))]


# for (target_gr in c("1_2_4", "5", "6", "7", "3", "9", "8")) {
for (target_gr in colnames(ctg_table)) {

    compl <- cbind(ctg_table[, target_gr], rowSums(ctg_table[, setdiff(colnames(ctg_table), target_gr)]))
    compl

    print(paste(target_gr, fisher.test(compl)$p.value, sep = ","))
}

