setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
source(file.path(general_scripts_folder, "RNA", "create_seurat_from_cr_h5.R"))
library(ggplot2)
library(Seurat)
library(qs)
library(tidyverse)

counts_folder <- file.path(project_folder, "preprocessing", "cellranger")
objects_folder <- file.path(project_folder, "objects", "R", "seurat")
if (!dir.exists(objects_folder)) {
    dir.create(objects_folder, recursive = TRUE)
}
metadata_folder <- file.path(project_folder, "metadata")

so <- create_so_default(
    h5_path = file.path(counts_folder, "105_2g_DMSO", "outs", "raw_feature_bc_matrix.h5"),
    sample_name = "105_2g_DMSO"
)
so$treatment <- "non_treated"
so$sample_type <- "PMS"
so$sample_name <- "105_2g_DMSO"

mt_regex <- "^MT-"
rp_regex <- "^RP[SL]\\d+"
mrp_regex <- "^MRP[SL]\\d+"
mt_genes <- grep(mt_regex, rownames(so), value = FALSE)
rp_genes <- grep(rp_regex, rownames(so), value = FALSE)
mrp_genes <- grep(mrp_regex, rownames(so), value = FALSE)

so$percent_mt <- PercentageFeatureSet(so, features = rownames(so)[mt_genes])
so$percent_rp <- PercentageFeatureSet(so, features = rownames(so)[c(rp_genes, mrp_genes)])
so$percent_mt[is.na(so$percent_mt)] <- 100
so$percent_rp[is.na(so$percent_rp)] <- 100
so <- subset(so, features = -c(mt_genes, rp_genes, mrp_genes))

so <- subset(so, nCount_RNA > 3e3 & nCount_RNA < 6e4 & nFeature_RNA > 1e3 & percent_mt < 50 & percent_rp < 50 & percent_mt < 10 & percent_mt > 2 & percent_rp < 40 & percent_rp > 15)
VlnPlot(so, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp"), pt.size = 0, log = TRUE, ncol = 4)
VlnPlot(so, features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp"), pt.size = 0, log = FALSE, ncol = 4)

n_median_genes <- floor(median(so$nFeature_RNA) / 500) * 500
cc_genes <- intersect(unlist(cc.genes.updated.2019), rownames(so))
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = n_median_genes + length(cc_genes), verbose = TRUE)
so <- ScaleData(so, features = rownames(so), verbose = TRUE)

var_features <- setdiff(VariableFeatures(so), cc_genes)[seq_len(n_median_genes)]
VariableFeatures(so) <- var_features
print(length(VariableFeatures(so)))
so <- CellCycleScoring(so, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

RhpcBLASctl::blas_set_num_threads(ncores)
so <- RunPCA(so, npcs = 30, approx = TRUE, verbose = TRUE)
so <- RunUMAP(so, reduction = "pca", dims = 1:30, verbose = TRUE, umap.method = "uwot")
qsave(so, file.path(objects_folder, paste0("105_2g_DMSO", "_norm.qs")), nthreads = ncores)
