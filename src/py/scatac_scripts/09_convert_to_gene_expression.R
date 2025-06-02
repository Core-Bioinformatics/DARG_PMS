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
library(future)
library(foreach)
library(future.apply)

seurat_path <- file.path(project_folder, "objects", "R", "seurat")
output_path <- file.path(project_folder, "output")
mt_path <- file.path(project_folder, "objects", "expr_matrix")
if (!dir.exists(mt_path)) {
    dir.create(mt_path, recursive = TRUE)
}

so <- qread(file.path(seurat_path, "aggregate_norm.qs"), nthreads = ncores)
# PEAK BASED
ann_peaks <- read.table(file.path(output_path, "peaks", "annotation_aggregate_peaks.txt"), header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
rownames(ann_peaks) <- paste0(ann_peaks$Chr, "-", ann_peaks$Start - 1, "-", ann_peaks$End)
nup <- 100
ndown <- 1000
promoter_peaks <- ann_peaks %>% filter(Distance.to.TSS <= nup, Distance.to.TSS >= -ndown)

subso <- so[rownames(promoter_peaks), ]
unique_genes <- unique(promoter_peaks$Gene.Name)
unique_genes <- setdiff(unique_genes, c("NA", ""))
gene_peak_mapping <- list()
for (pk in rownames(promoter_peaks)) {
    gname <- promoter_peaks[pk, "Gene.Name"]
    if (gname %in% names(gene_peak_mapping)) {
        gene_peak_mapping[[gname]] <- c(gene_peak_mapping[[gname]], pk)
    } else {
        gene_peak_mapping[[gname]] <- c(pk)
    }
}

count_matrix <- GetAssayData(subso, slot = "counts")

plan(multicore, workers = ncores)
options(future.globals.maxSize = 50 * 1024^3)
rows <- future_lapply(unique_genes, function(gname) {
    colSums(count_matrix[gene_peak_mapping[[gname]], , drop = FALSE])
})
gene_expr <- do.call(rbind, rows)

rownames(gene_expr) <- unique_genes
colnames(gene_expr) <- colnames(subso)

gene_expr <- NormalizeData(gene_expr, normalization.method = "LogNormalize", scale.factor = 10000)
qsave(gene_expr, file.path(mt_path, paste0("peak_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")),nthreads = ncores)

# FRAGMENT BASED
gene_expr <- GeneActivity(so, extend.upstream = nup, extend.downstream = ndown)
gene_expr <- NormalizeData(gene_expr, normalization.method = "LogNormalize", scale.factor = 10000)
qsave(gene_expr, file.path(mt_path, paste0("fragment_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")),nthreads = ncores)


##### ORIGINAL OBJECT #####
mt_path <- file.path(project_folder, "objects", "original_expr_matrix")
if (!dir.exists(mt_path)) {
    dir.create(mt_path, recursive = TRUE)
}

so_rna <- readRDS("metadata-integrated-2023-scRNAseq-harmony-exclude.rds")
so_atac <- readRDS("snATAC-seq-integrated-filter-Aug-2-2023.rds")
mtd_path <- read.csv("metadata-integrated-harmony-2023-matched-with-ATAC-cells.csv")

mtd_path <- mtd_path %>% filter(barcode %in% colnames(so_rna)) %>% filter(matched_cells %in% colnames(so_atac)) %>% filter(!(matched_cells %in% colnames(so_atac)[so_atac$seurat_clusters == 9]))

so_rna <- so_rna[, mtd_path$barcode]
so_atac <- so_atac[, mtd_path$matched_cells]

# PEAK BASED
ann_peaks <- read.table(file.path(output_path, "peaks", "annotation_old_aggregate_peaks.txt"), header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
rownames(ann_peaks) <- paste0(ann_peaks$Chr, "-", ann_peaks$Start - 1, "-", ann_peaks$End)
promoter_peaks <- ann_peaks %>% filter(Distance.to.TSS <= nup, Distance.to.TSS >= -ndown)

subso <- so_atac[rownames(promoter_peaks), ]
unique_genes <- unique(promoter_peaks$Gene.Name)
unique_genes <- setdiff(unique_genes, c("NA", ""))
gene_peak_mapping <- list()
for (pk in rownames(promoter_peaks)) {
    gname <- promoter_peaks[pk, "Gene.Name"]
    if (gname %in% names(gene_peak_mapping)) {
        gene_peak_mapping[[gname]] <- c(gene_peak_mapping[[gname]], pk)
    } else {
        gene_peak_mapping[[gname]] <- c(pk)
    }
}

count_matrix <- GetAssayData(subso, slot = "counts")

plan(multicore, workers = ncores)
options(future.globals.maxSize = 50 * 1024^3)
rows <- future_lapply(unique_genes, function(gname) {
    colSums(count_matrix[gene_peak_mapping[[gname]], , drop = FALSE])
})
gene_expr <- do.call(rbind, rows)

rownames(gene_expr) <- unique_genes
colnames(gene_expr) <- colnames(subso)

gene_expr <- NormalizeData(gene_expr, normalization.method = "LogNormalize", scale.factor = 10000)
qsave(gene_expr, file.path(mt_path, paste0("peak_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")),nthreads = ncores)

# FRAGMENT BASED
gene_expr <- GeneActivity(so_atac, extend.upstream = nup, extend.downstream = ndown)
gene_expr <- NormalizeData(gene_expr, normalization.method = "LogNormalize", scale.factor = 10000)
qsave(gene_expr, file.path(mt_path, paste0("fragment_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")),nthreads = ncores)


