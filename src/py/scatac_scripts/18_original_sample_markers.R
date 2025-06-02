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

# 50 GB
options(future.globals.maxSize = 50 * 1024^3)
plan(multicore, workers = ncores)

so_rna <- readRDS("metadata-integrated-2023-scRNAseq-harmony-exclude.rds")
so_atac <- readRDS("snATAC-seq-integrated-filter-Aug-2-2023.rds")
mtd_path <- read.csv("metadata-integrated-harmony-2023-matched-with-ATAC-cells.csv")

mtd_path <- mtd_path %>% filter(barcode %in% colnames(so_rna)) %>% filter(matched_cells %in% colnames(so_atac)) %>% filter(!(matched_cells %in% colnames(so_atac)[so_atac$seurat_clusters == 9]))

so_rna <- so_rna[, mtd_path$barcode]
so_atac <- so_atac[, mtd_path$matched_cells]
colnames(so_atac) <- mtd_path$barcode[match(colnames(so_atac), mtd_path$matched_cells)]
so_atac <- so_atac[, colnames(so_rna)]

so_atac$rna_clusters <- so_rna$seurat_clusters

rna_cl <- levels(so_rna$seurat_clusters)
marker_folder <- file.path(project_folder, "output", "gene_markers", "RNA_old_aggregate")
if (!dir.exists(marker_folder)) {
    dir.create(marker_folder, recursive = TRUE)
}
for (cl in rna_cl) {
    Idents(so_rna) <- ifelse(so_rna$seurat_clusters == cl, "cluster", "rest")
    output_path <- file.path(marker_folder, paste0("markers_", cl, ".csv"))
    markers <- FindMarkers(
        so_rna,
        ident.1 = "cluster",
        ident.2 = "rest",
        test.use = "wilcox_limma",
        min.pct = 0.1,
        logfc.threshold = 0.25,
        only.pos = TRUE
    ) %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))

    if (nrow(markers) == 0) {
        next
    }

    write.csv(
        markers,
        output_path,
        row.names = TRUE,
        quote = FALSE
    )
}
