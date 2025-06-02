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

objects_folder <- file.path(project_folder, "objects")
so_folder <- file.path(objects_folder, "R", "seurat")
ca_folder <- file.path(objects_folder, "R", "clustassess")
marker_folder <- file.path(project_folder, "output", "gene_markers")
mat_folder <- list(
    "new" = file.path(objects_folder, "expr_matrix"),
    "old" = file.path(objects_folder, "original_expr_matrix")
)
peak_folder <- file.path(project_folder, "output", "peaks")
if (!dir.exists(marker_folder)) {
    dir.create(marker_folder, recursive = TRUE)
}

log1pdata.mean.fxn <- function(x, pseudocount.use = 1, base = 2) {
    return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x), base = base))
}

##### NEW SAMPELS #####
so <- qread(file.path(so_folder, "aggregate_norm.qs"), nthreads = ncores)
smp_pth <- file.path(marker_folder, "new_aggregate")
if (!dir.exists(smp_pth)) {
    dir.create(smp_pth, recursive = TRUE)
}
clusters <- so$aggregate_stable_clusters
new_ann <- read.table(file.path(peak_folder, "annotation_aggregate_peaks.txt"), header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
pc_genes <- new_ann %>% filter(Gene.Type == "protein-coding") %>% pull(Gene.Name) %>% unique


for (nup_comb in list(
    c(100, 1000),
    c(3000, 3000),
    c(50000, 50000)
)) {
    nup <- nup_comb[1]
    ndown <- nup_comb[2]
    
    expr_matrix <- qread(file.path(mat_folder[["new"]], paste0("peak_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")), nthreads = ncores)

    for (cl in unique(clusters)) {
        gr1 <- colnames(expr_matrix)[clusters == cl]
        gr2 <- colnames(expr_matrix)[clusters != cl]

        marker_df <- FindMarkers(
            object = expr_matrix,
            cells.1 = gr1,
            cells.2 = gr2,
            logfc.threshold = 0.5,
            test.use = "wilcox_limma",
            fc.result = FoldChange(
                object = expr_matrix,
                cells.1 = gr1,
                cells.2 = gr2,
                mean.fxn = log1pdata.mean.fxn,
                fc.name = "avg_log2FC",
                norm.method = "LogNormalize"
            )
        ) %>% filter(.data$p_val_adj < 0.05) %>% arrange(desc(.data$avg_log2FC))

        write.csv(
            marker_df,
            file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, ".csv")),
            row.names = TRUE
        )

        marker_df <- marker_df[intersect(rownames(marker_df), pc_genes), ]
        if (nrow(marker_df) == 0) {
            next
        }

        write.csv(
            marker_df,
            file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, "_pc_genes.csv")),
            row.names = TRUE
        )
    }
}


##### OLD SAMPLES #####
so_rna <- readRDS("metadata-integrated-2023-scRNAseq-harmony-exclude.rds")
so_atac <- readRDS("snATAC-seq-integrated-filter-Aug-2-2023.rds")
mtd_path <- read.csv("metadata-integrated-harmony-2023-matched-with-ATAC-cells.csv")

mtd_path <- mtd_path %>% filter(barcode %in% colnames(so_rna)) %>% filter(matched_cells %in% colnames(so_atac)) %>% filter(!(matched_cells %in% colnames(so_atac)[so_atac$seurat_clusters == 9]))

so_rna <- so_rna[, mtd_path$barcode]
so_atac <- so_atac[, mtd_path$matched_cells]

smp_pth <- file.path(marker_folder, "old_aggregate")
if (!dir.exists(smp_pth)) {
    dir.create(smp_pth, recursive = TRUE)
}
clusters <- so_atac$seurat_clusters
new_ann <- read.table(file.path(peak_folder, "annotation_old_aggregate_peaks.txt"), header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
pc_genes <- new_ann %>% filter(Gene.Type == "protein-coding") %>% pull(Gene.Name) %>% unique


for (nup_comb in list(
    c(100, 1000),
    c(3000, 3000),
    c(50000, 50000)
)) {
    nup <- nup_comb[1]
    ndown <- nup_comb[2]
    
    expr_matrix <- qread(file.path(mat_folder[["old"]], paste0("peak_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")), nthreads = ncores)

    for (cl in unique(clusters)) {
        gr1 <- colnames(expr_matrix)[clusters == cl]
        gr2 <- colnames(expr_matrix)[clusters != cl]

        marker_df <- FindMarkers(
            object = expr_matrix,
            cells.1 = gr1,
            cells.2 = gr2,
            logfc.threshold = 0.5,
            test.use = "wilcox_limma",
            fc.result = FoldChange(
                object = expr_matrix,
                cells.1 = gr1,
                cells.2 = gr2,
                mean.fxn = log1pdata.mean.fxn,
                fc.name = "avg_log2FC",
                norm.method = "LogNormalize"
            )
        ) %>% filter(.data$p_val_adj < 0.05) %>% arrange(desc(.data$avg_log2FC))

        write.csv(
            marker_df,
            file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, ".csv")),
            row.names = TRUE
        )

        marker_df <- marker_df[intersect(rownames(marker_df), pc_genes), ]
        if (nrow(marker_df) == 0) {
            next
        }

        write.csv(
            marker_df,
            file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, "_pc_genes.csv")),
            row.names = TRUE
        )
    }
}
