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

objects_folder <- file.path(project_folder, "objects", "R")
so_folder <- file.path(objects_folder, "seurat")
ca_folder <- file.path(objects_folder, "clustassess")
peaks_folder <- file.path(project_folder, "output", "peaks") 
peak_ann <- read.table(file.path(peaks_folder, "annotation_aggregate_peaks.txt"), header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
peak_ann <- peak_ann %>% filter(!is.na(Distance.to.TSS))
rownames(peak_ann) <- paste0(peak_ann$Chr, "-", peak_ann$Start - 1, "-", peak_ann$End)
promoter_peaks <- peak_ann %>% filter(startsWith(Annotation, "promoter-TSS"))

old_peak_ann <- read.table(file.path(peaks_folder, "annotation_old_aggregate_peaks.txt"), header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
old_peak_ann <- old_peak_ann %>% filter(!is.na(Distance.to.TSS))
rownames(old_peak_ann) <- paste0(old_peak_ann$Chr, "-", old_peak_ann$Start - 1, "-", old_peak_ann$End)

if (!file.exists(file.path(peaks_folder, "annotation_aggregate_promoter_peaks.txt"))) {
    write.table(promoter_peaks, file.path(peaks_folder, "annotation_aggregate_promoter_peaks.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

peaks_folder <- file.path(peaks_folder, "markers")
# markers for new aggregate samples
current_pth <- file.path(peaks_folder, "new_aggregate")
if (!dir.exists(current_pth)) {
    dir.create(current_pth, recursive = TRUE)
}

so <- qread(file.path(so_folder, "aggregate_norm.qs"), nthreads = ncores)
Idents(so) <- "aggregate_stable_clusters"
for (i in unique(so$aggregate_stable_clusters)) {
    print(i)
    output_path <- file.path(current_pth, paste0("markers_", i, ".csv"))
    Idents(so) <- ifelse(so$aggregate_stable_clusters == i, "cluster", "rest")
    markers <- FindMarkers(
        so,
        ident.1 = "cluster",
        ident.2 = "rest",
        test.use = "wilcox_limma",
        min.pct = 0.1,
        logfc.threshold = 0.5,
    ) %>% filter(.data$p_val_adj < 0.05) %>% arrange(desc(.data$avg_log2FC))

    if (nrow(markers) == 0) {
        next
    }

    tmp_df <- peak_ann[rownames(markers), ]
    tmp_df[, colnames(markers)] <- markers

    write.table(
        tmp_df,
        output_path,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
    )

    # match the peaks with promoter_peaks
    tmp_df <- tmp_df %>% filter(startsWith(Annotation, "promoter-TSS"))

    if (nrow(tmp_df) == 0) {
        next
    }

    output_path <- file.path(current_pth, paste0("markers_promoter_", i, ".csv"))
    write.table(
        tmp_df,
        output_path,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
    )
}

# markers for the old aggregate
current_pth <- file.path(peaks_folder, "old_aggregate")
if (!dir.exists(current_pth)) {
    dir.create(current_pth, recursive = TRUE)
}

so_rna <- readRDS("/servers/iss-corescratch/lp488/bongsoo/scrna/metadata-integrated-2023-scRNAseq-harmony-exclude.rds")
so_atac <- readRDS("/servers/iss-corescratch/lp488/bongsoo/scrna/snATAC-seq-integrated-filter-Aug-2-2023.rds")
mtd_path <- read.csv("/servers/iss-corescratch/lp488/bongsoo/scrna/data/atac/atac/metadata-integrated-harmony-2023-matched-with-ATAC-cells.csv")

mtd_path <- mtd_path %>% filter(barcode %in% colnames(so_rna)) %>% filter(matched_cells %in% colnames(so_atac)) %>% filter(!(matched_cells %in% colnames(so_atac)[so_atac$seurat_clusters == 9]))

so_rna <- so_rna[, mtd_path$barcode]
so_atac <- so_atac[, mtd_path$matched_cells]

for (i in unique(so_atac$seurat_clusters)) {
    print(i)
    output_path <- file.path(current_pth, paste0("markers_", i, ".csv"))
    Idents(so_atac) <- ifelse(so_atac$seurat_clusters == i, "cluster", "rest")
    markers <- FindMarkers(
        so_atac,
        ident.1 = "cluster",
        ident.2 = "rest",
        test.use = "wilcox_limma",
        min.pct = 0.1,
        logfc.threshold = 0.5,
    ) %>% filter(.data$p_val_adj < 0.05) %>% arrange(desc(.data$avg_log2FC))

    if (nrow(markers) == 0) {
        next
    }

    tmp_df <- old_peak_ann[rownames(markers), ]
    tmp_df[, colnames(markers)] <- markers

    write.table(
        tmp_df,
        output_path,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
    )

    # match the peaks with promoter_peaks
    tmp_df <- tmp_df %>% filter(startsWith(Annotation, "promoter-TSS"))

    if (nrow(tmp_df) == 0) {
        next
    }

    output_path <- file.path(current_pth, paste0("markers_promoter_", i, ".csv"))
    write.table(
        tmp_df,
        output_path,
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
    )
}
