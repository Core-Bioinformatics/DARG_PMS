setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(future)
library(foreach)
library(dplyr)
library(gprofiler2)
library(foreach)
library(future.apply)

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
enrichment_folder <- file.path(project_folder, "output", "enrichment")
if (!dir.exists(enrichment_folder)) {
    dir.create(enrichment_folder, recursive = TRUE)
}

##### NEW SAMPELS #####
so <- qread(file.path(so_folder, "aggregate_norm.qs"), nthreads = ncores)
smp_pth <- file.path(marker_folder, "new_aggregate")
en_pth <- file.path(enrichment_folder, "new_aggregate")
if (!dir.exists(en_pth)) {
    dir.create(en_pth, recursive = TRUE)
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
    background_genes <- rownames(expr_matrix)
    background_pc_genes <- intersect(background_genes, pc_genes)

    future_lapply(unique(clusters), function(cl) {
        mk_path <- file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, ".csv"))
        mk_genes <- read.csv(mk_path, header = TRUE, row.names = 1)
        mk_genes <- mk_genes %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% rownames

        enriched_terms <- gost(
            query = mk_genes,
            organism = "hsapiens",
            sources = c("GO", "KEGG", "REAC"),
            evcodes = TRUE,
            user_threshold = 1
        )

        if (!is.null(enriched_terms)) {
            enriched_terms <- enriched_terms$result
            enriched_terms$parents <- sapply(enriched_terms$parents, toString)
            rownames(enriched_terms) <- enriched_terms$term_id

            enrich_path <- file.path(en_pth, paste0("enrichment_", cl, "_up_", nup, "_down_", ndown, ".csv"))
            write.csv(enriched_terms %>% arrange(p_value), enrich_path, row.names = TRUE)
        }

        mk_path <- file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, "_pc_genes.csv"))
        mk_genes <- read.csv(mk_path, header = TRUE, row.names = 1)
        mk_genes <- mk_genes %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% rownames

        enriched_terms <- gost(
            query = intersect(mk_genes, background_pc_genes),
            organism = "hsapiens",
            sources = c("GO", "KEGG", "REAC"),
            evcodes = TRUE,
            user_threshold = 1
        )

        if (!is.null(enriched_terms)) {
            enriched_terms <- enriched_terms$result
            enriched_terms$parents <- sapply(enriched_terms$parents, toString)
            rownames(enriched_terms) <- enriched_terms$term_id

            enrich_path <- file.path(en_pth, paste0("enrichment_", cl, "_up_", nup, "_down_", ndown, "_pc_genes.csv"))
            write.csv(enriched_terms %>% arrange(p_value), enrich_path, row.names = TRUE)
        }
    })
}

##### OLD SAMPELS #####
so_rna <- readRDS("metadata-integrated-2023-scRNAseq-harmony-exclude.rds")
so_atac <- readRDS("snATAC-seq-integrated-filter-Aug-2-2023.rds")
mtd_path <- read.csv("metadata-integrated-harmony-2023-matched-with-ATAC-cells.csv")

mtd_path <- mtd_path %>% filter(barcode %in% colnames(so_rna)) %>% filter(matched_cells %in% colnames(so_atac)) %>% filter(!(matched_cells %in% colnames(so_atac)[so_atac$seurat_clusters == 9]))

so_rna <- so_rna[, mtd_path$barcode]
so_atac <- so_atac[, mtd_path$matched_cells]

smp_pth <- file.path(marker_folder, "old_aggregate")
en_pth <- file.path(enrichment_folder, "old_aggregate")
if (!dir.exists(en_pth)) {
    dir.create(en_pth, recursive = TRUE)
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

    background_genes <- rownames(expr_matrix)
    background_pc_genes <- intersect(background_genes, pc_genes)

    future_lapply(unique(clusters), function(cl) {
        mk_path <- file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, ".csv"))
        mk_genes <- read.csv(mk_path, header = TRUE, row.names = 1)
        mk_genes <- mk_genes %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% rownames

        enriched_terms <- gost(
            query = mk_genes,
            organism = "hsapiens",
            sources = c("GO", "KEGG", "REAC"),
            evcodes = TRUE,
            user_threshold = 1
        )

        if (!is.null(enriched_terms)) {
            enriched_terms <- enriched_terms$result
            enriched_terms$parents <- sapply(enriched_terms$parents, toString)
            rownames(enriched_terms) <- enriched_terms$term_id

            enrich_path <- file.path(en_pth, paste0("enrichment_", cl, "_up_", nup, "_down_", ndown, ".csv"))
            write.csv(enriched_terms %>% arrange(p_value), enrich_path, row.names = TRUE)
        }

        mk_path <- file.path(smp_pth, paste0("markers_", cl, "_up_", nup, "_down_", ndown, "_pc_genes.csv"))
        mk_genes <- read.csv(mk_path, header = TRUE, row.names = 1)
        mk_genes <- mk_genes %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>% rownames

        enriched_terms <- gost(
            query = intersect(mk_genes, background_pc_genes),
            organism = "hsapiens",
            sources = c("GO", "KEGG", "REAC"),
            evcodes = TRUE,
            user_threshold = 1
        )

        if (!is.null(enriched_terms)) {
            enriched_terms <- enriched_terms$result
            enriched_terms$parents <- sapply(enriched_terms$parents, toString)
            rownames(enriched_terms) <- enriched_terms$term_id

            enrich_path <- file.path(en_pth, paste0("enrichment_", cl, "_up_", nup, "_down_", ndown, "_pc_genes.csv"))
            write.csv(enriched_terms %>% arrange(p_value), enrich_path, row.names = TRUE)
        }
    })
}
