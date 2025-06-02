setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(future)
library(foreach)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(gprofiler2)

objects_folder <- file.path(project_folder, "objects", "R")
so_folder <- file.path(objects_folder, "seurat")

mtd_folder <- file.path(project_folder, "metadata")
output_folder <- file.path(project_folder, "output", "RNA", "markers")
ca_folder <- file.path(project_folder, "shiny_apps", "clustassess")

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

min_pct <- 0.05
logfc_thresh <- 0.5
pval_thresh <- 0.05

plan(multicore, workers = ncores)

#### HEALTHY
so <- qread(file.path(so_folder, "healthy_norm.qs"), nthreads = ncores)

target_mtd <- "stable_clusters"
temp_output_folder <- file.path(output_folder, "healthy", target_mtd)
if (!dir.exists(temp_output_folder)) {
  dir.create(temp_output_folder, recursive = TRUE)
}

for (cl_val in unique(so[[target_mtd]][,1])) {
    group_mask <- rep("group2", ncol(so))
    group_mask[so[[target_mtd]][,1] == cl_val] <- "group1"

    if (sum(group_mask == "group1") < 5) {
        next
    }

    if (sum(group_mask == "group2") < 5) {
        next
    }
    print(paste("Identifying markers for cluster", cl_val))

    Idents(so) <- factor(group_mask)

    markers_df <- FindMarkers(
        so,
        ident.1 = "group1",
        ident.2 = "group2",
        logfc.threshold = logfc_thresh,
        test.use = "wilcox_limma",
        min.pct = min_pct,
        only.pos = FALSE
    ) %>% dplyr::filter(
        .data$p_val_adj < pval_thresh
    ) %>% arrange(
        desc(.data$avg_log2FC)
    )

    write.csv(
        markers_df,
        file.path(temp_output_folder, paste0("markers_", cl_val, ".csv")),
        row.names = TRUE
    )

    gene_names <- rownames(markers_df)
    gene_conv <- select(EnsDb.Hsapiens.v86, keys = gene_names, columns = c("GENENAME", "GENEBIOTYPE"), keytype = "GENENAME")
    pc_genes <- intersect(unique(gene_conv[gene_conv$GENEBIOTYPE == "protein_coding", "GENENAME"]), gene_names)

    write.csv(
        markers_df[pc_genes, ],
        file.path(temp_output_folder, paste0("PC_markers_", cl_val, ".csv")),
        row.names = TRUE
    )
}

#### PMS
so <- qread(file.path(so_folder, "pms_norm.qs"), nthreads = ncores)
target_mtd <- "clusters_from_samples"
temp_output_folder <- file.path(output_folder, "pms", target_mtd)
if (!dir.exists(temp_output_folder)) {
  dir.create(temp_output_folder, recursive = TRUE)
}

for (cl_val in unique(so[[target_mtd]][,1])) {
    group_mask <- rep("group2", ncol(so))
    group_mask[so[[target_mtd]][,1] == cl_val] <- "group1"

    if (sum(group_mask == "group1") < 5) {
        next
    }

    if (sum(group_mask == "group2") < 5) {
        next
    }
    print(paste("Identifying markers for cluster", cl_val))

    Idents(so) <- factor(group_mask)

    markers_df <- FindMarkers(
        so,
        ident.1 = "group1",
        ident.2 = "group2",
        logfc.threshold = logfc_thresh,
        test.use = "wilcox_limma",
        min.pct = min_pct,
        only.pos = FALSE
    ) %>% dplyr::filter(
        .data$p_val_adj < pval_thresh
    ) %>% arrange(
        desc(.data$avg_log2FC)
    )

    write.csv(
        markers_df,
        file.path(temp_output_folder, paste0("markers_", cl_val, ".csv")),
        row.names = TRUE
    )

    gene_names <- rownames(markers_df)
    gene_conv <- select(EnsDb.Hsapiens.v86, keys = gene_names, columns = c("GENENAME", "GENEBIOTYPE"), keytype = "GENENAME")
    pc_genes <- intersect(unique(gene_conv[gene_conv$GENEBIOTYPE == "protein_coding", "GENENAME"]), gene_names)

    write.csv(
        markers_df[pc_genes, ],
        file.path(temp_output_folder, paste0("PC_markers_", cl_val, ".csv")),
        row.names = TRUE
    )
}

plan(sequential)
