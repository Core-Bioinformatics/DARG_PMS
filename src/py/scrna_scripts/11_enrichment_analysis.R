setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(future)
library(foreach)
library(dplyr)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)
library(gprofiler2)


objects_folder <- file.path(project_folder, "objects", "R")
so_folder <- file.path(objects_folder, "seurat")

markers_folder <- file.path(project_folder, "output", "RNA", "markers")
enrich_folder <- file.path(project_folder, "output", "RNA", "enrichment")

dts <- list(
    "healthy" = c("stable_clusters"),
    "pms" = c("clusters_from_samples")
)

for (dts_name in names(dts)) {
    so <- qread(file.path(so_folder, paste0(dts_name, "_norm.qs")), nthreads = ncores)
    expr_matrix <- GetAssayData(so, slot = "data")
    genes <- unique(rownames(expr_matrix))
    genes <- genes[apply(expr_matrix, 1, function(x) min(x) != max(x))]

    gene_conv <- select(EnsDb.Hsapiens.v86, keys = genes, columns = c("GENENAME", "GENEBIOTYPE"), keytype = "GENENAME")
    pc_genes <- unique(gene_conv[gene_conv$GENEBIOTYPE == "protein_coding", "GENENAME"])

    for (mtd_name in dts[[dts_name]]) {
        tmp_input_folder <- file.path(markers_folder, dts_name, mtd_name)
        tmp_output_folder <- file.path(enrich_folder, dts_name, mtd_name)
        if (!dir.exists(tmp_output_folder)) {
            dir.create(tmp_output_folder, recursive = TRUE)
        }

        for (f in list.files(tmp_input_folder)) {
            if (!startsWith(f, "markers_")) {
                next
            }
            
            if (!endsWith(f, ".csv")) {
                next
            }

            print(paste("Perform enrichment on", f))
            new_name <- paste0("enrichment_", strsplit(f, "markers_")[[1]][2])
            marker_genes <- read.csv(file.path(tmp_input_folder, f), header = TRUE, row.names = 1) %>% dplyr::filter(avg_log2FC > 0) %>% rownames()

            gprof_res <- gost(
                query = marker_genes,
                organism = "hsapiens",
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF"),
                evcodes = TRUE,
                domain_scope = "custom",
                custom_bg = genes
            )

            if (!is.null(gprof_res)) {
                gprof_res <- gprof_res$result
                gprof_res$parents <- sapply(gprof_res$parents, toString)
                rownames(gprof_res) <- gprof_res$term_id

                write.csv(gprof_res %>% arrange(p_value), file.path(tmp_output_folder, new_name), row.names = TRUE)
            }

            new_name <- paste0("PC_enrichment_", strsplit(f, "markers_")[[1]][2])
            gprof_res <- gost(
                query = intersect(marker_genes, pc_genes),
                organism = "hsapiens",
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF"),
                evcodes = TRUE,
                domain_scope = "custom",
                custom_bg = pc_genes
            )

            if (!is.null(gprof_res)) {
                gprof_res <- gprof_res$result
                gprof_res$parents <- sapply(gprof_res$parents, toString)
                rownames(gprof_res) <- gprof_res$term_id

                write.csv(gprof_res %>% arrange(p_value), file.path(tmp_output_folder, new_name), row.names = TRUE)
            }
        }
    }
}
