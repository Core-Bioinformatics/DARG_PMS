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
library(ComplexHeatmap)


objects_folder <- file.path(project_folder, "objects", "R")
so_folder <- file.path(objects_folder, "seurat")

input_folders <- list(
    markers = file.path(project_folder, "output", "RNA", "markers"),
    enrichment = file.path(project_folder, "output", "RNA", "enrichment")
)

dts <- list(
    "healthy" = c("stable_clusters"),
    "pms" = c("clusters_from_samples")
)

n_top <- 100

for (dts_name in names(dts)) {
    for (mtd_name in dts[[dts_name]]) {
        for (comb in c("markers", "enrichment")) {
            tmp_input_folder <- file.path(input_folders[[comb]], dts_name, mtd_name)
            file_list <- list.files(tmp_input_folder)
            file_list <- file_list[endsWith(file_list, ".csv")]
            term_list <- list()
            for (pc in c("", "PC_")) {
                temp_fl <- file_list[startsWith(file_list, paste0(pc, comb))]
                temp_output <- file.path(tmp_input_folder, paste0(pc, comb, "_top_", n_top, "_jsi.pdf"))

                for (temp_f in temp_fl) {
                    cl_name <- strsplit(temp_f, paste0(pc, comb, "_"))[[1]][2]
                    cl_name <- strsplit(cl_name, ".csv")[[1]][1]

                    terms <- read.csv(file.path(tmp_input_folder, temp_f), header = TRUE, row.names = 1)
                    terms <- rownames(terms)[seq_len(n_top)]
                    term_list[[cl_name]] <- terms
                }

                jsi_matrix <- matrix(0, nrow = length(term_list), ncol = length(term_list))
                rownames(jsi_matrix) <- names(term_list)
                colnames(jsi_matrix) <- names(term_list)



                for (i in seq_along(term_list)) {
                    for (j in seq_along(term_list)) {
                        if (i == j) {
                            jsi_matrix[i, j] <- NA
                        } else {
                            jsi_matrix[i, j] <- length(intersect(term_list[[i]], term_list[[j]])) / length(union(term_list[[i]], term_list[[j]]))
                        }
                    }
                }

                if (dts_name == "pms") {
                    print(temp_fl)
                    jsi_matrix <- jsi_matrix[startsWith(rownames(jsi_matrix), "dmso"), startsWith(colnames(jsi_matrix), "abt")]
                    jsi_matrix <- t(jsi_matrix)
                }

                col_fun <- circlize::colorRamp2(c(0, 0.5, 1), RColorBrewer::brewer.pal(9, "Greens")[c(1, 5, 9)])
                htmp <- Heatmap(
                    jsi_matrix,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", jsi_matrix[i, j]), x, y, gp = gpar(fontsize = 13))

                    },
                    col = col_fun,
                    column_title = paste0("JSI - Top ", n_top, " ", pc, comb, " terms")
                )

                pdf(temp_output, width = nrow(jsi_matrix) * 0.9, height = ncol(jsi_matrix) * 0.9)
                print(htmp)
                dev.off()
            }
        }
    }

}
