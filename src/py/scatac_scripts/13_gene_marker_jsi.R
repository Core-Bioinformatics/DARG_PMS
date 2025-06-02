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

marker_path <- file.path(project_folder, "output", "gene_markers")
enrich_path <- file.path(project_folder, "output", "enrichment")
name1 <- "new_aggregate"
name2 <- "new_aggregate"

figure_path <- file.path(project_folder, "output", "figures", "new_sample_comparison")
if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}

clusters1 <- list.files(file.path(marker_path, name1))
clusters1 <- sapply(clusters1, function(x) {
    cl1_name <- strsplit(basename(x), "_up")[[1]][1]
    cl1_name <- strsplit(cl1_name, "_")[[1]][2]
    return(cl1_name)
})
names(clusters1) <- NULL
clusters1 <- unique(clusters1)

clusters2 <- list.files(file.path(marker_path, name2))
clusters2 <- sapply(clusters2, function(x) {
    cl2_name <- strsplit(basename(x), "_up")[[1]][1]
    cl2_name <- strsplit(cl2_name, "_")[[1]][2]
    return(cl2_name)
})
names(clusters2) <- NULL
clusters2 <- unique(clusters2)

for (pc_suf in c("", "_pc_genes")) {
    for (nup_comb in list(
        c(100, 1000),
        c(3000, 3000),
        c(50000, 50000)
    )) {
        nup <- nup_comb[1]
        ndown <- nup_comb[2]
        marker_list1 <- lapply(clusters1, function(i) {
            f <- file.path(marker_path, name1, paste0("markers_", i, "_up_", nup, "_down_", ndown, pc_suf, ".csv"))
            m <- read.csv(f, header = TRUE, row.names = 1, sep = ",", comment.char = "#") %>% arrange(desc(avg_log2FC))
            return(m)
        })
        names(marker_list1) <- clusters1
        tmp_clusters1 <- clusters1
        # combine cluster 1 and 2
        # marker_list1[["1_2"]] <- rbind(marker_list1[["1"]], marker_list1[["2"]])
        # marker_list1[["1_2"]] <- marker_list1[["1_2"]] %>% arrange(desc(avg_log2FC))
        # marker_list1[["1"]] <- NULL
        # marker_list1[["2"]] <- NULL
        # tmp_clusters1 <- setdiff(c(clusters1, "1_2"), c("1", "2"))

        marker_list2 <- lapply(clusters2, function(i) {
            f <- file.path(marker_path, name2, paste0("markers_", i, "_up_", nup, "_down_", ndown, pc_suf, ".csv"))
            m <- read.csv(f, header = TRUE, row.names = 1, sep = ",", comment.char = "#") %>% arrange(desc(avg_log2FC))
            return(m)
        })
        names(marker_list2) <- clusters2

        # MARKERS
        for (n_top in c(30, 50, 100, 150)) {
            jsi_mat <- matrix(0, nrow = length(tmp_clusters1), ncol = length(clusters2))
            rownames(jsi_mat) <- tmp_clusters1
            colnames(jsi_mat) <- clusters2

            for (i in tmp_clusters1) {
                for (j in clusters2) {
                    m1 <- marker_list1[[i]]
                    m2 <- marker_list2[[j]]

                    actual_n <- min(n_top, nrow(m1), nrow(m2))

                    if (actual_n == 0) {
                        jsi_mat[i, j] <- 0
                        jsi_mat[j, i] <- 0
                        next
                    }

                    m1 <- m1 %>% head(actual_n) %>% rownames
                    m2 <- m2 %>% head(actual_n) %>% rownames
                    jsi_mat[i, j] <- length(intersect(m1, m2)) / length(union(m1, m2))
                }
            }

            col_fun <- circlize::colorRamp2(c(0, 0.5, 1), RColorBrewer::brewer.pal(9, "Greens")[c(1, 5, 9)])
            pdf(file.path(figure_path, paste0("marker_jsi_", nup, "_", ndown, pc_suf, "_top_", n_top, ".pdf")), width = 10, height = 10)
            print(
                Heatmap(
                    jsi_mat,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", jsi_mat[i, j]), x, y, gp = gpar(fontsize = 13))
                    },
                    col = col_fun,
                    column_title = paste0("JSI - Top ", n_top, " genes - expression derived from region +", nup, "/-", ndown, pc_suf)
                )
            )
            dev.off()
        }

        # ENRICHMENT
        enrich_list1 <- lapply(clusters1, function(i) {
            f <- file.path(enrich_path, name1, paste0("enrichment_", i, "_up_", nup, "_down_", ndown, ".csv"))
            m <- read.csv(f, header = TRUE, row.names = 1, sep = ",", comment.char = "#") %>% arrange(p_value)
            # m <- m %>% filter(p_value < 0.05)
            return(m)
        })
        names(enrich_list1) <- clusters1
        temp_clusters1 <- clusters1

        # combine cluster 1 and 2
        # if (!is.null(enrich_list1[["1"]]) && !is.null(enrich_list1[["2"]])) {
        #     enrich_list1[["1_2"]] <- rbind(enrich_list1[["1"]], enrich_list1[["2"]])
        #     enrich_list1[["1_2"]] <- enrich_list1[["1_2"]] %>% arrange(p_value)
        #     enrich_list1[["1"]] <- NULL
        #     enrich_list1[["2"]] <- NULL
        #     temp_clusters1 <- setdiff(c(clusters1, "1_2"), c("1", "2"))
        # } else {
        #     temp_clusters1 <- clusters1
        # }

        enrich_list2 <- lapply(clusters2, function(i) {
            f <- file.path(enrich_path, name2, paste0("enrichment_", i, "_up_", nup, "_down_", ndown, ".csv"))
            m <- read.csv(f, header = TRUE, row.names = 1, sep = ",", comment.char = "#") %>% arrange(p_value)
            # m <- m %>% filter(p_value < 0.05)
            return(m)
        })
        names(enrich_list2) <- clusters2

        jsi_mat <- matrix(0, nrow = length(temp_clusters1), ncol = length(clusters2))
        rownames(jsi_mat) <- temp_clusters1
        colnames(jsi_mat) <- clusters2

        for (i in temp_clusters1) {
            for (j in clusters2) {
                m1 <- enrich_list1[[i]]
                m2 <- enrich_list2[[j]]
                if (is.null(m1) || is.null(m2)) {
                    jsi_mat[i, j] <- 0
                    jsi_mat[j, i] <- 0
                    next
                }
                n_top <- min(30, nrow(m1), nrow(m2))

                if (actual_n == 0) {
                    jsi_mat[i, j] <- 0
                    jsi_mat[j, i] <- 0
                    next
                }

                m1 <- m1 %>% head(n_top) %>% rownames
                m2 <- m2 %>% head(n_top) %>% rownames

                jsi_mat[i, j] <- length(intersect(m1, m2)) / length(union(m1, m2))
            }
        }

        col_fun <- circlize::colorRamp2(c(0, 0.5, 1), RColorBrewer::brewer.pal(9, "Greens")[c(1, 5, 9)])
        pdf(file.path(figure_path, paste0("enrich_jsi_", nup, "_", ndown, pc_suf, "_top_", n_top, ".pdf")), width = 10, height = 10)
        print(
            Heatmap(
                jsi_mat,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", jsi_mat[i, j]), x, y, gp = gpar(fontsize = 13))
                },
                col = col_fun,
                column_title = paste0("JSI - Top ", n_top, " genes - enrichment derived from region +", nup, "/-", ndown, pc_suf)
            )
        )
        dev.off()
    }
}
