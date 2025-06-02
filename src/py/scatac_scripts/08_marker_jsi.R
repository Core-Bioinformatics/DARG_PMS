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

first_group <- "new_aggregate"
second_group <- "old_aggregate"

marker_folder <- file.path(project_folder, "output", "peaks", "markers")
prefix <- "markers_all_"

list_files_1 <- list.files(file.path(marker_folder, first_group), pattern = paste0(prefix, ".*\\.csv$"), full.names = TRUE)
list_markers1 <- lapply(list_files_1, function(x) {
    read.table(x, header = TRUE, sep = "\t", comment.char = "#", fill = TRUE, quote = "")
})
names(list_markers1) <- sapply(list_files_1, function(x) {
    cl1_name <- strsplit(basename(x), "_")[[1]][3]
    cl1_name <- strsplit(cl1_name, "\\.")[[1]][1]
    return(cl1_name)
})
list_files_2 <- list.files(file.path(marker_folder, second_group), pattern = paste0(prefix, ".*\\.csv$"), full.names = TRUE)
list_markers2 <- lapply(list_files_2, function(x) {
    read.table(x, header = TRUE, sep = "\t", comment.char = "#", fill = TRUE, quote = "")
})
names(list_markers2) <- sapply(list_files_2, function(x) {
    cl2_name <- strsplit(basename(x), "_")[[1]][3]
    cl2_name <- strsplit(cl2_name, "\\.")[[1]][1]
    return(cl2_name)
})

peak_overlap_relaxed <- function(pk1, pk2, coverage_thresh = 0.00001) {
    chr1 <- strsplit(pk1, "-")[[1]][1]
    st1 <- as.integer(strsplit(pk1, "-")[[1]][2])
    end1 <- as.integer(strsplit(pk1, "-")[[1]][3])
    l1 <- end1 - st1
    m1 <- (end1 + st1) %/% 2

    chr2 <- strsplit(pk2, "-")[[1]][1]
    st2 <- as.integer(strsplit(pk2, "-")[[1]][2])
    end2 <- as.integer(strsplit(pk2, "-")[[1]][3])
    l2 <- end2 - st2
    m2 <- (end2 + st2) %/% 2
    if (chr1 != chr2) {
        return(FALSE)
    }

    if (st1 <= st2) {
        if (st2 >= end1) {
            return(FALSE)
        }

        min_l <- min(l1, l2)
        intersect_l <- min(end1, end2) - st2
        if (intersect_l / min_l >= coverage_thresh) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    return(peak_overal_relaxed(pk2, pk1, coverage_thresh))
}

peak_overlap_strict <- function(pk1, pk2) {
    chr1 <- strsplit(pk1, "-")[[1]][1]
    st1 <- as.integer(strsplit(pk1, "-")[[1]][2])
    end1 <- as.integer(strsplit(pk1, "-")[[1]][3])
    l1 <- end1 - st1
    m1 <- (end1 + st1) %/% 2

    chr2 <- strsplit(pk2, "-")[[1]][1]
    st2 <- as.integer(strsplit(pk2, "-")[[1]][2])
    end2 <- as.integer(strsplit(pk2, "-")[[1]][3])
    l2 <- end2 - st2
    m2 <- (end2 + st2) %/% 2
    if (chr1 != chr2) {
        return(FALSE)
    }

    if (l1 > l2) {
        if (st1 <= m2 && end1 >= m2) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }

    if (st2 <= m1 && end2 >= m1) {
        return(TRUE)
    }

    return(FALSE)
}

check_peak_overlap <- function(pklist1, pklist2) {
    i1 <- seq_along(pklist1)
    i2 <- c()

    for (i in seq_along(pklist2)) {
        pk <- pklist2[[i]]
        pk_matches <- sapply(pklist1, function(x) peak_overlap_relaxed(pk, x))
        if (sum(pk_matches) == 0) {
            i2 <- c(i2, i + length(i1))
            next
        }
        pk_matches <- which(pk_matches == TRUE)

        if (length(intersect(pk_matches, i2)) == length(pk_matches)) {
            i2 <- c(i2, pk_matches[1])
            next
        }

        i2 <- c(i2, setdiff(pk_matches, i2)[1])
    }

    return(list(i1, i2))
}

for (top_n_peaks in c(30, 50, 100, 150, 200, 250, 300)) {
    for (gene_type in c("pc", "all")) {
        jsi_matrix <- matrix(0, nrow = length(list_files_1), ncol = length(list_files_2))
        colnames(jsi_matrix) <- names(list_markers2)
        rownames(jsi_matrix) <- names(list_markers1)

        for (i in seq_along(list_markers1)) {
            current_mk1 <- list_markers1[[i]]
            if (gene_type == "pc") {
                current_mk1 <- current_mk1 %>% filter(Gene.Type == "protein-coding")
            }
            n1 <- nrow(current_mk1)
            current_mk1 <- sapply(seq_len(nrow(current_mk1)), function(x) paste(current_mk1$Chr[x], current_mk1$Start[x], current_mk1$End[x], sep = "-"))[seq_len(top_n_peaks)]

            for (j in seq_along(list_markers2)) {
                current_mk2 <- list_markers2[[j]]
                if (gene_type == "pc") {
                    current_mk2 <- current_mk2 %>% filter(Gene.Type == "protein-coding")
                }
                n2 <- nrow(current_mk2)
                current_mk2 <- sapply(seq_len(nrow(current_mk2)), function(x) paste(current_mk2$Chr[x], current_mk2$Start[x], current_mk2$End[x], sep = "-"))[seq_len(top_n_peaks)]
                n_actual_peaks <- min(n1, n2)
                n_actual_peaks <- min(n_actual_peaks, top_n_peaks)

                overlap_index_list <- check_peak_overlap(current_mk1[seq_len(n_actual_peaks)], current_mk2[seq_len(n_actual_peaks)])
                jsi_matrix[i, j] <- length(intersect(overlap_index_list[[1]], overlap_index_list[[2]])) / length(union(overlap_index_list[[1]], overlap_index_list[[2]]))
            }
        }


        col_fun <- circlize::colorRamp2(c(0, 0.5, 1), RColorBrewer::brewer.pal(9, "Greens")[c(1, 5, 9)])
        pdf(file.path(marker_folder, paste0("jsi_", prefix, first_group, "_", second_group, "_", gene_type, "_genes_", top_n_peaks, ".pdf")), width = 10, height = 10)
        print(Heatmap(
            jsi_matrix,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", jsi_matrix[i, j]), x, y, gp = gpar(fontsize = 13))
            },
            col = col_fun,
            column_title = paste0("JSI - Top ", top_n_peaks, " ", gene_type, " genes")
        ))
        dev.off()
    }
}
