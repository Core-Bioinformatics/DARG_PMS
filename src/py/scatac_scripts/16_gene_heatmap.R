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

figures_path <- file.path(project_folder, "output", "figures")
expr_path <- file.path(project_folder, "objects", "expr_matrix")
peak_path <- file.path(project_folder, "output", "peaks")
seurat_path <- file.path(project_folder, "objects", "R", "seurat")

so <- qread(file.path(seurat_path, "aggregate_norm.qs"), nthreads = ncores)

bgsoo_ga_genes <- read.table(file.path(peak_path, "bongsoo_DARs_cluster25_vs_013_GA_3kTSS.txt"), sep = "\t", header = FALSE, comment.char = "#")
bgsoo_ga_genes <- unique(bgsoo_ga_genes$V5)

bgsoo_la_genes <- read.table(file.path(peak_path, "bongsoo_DARs_cluster25_vs_013_LA_3kTSS.txt"), sep = "\t", header = FALSE, comment.char = "#")
bgsoo_la_genes <- unique(bgsoo_la_genes$V5)


###### Heatmap with list of DE genes from Bongsoo ######
for (nup_comb in list(
    c(100, 1000),
    c(3000, 3000),
    c(50000, 50000)
)) {
    nup <- nup_comb[1]
    ndown <- nup_comb[2]
    expr_matrix <- qread(file.path(expr_path, paste0("peak_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")), nthreads = ncores)
    new_expr_matrix <- expr_matrix[intersect(bgsoo_ga_genes, rownames(expr_matrix)), ]

    uniq_cl <- stringr::str_sort(unique(so$stable_clusters_9), numeric = TRUE)
    psdb_mat <- matrix(0, nrow = nrow(new_expr_matrix), ncol = length(uniq_cl))
    rownames(psdb_mat) <- NULL
    colnames(psdb_mat) <- uniq_cl

    
    for (i in uniq_cl) {
        psdb_mat[, i] <- rowMeans(new_expr_matrix[, so$stable_clusters_9 == i, drop = FALSE])
    }
    # combine clusters 1, 4 and 2
    psdb_mat <- cbind(psdb_mat, rowMeans(psdb_mat[, c("1", "4", "2")]))
    colnames(psdb_mat)[ncol(psdb_mat)] <- "1_4_2"
    psdb_mat <- psdb_mat[ , setdiff(colnames(psdb_mat), c("1", "4", "2"))]

    # scale row-wise
    psdb_mat <- t(scale(t(psdb_mat)))


    pdf(file.path(figures_path, paste0("heatmap_de_genes_bongsoo_ga_ATAC_up_", nup, "_down_", ndown, ".pdf")), width = 10, height = 10)
    print(
        Heatmap(
            psdb_mat,
            cluster_columns = FALSE,
            show_row_dend = FALSE
        )
    )
    dev.off()

    new_expr_matrix <- expr_matrix[intersect(bgsoo_la_genes, rownames(expr_matrix)), ]

    uniq_cl <- stringr::str_sort(unique(so$stable_clusters_9), numeric = TRUE)
    psdb_mat <- matrix(0, nrow = nrow(new_expr_matrix), ncol = length(uniq_cl))
    rownames(psdb_mat) <- NULL
    colnames(psdb_mat) <- uniq_cl

    
    for (i in uniq_cl) {
        psdb_mat[, i] <- rowMeans(new_expr_matrix[, so$stable_clusters_9 == i, drop = FALSE])
    }
    # combine clusters 1,4 and 2
    psdb_mat <- cbind(psdb_mat, rowMeans(psdb_mat[, c("1", "4", "2")]))
    colnames(psdb_mat)[ncol(psdb_mat)] <- "1_4_2"
    psdb_mat <- psdb_mat[ , setdiff(colnames(psdb_mat), c("1", "4", "2"))]

    # scale row-wise
    psdb_mat <- t(scale(t(psdb_mat)))


    pdf(file.path(figures_path, paste0("heatmap_de_genes_bongsoo_la_ATAC_up_", nup, "_down_", ndown, ".pdf")), width = 10, height = 10)
    print(
        Heatmap(
            psdb_mat,
            cluster_columns = FALSE,
            show_row_dend = FALSE
        )
    )
    dev.off()
}

custom_gene_list <- c("NOD1", "WDR49", "PLP1", "RBMS3", "TBC1D1", "EIF3G", "IGSF11", "C5orf15", "PARP14", "OAF", "TENT5A", "DDX60", "ZBTB32", "PLEKHA4", "SULF1", "IFIT3", "CEACAM1", "BST2", "LACTB", "MYO1E", "CALHM5", "TRIM26", "CLMN", "IL17RA", "IFI6", "IFIT5", "NREP", "KIAA1217", "APOL6", "TDRD7", "SMARCA5", "IFIT1", "RTL9", "RANBP3L", "SH2B3", "APOL2", "USP18", "IFI44", "OAS2", "WDR25", "TMEM219", "NR6A1", "GCNT1", "GTPBP1", "HAX1", "SYN3", "SLC6A15", "UBFD1", "SERPINH1", "PRKAR1B", "NQO2")
for (nup_comb in list(
    c(100, 1000),
    c(3000, 3000),
    c(50000, 50000)
)) {
    nup <- nup_comb[1]
    ndown <- nup_comb[2]
    expr_matrix <- qread(file.path(expr_path, paste0("peak_based_gene_matrix_up_", nup, "_down_", ndown, ".qs")), nthreads = ncores)
    new_expr_matrix <- expr_matrix[intersect(custom_gene_list, rownames(expr_matrix)), ]
    print(setdiff(custom_gene_list, rownames(new_expr_matrix)))

    uniq_cl <- stringr::str_sort(unique(so$stable_clusters_9), numeric = TRUE)
    psdb_mat <- matrix(0, nrow = nrow(new_expr_matrix), ncol = length(uniq_cl))
    rownames(psdb_mat) <- rownames(new_expr_matrix)
    colnames(psdb_mat) <- uniq_cl

    
    for (i in uniq_cl) {
        psdb_mat[, i] <- rowMeans(new_expr_matrix[, so$stable_clusters_9 == i, drop = FALSE])
    }

    # combine clusters 1,4 and 2
    psdb_mat <- cbind(psdb_mat, rowMeans(psdb_mat[, c("1", "4", "2")]))
    colnames(psdb_mat)[ncol(psdb_mat)] <- "1_4_2"
    psdb_mat <- psdb_mat[ , setdiff(colnames(psdb_mat), c("1", "4", "2"))]

    # scale row-wise
    psdb_mat <- t(scale(t(psdb_mat)))

    pdf(file.path(figures_path, paste0("heatmap_custom_genes_ATAC_up_", nup, "_down_", ndown, ".pdf")), width = 10, height = 10)
    print(
        Heatmap(
            psdb_mat,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_row_dend = FALSE
        )
    )
    dev.off()
}
