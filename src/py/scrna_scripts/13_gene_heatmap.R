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

figures_path <- file.path(project_folder, "output", "RNA")
seurat_path <- file.path(project_folder, "objects", "R", "seurat")
peak_path <- "/servers/iss-corescratch/lp488/bongsoo/snatac_analysis_andi/output/peaks"

so <- list(
    "healthy" = qread(file.path(seurat_path, "healthy_norm.qs"), nthreads = ncores),
    "pms" = qread(file.path(seurat_path, "pms_norm.qs"), nthreads = ncores)
)
so$healthy <- so$healthy[, so$healthy$treatment == "non_treated"]
so$pms <- so$pms[, so$pms$treatment == "non_treated"]
so$pms$stable_clusters <- so$pms$clusters_from_samples

bgsoo_ga_genes <- read.table(file.path(peak_path, "bongsoo_DARs_cluster25_vs_013_GA_3kTSS.txt"), sep = "\t", header = FALSE, comment.char = "#")
bgsoo_ga_genes <- unique(bgsoo_ga_genes$V5)

bgsoo_la_genes <- read.table(file.path(peak_path, "bongsoo_DARs_cluster25_vs_013_LA_3kTSS.txt"), sep = "\t", header = FALSE, comment.char = "#")
bgsoo_la_genes <- unique(bgsoo_la_genes$V5)

custom_gene_list <- c("NOD1", "WDR49", "PLP1", "RBMS3", "TBC1D1", "EIF3G", "IGSF11", "C5orf15", "PARP14", "OAF", "TENT5A", "DDX60", "ZBTB32", "PLEKHA4", "SULF1", "IFIT3", "CEACAM1", "BST2", "LACTB", "MYO1E", "CALHM5", "TRIM26", "CLMN", "IL17RA", "IFI6", "IFIT5", "NREP", "KIAA1217", "APOL6", "TDRD7", "SMARCA5", "IFIT1", "RTL9", "RANBP3L", "SH2B3", "APOL2", "USP18", "IFI44", "OAS2", "WDR25", "TMEM219", "NR6A1", "GCNT1", "GTPBP1", "HAX1", "SYN3", "SLC6A15", "UBFD1", "SERPINH1", "PRKAR1B", "NQO2")

###### Heatmap with list of DE genes from Bongsoo ######
for (stype in names(so)) {
    expr_matrix <- GetAssayData(so[[stype]], layer = "data")
    # GA
    new_expr_matrix <- expr_matrix[intersect(bgsoo_ga_genes, rownames(expr_matrix)), ]
    new_expr_matrix <- new_expr_matrix[rowSums(new_expr_matrix) > 0, ]

    uniq_cl <- stringr::str_sort(unique(so[[stype]]$stable_clusters), numeric = TRUE)
    psdb_mat <- matrix(0, nrow = nrow(new_expr_matrix), ncol = length(uniq_cl))
    rownames(psdb_mat) <- NULL
    colnames(psdb_mat) <- uniq_cl

    
    for (i in uniq_cl) {
        psdb_mat[, i] <- rowMeans(new_expr_matrix[, so[[stype]]$stable_clusters == i, drop = FALSE])
    }

    psdb_mat <- t(scale(t(psdb_mat)))

    pdf(file.path(figures_path, paste0("heatmap_de_genes_bongsoo_ga_RNA_", stype, ".pdf")), width = 10, height = 10)
    print(
        Heatmap(
            psdb_mat,
            cluster_columns = FALSE,
            show_row_dend = FALSE
        )
    )
    dev.off()

    # LA
    new_expr_matrix <- expr_matrix[intersect(bgsoo_la_genes, rownames(expr_matrix)), ]
    new_expr_matrix <- new_expr_matrix[rowSums(new_expr_matrix) > 0, ]

    uniq_cl <- stringr::str_sort(unique(so[[stype]]$stable_clusters), numeric = TRUE)
    psdb_mat <- matrix(0, nrow = nrow(new_expr_matrix), ncol = length(uniq_cl))
    rownames(psdb_mat) <- NULL
    colnames(psdb_mat) <- uniq_cl

    
    for (i in uniq_cl) {
        psdb_mat[, i] <- rowMeans(new_expr_matrix[, so[[stype]]$stable_clusters == i, drop = FALSE])
    }

    psdb_mat <- t(scale(t(psdb_mat)))

    pdf(file.path(figures_path, paste0("heatmap_de_genes_bongsoo_la_RNA_", stype, ".pdf")), width = 10, height = 10)
    print(
        Heatmap(
            psdb_mat,
            cluster_columns = FALSE,
            show_row_dend = FALSE
        )
    )
    dev.off()

    # custom
    new_expr_matrix <- expr_matrix[intersect(custom_gene_list, rownames(expr_matrix)), ]
    new_expr_matrix <- new_expr_matrix[rowSums(new_expr_matrix) > 0, ]
    print(setdiff(custom_gene_list, rownames(new_expr_matrix)))

    uniq_cl <- stringr::str_sort(unique(so[[stype]]$stable_clusters), numeric = TRUE)
    psdb_mat <- matrix(0, nrow = nrow(new_expr_matrix), ncol = length(uniq_cl))
    rownames(psdb_mat) <- rownames(new_expr_matrix)
    colnames(psdb_mat) <- uniq_cl
    
    for (i in uniq_cl) {
        psdb_mat[, i] <- rowMeans(new_expr_matrix[, so[[stype]]$stable_clusters == i, drop = FALSE])
    }

    psdb_mat <- t(scale(t(psdb_mat)))

    pdf(file.path(figures_path, paste0("heatmap_custom_genes_RNA_", stype, ".pdf")), width = 10, height = 10)
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
