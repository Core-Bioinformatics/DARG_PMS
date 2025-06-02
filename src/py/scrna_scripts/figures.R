setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(ClustAssess)
library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)

ca_folder <- file.path(project_folder, "objects", "R", "clustassess")
ca_app_folder <- file.path(project_folder, "objects", "shiny_apps", "clustassess")
so_folder <- file.path(project_folder, "objects", "R", "seurat")
figure_folder <- file.path(project_folder, "figures")

if (!dir.exists(figure_folder)) {
  dir.create(figure_folder, recursive = TRUE)
}

pms_ca <- qs::qread(file.path(ca_folder, "pms.qs"), nthreads = ncores)
abt_pms_ca <- qs::qread(file.path(ca_folder, "105_2g_ABT_263.qs"), nthreads = ncores)
dmso_pms_ca <- qs::qread(file.path(ca_folder, "105_2g_DMSO_sub_0_4.qs"), nthreads = ncores)

pms <- qs::qread(file.path(so_folder, "pms_norm.qs"), nthreads = ncores)
pms_ecc <- c(
    get_clusters_from_clustassess_object(
        clustassess_object = abt_pms_ca,
        feature_type = "Highly_Variable",
        feature_size = 2000,
        clustering_method = "SLM",
        nclusters = 11
    )[[1]]$ecc,
    get_clusters_from_clustassess_object(
        clustassess_object = dmso_pms_ca,
        feature_type = "Highly_Variable",
        feature_size = 2500,
        clustering_method = "SLM",
        nclusters = 9 
    )[[1]]$ecc
)
names(pms_ecc) <- c(
    colnames(pms)[pms$sample_name == "105_2g_ABT_263"],
    colnames(pms)[pms$sample_name == "105_2g_DMSO_sub_0_4"]
)
pms_ecc <- pms_ecc[colnames(pms)]
pms$ecc <- pms_ecc

ctrl_ca <- qs::qread(file.path(ca_folder, "healthy.qs"), nthreads = ncores)
ctrl <- qs::qread(file.path(so_folder, "healthy_norm.qs"), nthreads = ncores)

pms@reductions$umap@cell.embeddings <- pms_ca$Highly_Variable$"2000"$umap
colnames(pms@reductions$umap@cell.embeddings) <- paste0("umap_", 1:2)

ctrl@reductions$umap@cell.embeddings <- ctrl_ca$Highly_Variable$"3000"$umap
colnames(ctrl@reductions$umap@cell.embeddings) <- paste0("umap_", 1:2)
ctrl$stable_clusters <- get_clusters_from_clustassess_object(
    clustassess_object = ctrl_ca,
    feature_type = "Highly_Variable",
    feature_size = 3000,
    clustering_method = "SLM",
    nclusters = 13
)[[1]]$partitions[[1]]$mb - 1

ctrl$ecc <- get_clusters_from_clustassess_object(
    clustassess_object = ctrl_ca,
    feature_type = "Highly_Variable",
    feature_size = 3000,
    clustering_method = "SLM",
    nclusters = 13
)[[1]]$ecc

##### A FIGURES
p <- DimPlot(pms, group.by = "clusters_from_samples") +
    ggtitle("") +
    theme(
        aspect.ratio = 1
    )

pdf(file.path(figure_folder, "1_a_pms.pdf"))
print(p)
dev.off()

p <- DimPlot(ctrl, group.by = "stable_clusters") +
    ggtitle("") +
    theme(
        aspect.ratio = 1
    )

pdf(file.path(figure_folder, "1_a_ctrl.pdf"))
print(p)
dev.off()

##### B FIGURES
p <- FeaturePlot(ctrl, features = "ecc") +
    scale_colour_viridis_c() +
    ggtitle("")

pdf(file.path(figure_folder, "1_b_ctrl.pdf"))
print(p)
dev.off()

p <- FeaturePlot(pms, features = "ecc") +
    scale_colour_viridis_c() +
    ggtitle("")

pdf(file.path(figure_folder, "1_b_pms.pdf"))
print(p)
dev.off()

##### C FIGURES
type1 = c( 'SOX2', 'NES', 'PAX6', 'SOX4', 'SOX11', 'PTPRZ1', 'FABP7', 'TUBB2B', 'YBX1')
type2 = c( 'HES1', 'SLC1A3', 'VIM', 'CKB', 'MEG3', 'PTN', 'CLU', 'MT1X', 'CST3', 'ID4')
type3 = c( 'SOX4', 'SOX11', 'ASCL1', 'ELAVL4', 'NHLH1', 'CBFA2T2', 'DCX')

set1 = c( 'OLIG1', 'OLIG2', 'MBP', 'PLP', 'MOG', 'PDGRA', 'CSPG4')
set2 = c( 'AQP4', 'ALDH1L1', 'GFAP', 'S100B')
set3 = c( 'CALB', 'CCK', 'PVALB', 'VIP')

plot_genes <- function(so, list_of_genes, title, expr_threshold=0, relaxation=0) {
    expr_matrix <- GetAssayData(so, layer = "data")
    expr_matrix <- as.matrix(expr_matrix[rownames(expr_matrix) %in% list_of_genes,])
    new_matrix <- matrixStats::colSums2(expr_matrix > expr_threshold) >= (length(list_of_genes) - relaxation)
    a = so@reductions$umap@cell.embeddings
    a_2 = cbind(a, list(new_matrix)[[1]])
    colnames(a_2) = c('UMAP_1', 'UMAP_2', 'Keep')
    
    df = as.data.frame(a_2)
    x = ggplot(df %>% arrange(Keep), aes(x=UMAP_1, y=UMAP_2, col)) + geom_point(size=3, aes(color=factor(Keep))) +
        ggtitle(title) + theme_classic() + scale_color_manual(values=c('#132B43', '#56B1F7')) +
        labs(color='In gene set')
    print(x)
}

plot_genes_average <- function(so, list_of_genes, title, expr_threshold=0, relaxation=0, midpoint=0) {
    expr_matrix <- GetAssayData(so, layer = "data")
    expr_matrix <- as.matrix(expr_matrix[rownames(expr_matrix) %in% list_of_genes,])
    new_matrix = matrixStats::colMeans2(expr_matrix)
    a = so@reductions$umap@cell.embeddings
    a_2 = cbind(a, list(new_matrix)[[1]])
    colnames(a_2) = c('umap_1', 'umap_2', 'Mean gene set expression')
    
    df = as.data.frame(a_2)
    x = ggplot(df %>% arrange(`Mean gene set expression`), aes(x=umap_1, y=umap_2)) + 
        geom_point(aes(color=`Mean gene set expression`, alpha=`Mean gene set expression`), size=4) +
        scale_color_gradient2(low='#450d59', mid='#1f8f8d', high='#fde726', midpoint=midpoint) +
        ggtitle(title) + theme_classic()
    return(list("df" = df, "g" = x))
}

ggsave(
    file.path(figure_folder, "1_c_ctrl.pdf"),
    plot_genes(ctrl, type1, "Type 1", expr_threshold=1, relaxation=2)
)

ggsave(
    file.path(figure_folder, "1_d_ctrl.pdf"),
    plot_genes(ctrl, type2, "Type 2", expr_threshold=0.5, relaxation=2)
)

ggsave(
    file.path(figure_folder, "1_e_ctrl.pdf"),
    plot_genes(ctrl, type3, "Type 3", expr_threshold=0.5, relaxation=1)
)


ggsave(
    file.path(figure_folder, "1_c_pms.pdf"),
    plot_genes(pms, type1, "Type 1", expr_threshold=1, relaxation=2)
)

ggsave(
    file.path(figure_folder, "1_d_pms.pdf"),
    plot_genes(pms, type2, "Type 2", expr_threshold=0.5, relaxation=2)
)

ggsave(
    file.path(figure_folder, "1_e_pms.pdf"),
    plot_genes(pms, type3, "Type 3", expr_threshold=0.5, relaxation=1)
)

##### F FIGURES
cscheme <- c("#A87B95", "#59BFC7")
p <- DimPlot(pms, group.by = "treatment") +
    ggtitle("") +
    scale_colour_manual(values = setNames(cscheme, c("non_treated", "senolytic"))) +
    theme(
        aspect.ratio = 1
    )

pdf(file.path(figure_folder, "1_f1_pms.pdf"))
print(p)
dev.off()

df <- pms[[c("treatment", "clusters_from_samples")]]
df <- melt(table(df[, 1], df[, 2]))
df <- df %>%
    group_by(Var2) %>%
    mutate(percentage = value / sum(value))

p <- ggplot(df, aes(fill = Var1, y = Var2, x = percentage)) +
    geom_bar(stat = "identity", colour = "black", width = 0.5) +
    scale_fill_manual(values = setNames(cscheme, c("non_treated", "senolytic"))) +
    theme_classic() +
    labs(x = "Cell Proportion", y = "Clusters", fill = "Treatment") +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20)
    )

pdf(file.path(figure_folder, "1_f2_pms.pdf"), height = 10)
print(p)
dev.off()

p <- DimPlot(ctrl, group.by = "treatment") +
    ggtitle("") +
    scale_colour_manual(values = setNames(cscheme, c("non_treated", "senolytic"))) +
    theme(
        aspect.ratio = 1
    )

pdf(file.path(figure_folder, "1_f1_ctrl.pdf"))
print(p)
dev.off()

df <- ctrl[[c("treatment", "stable_clusters")]]
df <- melt(table(df[, 1], df[, 2]))
df$Var2 <- factor(df$Var2, levels = rev(stringr::str_sort(unique(df$Var2), numeric = TRUE)))
df <- df %>%
    group_by(Var2) %>%
    mutate(percentage = value / sum(value))

p <- ggplot(df, aes(fill = Var1, y = Var2, x = percentage)) +
    geom_bar(stat = "identity", colour = "black", width = 0.5) +
    scale_fill_manual(values = setNames(cscheme, c("non_treated", "senolytic"))) +
    theme_classic() +
    labs(x = "Cell Proportion", y = "Clusters", fill = "Treatment") +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20)
    )
pdf(file.path(figure_folder, "1_f2_ctrl.pdf"), height = 10)
print(p)
dev.off()

##### SLIDE 3

p <- DimPlot(pms, group.by = "Phase") +
    ggtitle("") +
    theme(
        aspect.ratio = 1
    )

pdf(file.path(figure_folder, "3_a_pms.pdf"))
print(p)
dev.off()

p <- DimPlot(ctrl, group.by = "Phase") +
    ggtitle("") +
    theme(
        aspect.ratio = 1
    )

pdf(file.path(figure_folder, "3_a_ctrl.pdf"))
print(p)
dev.off()

df <- pms[[c("Phase", "clusters_from_samples")]]
df <- melt(table(df[, 1], df[, 2]))
df <- df %>%
    group_by(Var2) %>%
    mutate(percentage = value / sum(value))

p <- ggplot(df, aes(fill = Var1, x = Var2, y = percentage)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    labs(y = "Cell Proportion", x = "Clusters", fill = "Phase") +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

pdf(file.path(figure_folder, "3_b_pms.pdf"), width = 10)
print(p)
dev.off()

df <- ctrl[[c("Phase", "stable_clusters")]]
df <- melt(table(df[, 1], df[, 2]))
df$Var2 <- factor(df$Var2, levels = stringr::str_sort(unique(df$Var2), numeric = TRUE))
df <- df %>%
    group_by(Var2) %>%
    mutate(percentage = value / sum(value))

p <- ggplot(df, aes(fill = Var1, x = Var2, y = percentage)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    labs(y = "Cell Proportion", x = "Clusters", fill = "Phase") +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

pdf(file.path(figure_folder, "3_b_ctrl.pdf"), width = 10)
print(p)
dev.off()

##### SLIDE 4

ggsave(
    file.path(figure_folder, "4_a_ctrl.pdf"),
    plot_genes_average(ctrl, set1, "Set 1", midpoint=0.15)$g
)

ggsave(
    file.path(figure_folder, "4_b_ctrl.pdf"),
    plot_genes_average(ctrl, set2, "Set 2", midpoint=0.2)$g
)

ggsave(
    file.path(figure_folder, "4_c_ctrl.pdf"),
    plot_genes_average(ctrl, set3, "Set 3", midpoint=0.15)$g
)

ggsave(
    file.path(figure_folder, "4_a_pms.pdf"),
    plot_genes_average(pms, set1, "Set 1", midpoint=0.15)$g
)

ggsave(
    file.path(figure_folder, "4_b_pms.pdf"),
    plot_genes_average(pms, set2, "Set 2", midpoint=0.2)$g
)

ggsave(
    file.path(figure_folder, "4_c_pms.pdf"),
    plot_genes_average(pms, set3, "Set 3", midpoint=0.2)$g
)

##### SLIDE 5

##### SLIDE 6
ifn <- c("BST2", "SAMHD1", "IFIT1", "IFIT2", "IFIT3", "IFITM3")
notch1 <- c("CDK8", "ADAM10", "MAML2", "NEURL1B", "HDAC9", "TBL1XR1", "TLE4", "MIB1", "JAG1", "NOTCH1", "DTX4", "MAMLD1", "TLE1", "MAML3", "HES5", "MAMLD1")

ggsave(
    file.path(figure_folder, "6_ifn_ctrl.pdf"),
    plot_genes(
        ctrl,
        ifn,
        "IFN",
        expr_threshold = 0,
        relaxation = 2
    )
)

ggsave(
    file.path(figure_folder, "6_ifn_pms.pdf"),
    plot_genes(
        pms,
        ifn,
        "IFN",
        expr_threshold = 0,
        relaxation = 2
    )
)

ggsave(
    file.path(figure_folder, "6_notch1_ctrl.pdf"),
    plot_genes(
        ctrl,
        notch1,
        "Notch1",
        expr_threshold = 0.6,
        relaxation = 8
    )
)

ggsave(
    file.path(figure_folder, "6_notch1_pms.pdf"),
    plot_genes(
        pms,
        notch1,
        "Notch1",
        expr_threshold = 0.6,
        relaxation = 8
    )
)


#### QC Figures #####
aggr <- qs::qread(file.path(so_folder, "aggregate_norm.qs"), nthreads = ncores)

pdf(file.path(figure_folder, "qc_aggregate.pdf"), width = 10, height = 10)
print(VlnPlot(
    aggr,
    features = c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_rp"),
    ncol = 4,
    pt.size = 0,
    log = FALSE, 
    cols = c("#f9af49", "#768b57", "#fe3619", "#01b876"),
) + NoLegend())
dev.off()

#### bubbleplot
list_genes <- c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2", "ISG15", "IFIT1", "IFIT2", "CDKN1A", "CDKN2A", "TP53", "B2M", "STAT3", "RELA", "NFKB1", "SP1", "TGFB1", "BMP1", "BST2", "SAMHD1", "IFITM3", "USP18", "HLA-E", "TRIM21", "OASL", "DDX60", "DHX58", "FN1")

list_genes_list <- list(
    "fn1" = c("FN1", "TNF", "CDKN1A", "CDKN2A", "TP53", "TGFB1", "GDF15", "HMGB1", "TNFRSF1B", "TIMP2", "MMP2", "IL15", "IL18", "IL32", "ITGA2", "IGFBP2", "IGFBP4", "IGFBP5", "IGF1"),
    "usp18" = c("USP18", "IFITM3", "SAMHD1", "BST2", "IFIT1", "IFIT2", "ISG15", "NFKB1"),
    "trim21" = c("TRIM21", "OASL", "DDX60", "DHX58"),
    "cxcl12" = c("CXCL12", "CD9", "C3")
)
list_genes_list[["fn1_usp18_trim21_cxcl12"]] <- unique(unlist(list_genes_list))

expr_mat <- GetAssayData(pms, layer = "data")

for (prefix in names(list_genes_list)) {
    list_genes <- intersect(rownames(expr_mat), list_genes_list[[prefix]])
    print(paste(prefix, length(list_genes)))
    temp_expr_mat <- as.matrix(expr_mat[list_genes,,drop = FALSE])

    dir.create(file.path(figure_folder, paste0(prefix, "_bubbleplot")), recursive = TRUE)

    # pseudobulk by clusters from samples
    psdb_expr_mat <- matrix(0, nrow = length(list_genes), ncol = length(unique(pms$clusters_from_samples)))
    freq_mat <- matrix(0, nrow = length(list_genes), ncol = length(unique(pms$clusters_from_samples)))
    rownames(psdb_expr_mat) <- list_genes
    rownames(freq_mat) <- list_genes
    colnames(psdb_expr_mat) <- unique(pms$clusters_from_samples)
    colnames(freq_mat) <- unique(pms$clusters_from_samples)

    for (i in unique(pms$clusters_from_samples)) {
        cells <- colnames(pms)[pms$clusters_from_samples == i]
        psdb_expr_mat[, i] <- rowMeans(temp_expr_mat[, cells, drop = FALSE])
        freq_mat[, i] <- rowSums(temp_expr_mat[, cells, drop = FALSE] > 0) / length(cells)
    }
    total_freq <- rowSums(psdb_expr_mat > 0)
    list_genes <- names(total_freq)[total_freq > 0]
    psdb_expr_mat <- psdb_expr_mat[list_genes, ]
    freq_mat <- freq_mat[list_genes, ]

    exp_df <- melt(psdb_expr_mat)
    freq_df <- melt(freq_mat)

    colnames(exp_df) <- c("gene", "cluster", "avg_expression")
    exp_df$freq <- freq_df$value

    head(exp_df)
    exp_df$cluster <- factor(exp_df$cluster, levels = stringr::str_sort(unique(exp_df$cluster), numeric = TRUE))
    exp_df$gene <- factor(exp_df$gene, levels = list_genes)


    pdf(file.path(figure_folder, paste0(prefix, "_bubbleplot"), "bubbleplot_all_norm.pdf"), width = 5 + nlevels(exp_df$cluster) * 0.1, height = 5 + nlevels(exp_df$gene) * 0.1)
    print(ggplot(exp_df, aes(x = cluster, y = gene, colour = avg_expression, size = freq)) +
        geom_point() +
        theme_classic() +
        scale_colour_viridis_c() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
    dev.off()

    exp_df <- exp_df %>% filter(cluster %in% c("dmso_cluster_8", "abt_cluster_9"))
    pdf(file.path(figure_folder, paste0(prefix, "_bubbleplot"), "bubbleplot_specific_norm.pdf"), width = 5 + nlevels(exp_df$cluster) * 0.1, height = 5 + nlevels(exp_df$gene) * 0.1)
    print(ggplot(exp_df, aes(x = cluster, y = gene, colour = avg_expression, size = freq)) +
        geom_point() +
        theme_classic() +
        scale_colour_viridis_c() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
    dev.off()

    # scale psdb mat
    psdb_expr_mat_temp <- t(scale(t(psdb_expr_mat)))
    psdb_expr_mat_temp[psdb_expr_mat > 3] <- 3
    psdb_expr_mat_temp[psdb_expr_mat < -3] <- -3
    exp_df <- melt(psdb_expr_mat_temp)

    colnames(exp_df) <- c("gene", "cluster", "avg_expression")
    exp_df$freq <- freq_df$value

    head(exp_df)
    exp_df$cluster <- factor(exp_df$cluster, levels = stringr::str_sort(unique(exp_df$cluster), numeric = TRUE))
    exp_df$gene <- factor(exp_df$gene, levels = list_genes)


    pdf(file.path(figure_folder, paste0(prefix, "_bubbleplot"), "bubbleplot_all_scaled.pdf"), width = 5 + nlevels(exp_df$cluster) * 0.1, height = 5 + nlevels(exp_df$gene) * 0.1)
    print(ggplot(exp_df, aes(x = cluster, y = gene, colour = avg_expression, size = freq)) +
        geom_point() +
        theme_classic() +
        scale_colour_gradient2(breaks = c(-3, 0, 3), low = "#4575b4", mid =  "#ffffff", high = "#d73027") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
    dev.off()


    target_clusters <-  c("dmso_cluster_8", "abt_cluster_9")
    psdb_expr_mat <- matrix(0, nrow = length(list_genes), ncol = length(target_clusters))
    freq_mat <- matrix(0, nrow = length(list_genes), ncol = length(target_clusters))
    rownames(psdb_expr_mat) <- list_genes
    rownames(freq_mat) <- list_genes
    colnames(psdb_expr_mat) <- target_clusters
    colnames(freq_mat) <- target_clusters

    for (i in colnames(psdb_expr_mat)) {
        cells <- colnames(pms)[pms$clusters_from_samples == i]
        psdb_expr_mat[, i] <- rowMeans(temp_expr_mat[list_genes, cells, drop = FALSE])
        freq_mat[, i] <- rowSums(temp_expr_mat[list_genes, cells, drop = FALSE] > 0) / length(cells)
    }
    total_freq <- rowSums(psdb_expr_mat > 0)
    list_genes <- names(total_freq)[total_freq > 0]
    psdb_expr_mat <- psdb_expr_mat[list_genes, ]
    freq_mat <- freq_mat[list_genes, ]

    psdb_expr_mat <- t(scale(t(psdb_expr_mat)))
    psdb_expr_mat[psdb_expr_mat > 3] <- 3
    psdb_expr_mat[psdb_expr_mat < -3] <- -3

    exp_df <- melt(psdb_expr_mat)
    freq_df <- melt(freq_mat)

    colnames(exp_df) <- c("gene", "cluster", "avg_expression")
    exp_df$freq <- freq_df$value

    head(exp_df)
    exp_df$cluster <- factor(exp_df$cluster, levels = stringr::str_sort(unique(exp_df$cluster), numeric = TRUE))
    exp_df$gene <- factor(exp_df$gene, levels = list_genes)
    pdf(file.path(figure_folder, paste0(prefix, "_bubbleplot"), "bubbleplot_specific_scaled.pdf"), width = 5, height = 5 + nlevels(exp_df$gene) * 0.1)
    print(ggplot(exp_df, aes(x = cluster, y = gene, colour = avg_expression, size = freq)) +
        geom_point() +
        theme_classic() +
        scale_colour_gradient2(breaks = c(-3, 0, 3), low = "#4575b4", mid =  "#ffffff", high = "#d73027") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
    )
    dev.off()
}

#### markers crossplot
target_cl1 <- "dmso_cluster_8"
target_cl2 <- "abt_cluster_9"

subso <- pms[, pms$clusters_from_samples %in% c(target_cl1, target_cl2)]
Idents(subso) <- ifelse(subso$clusters_from_samples == target_cl1, "cluster_1", "cluster_2")
mk <- FindMarkers(
    subso,
    ident.1 = "cluster_1",
    ident.2 = "cluster_2",
    logfc.threshold = 0.25,
    min.pct = 0.01,
    test.use = "wilcox_limma",
    only.pos = FALSE
)

actual_markers <- mk %>% dplyr::filter(p_val_adj < 0.05) %>% rownames
mk$gene_name <- rownames(mk)
mk$gene_name[!(mk$gene_name %in% actual_markers)] <- ""

pdf(file.path(figure_folder, "marker_crossplot.pdf"), width = 10, height = 10)
ggplot(mk, aes(x = pct.1, y = pct.2, size = -log10(p_val_adj), colour = avg_log2FC, label = gene_name)) +
    theme_classic() +
    geom_point() +
    ggrepel::geom_label_repel(force_pull = 0.1, size = 3.8, fill = "white", colour = "black") +
    scale_colour_gradient2(breaks = c(-3, 0, 3), low = "#4575b4", mid =  "#ffffff", high = "#d73027")
dev.off()

mk$gene_name <- rownames(mk)
mk$gene_name[!(mk$gene_name %in% list_genes)] <- ""

pdf(file.path(figure_folder, "marker_crossplot_given_list.pdf"), width = 10, height = 10)
ggplot(mk, aes(x = pct.1, y = pct.2, size = -log10(p_val_adj), colour = avg_log2FC, label = gene_name)) +
    theme_classic() +
    geom_point() +
    ggrepel::geom_label_repel(force_pull = 0.1, size = 3.8, fill = "white", colour = "black", nudge_x = 0.01, nudge_y = 0.01, max.overlaps = Inf) +
    scale_colour_gradient2(breaks = c(-3, 0, 3), low = "#4575b4", mid =  "#ffffff", high = "#d73027")
dev.off()


#### k stab plots ####
pdf(file.path(figure_folder, "k_stab_healthy.pdf"), width = 17, height = 7)
plot_k_n_partitions(ctrl_ca$"Highly_Variable"$"3000"$"clustering_stability", y_step = 50, dodge_width = 0.5) +
    scale_fill_viridis_c(limits = c(0.7, 1)) +
    scale_size_continuous(limits = c(0, 0.1)) +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)
    )
dev.off()

pdf(file.path(figure_folder, "k_stab_pms_dmso.pdf"), width = 17, height = 7)
plot_k_n_partitions(dmso_pms_ca$"Highly_Variable"$"2500"$"clustering_stability", y_step = 50, dodge_width = 0.5) +
    scale_fill_viridis_c(limits = c(0.7, 1)) +
    scale_size_continuous(limits = c(0, 0.2)) +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)
    )
dev.off()

pdf(file.path(figure_folder, "k_stab_pms_abt.pdf"), width = 17, height = 7)
plot_k_n_partitions(abt_pms_ca$"Highly_Variable"$"2000"$"clustering_stability", y_step = 50, dodge_width = 0.5) +
    scale_fill_viridis_c(limits = c(0.7, 1)) +
    scale_size_continuous(limits = c(0, 0.2)) +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)
    )
dev.off()
