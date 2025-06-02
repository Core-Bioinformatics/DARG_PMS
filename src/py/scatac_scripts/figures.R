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
figure_folder <- file.path(project_folder, "output", "figures")

if (!dir.exists(figure_folder)) {
  dir.create(figure_folder, recursive = TRUE)
}

so <- qs::qread(file.path(so_folder, "aggregate_norm.qs"), nthreads = ncores)

cscheme <- c("#f9af49", "#768b57", "#fe3619", "#01b876")
pdf(file.path(figure_folder, "umap_sample_name.pdf"), width = 8, height = 8)
DimPlot(so, group.by = "sample_name") +
    scale_color_manual(values = cscheme) +
    theme(
        aspect.ratio = 1
    ) +
    ggtitle("") 
dev.off()

cscheme <- c("#A87B95", "#59BFC7")
pdf(file.path(figure_folder, "umap_sample_type.pdf"), width = 8, height = 8)
DimPlot(so, group.by = "sample_type") +
    scale_color_manual(values = cscheme) +
    theme(
        aspect.ratio = 1
    ) +
    ggtitle("")
dev.off()

pdf(file.path(figure_folder, "umap_cluster.pdf"), width = 8, height = 8)
DimPlot(so, group.by = "stable_clusters_9") +
    theme(
        aspect.ratio = 1
    ) +
    ggtitle("")
dev.off()

pdf(file.path(figure_folder, "umap_ecc.pdf"), width = 8, height = 8)
FeaturePlot(so, features = "ecc_9") +
    scale_colour_viridis_c() +
    theme(
        aspect.ratio = 1
    ) +
    ggtitle("")
dev.off()


ca <- qs::qread(file.path(ca_folder, "aggregate.qs"), nthreads = ncores)
pdf(file.path(figure_folder, "k_stab.pdf"), width = 17, height = 7)
plot_k_n_partitions(ca$Highly_Variable_peaks$"55844"$clustering_stability, y_step = 50, dodge_width = 0.5) +
    scale_size_continuous(limits = c(0, 0.1)) +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)
    )
dev.off()

cscheme <- c("#A87B95", "#59BFC7")
cont_table <- table(so$sample_type, so$stable_clusters_9)

cont_table

# stacked histogram 
cont_table <- cbind(rowSums(cont_table[, c("1", "2", "4")]), cont_table)
colnames(cont_table)[1] <- "1_2_4"
cont_table <- cont_table[, setdiff(colnames(cont_table), c("1", "2", "4"))]

cont_table_df <- reshape2::melt(cont_table)
cont_table_df$Var2 <- factor(cont_table_df$Var2, levels = colnames(cont_table))
pdf(file.path(figure_folder, "barplot_cluster_sample_type_perc.pdf"), width = 8, height = 8)
ggplot(cont_table_df, aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = cscheme) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()
    ) +
    ggtitle("") +
    ylab("Proportion") +
    xlab("Clusters")
dev.off()

# stacked histogram with n cells on y axis
pdf(file.path(figure_folder, "barplot_cluster_sample_type_counts.pdf"), width = 8, height = 8)
ggplot(cont_table_df, aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cscheme) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()
    ) +
    ggtitle("") +
    ylab("Number of cells") +
    xlab("Clusters")
dev.off()
