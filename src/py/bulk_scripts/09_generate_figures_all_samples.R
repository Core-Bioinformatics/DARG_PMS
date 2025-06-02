library(reticulate)
library(ggplot2)
library(Seurat)
library(bulkAnalyseR)
library(ComplexHeatmap)
library(dplyr)
library(plotly)
library(MASS)
library(rgl)
library(processx)

source_folder <- 
expr_folder <- file.path(source_folder, "07.Expr_matrix")
bulk_folder <- file.path(source_folder, "08.Bulkanalyser")
output_folder <- file.path(source_folder, "09.Figures_all_samples")

if (!dir.exists(output_folder)) {
    dir.create(output_folder)
}

load(file.path(bulk_folder, "alex_feb_protein_coding", "expression_matrix.rda"))
load(file.path(bulk_folder, "alex_feb_protein_coding", "metadata.rda"))

##### PCA plot #####
generate_ellipsoid <- function(data, n_points = 200) {
  mu <- colMeans(data[, 1:3])    # Mean of x, y, z
  sigma <- cov.trob(data[, 1:3])$cov  # Robust covariance matrix
  
  eig <- eigen(sigma)  # Eigen decomposition
  values <- sqrt(eig$values)  # Radii of the ellipsoid
  vectors <- eig$vectors  # Principal axes
  
  # Create sphere points
  u <- seq(0, pi, length.out = n_points)
  v <- seq(0, 2 * pi, length.out = n_points)
  u_grid <- rep(u, each = n_points)
  v_grid <- rep(v, times = n_points)
  
  x_sphere <- cos(v_grid) * sin(u_grid)
  y_sphere <- sin(v_grid) * sin(u_grid)
  z_sphere <- cos(u_grid)
  
  sphere <- cbind(x_sphere, y_sphere, z_sphere)
  
  # Transform sphere into an ellipsoid
  ellipsoid <- sphere %*% diag(values) %*% t(vectors)
  ellipsoid <- sweep(ellipsoid, 2, mu, "+")  # Adjust position
  
  return(as.data.frame(ellipsoid))
}

metadata <- metadata[[1]]
expr_matrix <- expression.matrix[[1]]
metadata$"Sample.type" <- "PMS"
metadata$"Sample.type"[1:3] <- "Ctrl"
metadata$Sample.condition <- as.character(metadata$Sample.condition)
metadata$"Sample.condition"[1:3] <- paste0(metadata$Sample.condition[1:3], " Ctrl")
metadata$"Sample.condition"[4:6] <- paste0(metadata$Sample.condition[1:3], " PMS")

expr_matrix <- t(expr_matrix)

identical_genes <- apply(expr_matrix, 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))
expr_matrix <- expr_matrix[, !identical_genes]
most_abundant_genes <- colnames(expr_matrix)[order(colMeans(expr_matrix, na.rm = TRUE), decreasing = TRUE)]

pca_emb <- prcomp(expr_matrix[, most_abundant_genes[seq_len(500)]], scale = TRUE, center = TRUE)
pca_emb2 <- data.frame(pca_emb$x[, 1:3])
pca_emb2$sample_names <- rownames(pca_emb2)
pca_emb2$sample_types <- metadata$"Sample.condition"
# pca_emb2$sample_types <- factor(pca_emb2$sample_types, levels = c("Ctrl", "PMS"))


p1 <- ggplot(pca_emb2, aes(x = PC1, y = PC2, colour = .data$sample_types)) +
    geom_point(size = 5) +
    # scale_fill_manual(values = setNames(c("#4289c8", "#EF7E8F"), c("PMS", "Ctrl"))) +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18)
    ) +
    labs(
        x = paste0("PC1 (", summary(pca_emb)$importance[2, 1] * 100, "%)"),
        y = paste0("PC2 (", summary(pca_emb)$importance[2, 2] * 100, "%)")
    ) +
    ggforce::geom_mark_ellipse(
        aes(fill = .data$sample_types, colour = .data$sample_types),
        show.legend = FALSE
    )
p1

pdf(file.path(output_folder, "pca_all_by_type.pdf"), height = 15, width = 15)
print(p1)
dev.off()

scope <- kaleido()
fig <- plot_ly() %>% add_trace(
    data =pca_emb2,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~sample_types,
    # colors = c("#4289c8", "#EF7E8F"),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 25)
) %>% layout(
    scene = list(
        camera = list(eye = list(x = -2, y = 0.3, z = 1.1))
    )
)
scope$transform(fig, file.path(output_folder, "pac_3d_all_by_type.pdf"), width = 1500, height = 1500)

sub_expr_matrix <- expr_matrix[seq(from = 7, to = nrow(expr_matrix)), ]
most_abundant_genes <- colnames(sub_expr_matrix)[order(colMeans(sub_expr_matrix, na.rm = TRUE), decreasing = TRUE)]

pca_emb <- prcomp(sub_expr_matrix[, most_abundant_genes[seq_len(500)]], scale = TRUE, center = TRUE)
pca_emb2 <- data.frame(pca_emb$x[, 1:3])
pca_emb2$sample_types <- metadata$"Sample.condition"[seq(from = 7, to = nrow(metadata))]
pca_emb2$sample_types <- factor(pca_emb2$sample_types)

p2 <- ggplot(pca_emb2, aes(x = PC1, y = PC2, colour = .data$sample_types)) +
    geom_point(size = 5) +
    scale_fill_manual(values = setNames(c
("#4289c8", "#EF7E8F"), rev(levels(pca_emb2$sample_types)))) +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18)
    ) +
    labs(
        x = paste0("PC1 (", summary(pca_emb)$importance[2, 1] * 100, "%)"),
        y = paste0("PC2 (", summary(pca_emb)$importance[2, 2] * 100, "%)")
    ) +
    ggforce::geom_mark_ellipse(
        aes(fill = .data$sample_types, colour = .data$sample_types),
        show.legend = FALSE
    )
p2

pdf(file.path(output_folder, "pca_cm_pms_senloytic_by_type.pdf"), height = 15, width = 15)
print(p2)
dev.off()

ellipsoid <- pca_emb2 %>%
    group_by(sample_types) %>%
    group_split() %>%
    lapply(generate_ellipsoid) %>%
    bind_rows(.id = "sample_types")


ellipse <- ellipse3d(
    cov(pca_emb2[1:9, 1:3]),
    centre = colMeans(pca_emb2[1:9, 1:3]),
    level = 0.9999
)

plot_ly() %>% add_trace(
    data =pca_emb2,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~sample_types,
    colors = c("#4289c8", "#EF7E8F"),
    type = "scatter3d",
    mode = "markers"
) 

##### Jaccard heatmap #####
expression.matrix <- t(expr_matrix)
# p6 <- jaccard_heatmap(
#     expression.matrix = t(expr_matrix),
#     metadata = metadata,
#     n.abundant = 500,
#     show.values = FALSE
# )
# metadata$Sample_type <- c(rep("Ctrl", 3), rep("PMS", 3))
# metadata$id <- c(paste0("C", seq_len(3)), paste0("P", seq_len(3)))
# top_ha <- ComplexHeatmap::HeatmapAnnotation(
#     sample_type = metadata$Sample_type,
#     sample_name = metadata$id,
#     col = list(
#         sample_type = setNames(c("#4289c8", "#ef7a7c"), c("Ctrl", "PMS")),
#         sample_name = setNames(c("#2EAC63", "#95C11F", "#677935", "#F7A41B", "#E51B30", "#E9551F"), unique(metadata$id))
#     )
# )

# p6 + columnAnnotation(
#     sample_type = metadata$Sample_type
# )

# pdf(file.path(output_folder, "6_jaccard_heatmap.pdf"), width = 7, height = 5.5)
# print(top_ha %v% p6)
# dev.off()

##### Enrichment analysis #####
library(gprofiler2)
library(org.Hs.eg.db)
library(dplyr)
all_gene_names <- mapIds(org.Hs.eg.db, keys = colnames(expr_matrix), column = "SYMBOL", keytype = "ENSEMBL")

comparison_groups <- list(
    "pms_cm_vs_pms_senolytic_cm" = list(
        first = paste0("C", 7:15),
        second = paste0("C", 16:24)
    ),
    "base_ctrl_insc_vs_pms_cm" = list(
        first = paste0("C", 1:3),
        second = paste0("C", 7:15)
    ),
    "base_ctrl_insc_vs_pms_senolytic_cm" = list(
        first = paste0("C", 1:3),
        second = paste0("C", 16:24)
    ),
    "base_psm_insc_vs_pms_cm" = list(
        first = paste0("C", 4:6),
        second = paste0("C", 7:15)
    ),
    "base_psm_insc_vs_pms_senolytic_cm" = list(
        first = paste0("C", 4:6),
        second = paste0("C", 16:24)
    )
)
logfc <- 0.5

anno <- AnnotationDbi::select(
   org.Hs.eg.db::org.Hs.eg.db,
   keys = rownames(expression.matrix),
   keytype = 'ENSEMBL',
   columns = 'SYMBOL'
 ) %>%
   distinct(ENSEMBL, .keep_all = TRUE) %>%
   mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

for (comp_group in names(comparison_groups)) {
    mtd_temp <- metadata
    mtd_temp$new_column <- "other"
    frst <- comparison_groups[[comp_group]]$first
    scnd <- comparison_groups[[comp_group]]$second

    mtd_temp[frst, "new_column"] <- "first"
    mtd_temp[scnd, "new_column"] <- "second"

    de_table <- bulkAnalyseR::DEanalysis_edger(
        expression.matrix = expression.matrix[, c(frst, scnd)],
        condition = mtd_temp[c(frst, scnd), "new_column"],
        var1 = "first",
        var2 = "second",
        anno = anno
    ) %>% filter(abs(.data$log2FC) > logfc & .data$pvalAdj < 0.05) %>% arrange(desc(abs(.data$log2FC)))

    write.csv(de_table, file.path(output_folder, paste0(comp_group, ".csv")))

    enriched_terms <- de_table %>% filter(log2FC < 0)
    enriched_terms <- gost(
        query = enriched_terms$gene_name,
        ordered_query = TRUE,
        organism = "hsapiens",
        correction_method = "fdr",
        custom_bg = all_gene_names,
        sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF"),
        evcodes = TRUE
    )$result

    enriched_terms$intersection <- sapply(enriched_terms$intersection, function(x) paste(x, collapse = ";"))
    enriched_terms$parents <- sapply(enriched_terms$parents, function(x) paste(x, collapse = ";"))
    enriched_terms$evidence_codes <- sapply(enriched_terms$evidence_codes, function(x) paste(x, collapse = ";"))

    write.csv(enriched_terms, file.path(output_folder, paste0(comp_group, "_enriched_terms.csv")))

    enriched_terms <- de_table %>% filter(log2FC > 0)
    enriched_terms <- gost(
        query = enriched_terms$gene_name,
        ordered_query = TRUE,
        organism = "hsapiens",
        correction_method = "fdr",
        custom_bg = all_gene_names,
        sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF"),
        evcodes = TRUE
    )$result

    enriched_terms$intersection <- sapply(enriched_terms$intersection, function(x) paste(x, collapse = ";"))
    enriched_terms$parents <- sapply(enriched_terms$parents, function(x) paste(x, collapse = ";"))
    enriched_terms$evidence_codes <- sapply(enriched_terms$evidence_codes, function(x) paste(x, collapse = ";"))

    write.csv(enriched_terms, file.path(output_folder, paste0(comp_group, "_depleted_terms.csv")))
}

##### Volcano Heatmap #####
markers <- read.csv(file.path(output_folder, "pms_cm_vs_pms_senolytic_cm.csv"), row.names = 1)

p <- volcano_plot(
    genes.de.results = markers,
    pval.threshold = 0.05,
    lfc.threshold = 1,
    add.labels.custom = TRUE,
    add.labels.auto = FALSE,
    genes.to.label = c("CDKN1A", "SERPINE1", "ISG15", "IFIT1")
) + theme_bw() + theme(aspect.ratio = 1)

pdf(file.path(output_folder, "volcano_plot_pms_cm_vs_pms_senolytic_cm.pdf"), width = 7, height = 7)
print(p)
dev.off()

##### Gene heatmap ####
load(file.path(bulk_folder, "alex_feb_protein_coding", "expression_matrix.rda"))
expr_matrix <- expression.matrix[[1]]
gene_families <- list(
    "Senescence" = c("CDKN1A", "SERPINE1", "B2M", "PLK2", "BCL6", "PNPT1", "AKT3", "GDF15", "GADD45A", "RRAS", "TP53INP1", "TGFB1"),
    "Inflammation" = c("ISG15", "SPP1", "STAT1", "TNC", "CLU", "ANXA1", "HK1", "IFIT1", "IFI6", "IFITM2", "CD44", "TNFSF12", "IL21R", "TRIM25", "IFI16"),
    "Lipids" = c("DKK3", "TRIB3", "EGR1", "INSIG1", "STARD4", "CCN1", "ACSL3", "FABP3", "ABCA1", "CAV1", "ANXA2", "PLSCR1", "ABCB1", "IRS2", "HMGCS1", "SREBF1")
)
all_genes <- unique(unlist(gene_families))
anno <- AnnotationDbi::select(
   org.Hs.eg.db::org.Hs.eg.db,
   keys = rownames(expr_matrix),
   keytype = 'ENSEMBL',
   columns = 'SYMBOL'
 ) %>%
   distinct(ENSEMBL, .keep_all = TRUE) %>%
   mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

all_ens <- sapply(all_genes, function(x) {
    if (x %in% anno$SYMBOL) {
        matches <- anno[anno$SYMBOL == x, ]
        matches <- matches[!is.na(matches$SYMBOL), ]
        print(matches)

        return(matches$ENSEMBL[1])
    } else {
        return(NA)
    }
})
all_ens <- all_ens[!is.na(all_ens)]

for (gf in names(gene_families)) {
    gene_families[[gf]] <- intersect(gene_families[[gf]], names(all_ens))
}

heatmat <- expr_matrix[all_ens, ]
heatmat <- t(scale(t(heatmat)))
heatmat[heatmat > 3] <- 3
heatmat[heatmat < -3] <- -3
rownames(heatmat) <- names(all_ens)

metadata$Sample.condition.short <- c(
    rep("Ctrl iNSC", 3),
    rep("PMS iNSC", 3),
    rep("Ctrl iNSC and PMS CM", 9),
    rep("Ctrl iNSC and PMS Senolytic CM", 9)
)
metadata$Sample.condition.short <- factor(
    metadata$Sample.condition.short,
    levels = c("Ctrl iNSC", "PMS iNSC", "Ctrl iNSC and PMS CM", "Ctrl iNSC and PMS Senolytic CM")
)

htmp <- Heatmap(
    heatmat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(unlist(sapply(names(gene_families), function(x) { rep(x, length(gene_families[[x]]))})), levels = names(gene_families)),
    column_split = metadata$Sample.condition.short
)

pdf(file.path(output_folder, "all_gene_heatmap.pdf"), width = 15, height = 7)
print(htmp)
dev.off()

