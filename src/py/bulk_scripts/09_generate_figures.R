library(ggplot2)
library(Seurat)
library(bulkAnalyseR)
library(ComplexHeatmap)
library(dplyr)

source_folder <-
expr_folder <- file.path(source_folder, "07.Expr_matrix")
bulk_folder <- file.path(source_folder, "08.Bulkanalyser")
output_folder <- file.path(source_folder, "09.Figures")

if (!dir.exists(output_folder)) {
    dir.create(output_folder)
}

load(file.path(bulk_folder, "alex_feb_first_6_protein_coding", "expression_matrix.rda"))
load(file.path(bulk_folder, "alex_feb_first_6_protein_coding", "metadata.rda"))
##### PCA plot #####
metadata <- metadata[[1]]
expr_matrix <- expression.matrix[[1]]
colnames(expr_matrix) <- c(paste0("C", seq_len(3)), paste0("P", seq_len(3)))
expr_matrix <- t(expr_matrix)

identical_genes <- apply(expr_matrix, 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))
expr_matrix <- expr_matrix[, !identical_genes]
most_abundant_genes <- colnames(expr_matrix)[order(colMeans(expr_matrix, na.rm = TRUE), decreasing = TRUE)]

pca_emb <- prcomp(expr_matrix[, most_abundant_genes[seq_len(500)]], scale = TRUE, center = TRUE)
pca_emb2 <- data.frame(pca_emb$x[, 1:2])
pca_emb2$sample_names <- rownames(pca_emb2)

p1 <- ggplot(pca_emb2, aes(x = PC1, y = PC2, fill = sample_names)) +
    geom_point(size = 10, colour = "black", shape = 21) +
    scale_fill_manual(values = c("#2EAC63", "#95C11F", "#677935", "#F7A41B", "#E51B30", "#E9551F")) +
    theme_bw() +
    theme(
        aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(
        x = paste0("PC1 (", summary(pca_emb)$importance[2, 1] * 100, "%)"),
        y = paste0("PC2 (", summary(pca_emb)$importance[2, 2] * 100, "%)")
    )
p1

pdf(file.path(output_folder, "1_pca_plot.pdf"))
print(p1)
dev.off()

##### Jaccard heatmap #####
p6 <- jaccard_heatmap(
    expression.matrix = t(expr_matrix),
    metadata = metadata,
    n.abundant = 500,
    show.values = FALSE
)
metadata$Sample_type <- c(rep("Ctrl", 3), rep("PMS", 3))
metadata$id <- c(paste0("C", seq_len(3)), paste0("P", seq_len(3)))
top_ha <- ComplexHeatmap::HeatmapAnnotation(
    sample_type = metadata$Sample_type,
    sample_name = metadata$id,
    col = list(
        sample_type = setNames(c("#4289c8", "#ef7a7c"), c("Ctrl", "PMS")),
        sample_name = setNames(c("#2EAC63", "#95C11F", "#677935", "#F7A41B", "#E51B30", "#E9551F"), unique(metadata$id))
    )
)

p6 + columnAnnotation(
    sample_type = metadata$Sample_type
)

pdf(file.path(output_folder, "6_jaccard_heatmap.pdf"), width = 7, height = 5.5)
print(top_ha %v% p6)
dev.off()

##### Enrichment analysis #####
library(gprofiler2)
library(org.Hs.eg.db)
library(dplyr)
all_gene_names <- mapIds(org.Hs.eg.db, keys = colnames(expr_matrix), column = "SYMBOL", keytype = "ENSEMBL")

markers <- read.csv(file.path(output_folder, "DEset.csv"))

enriched_terms <- markers %>% filter(log2FC < 0)
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

write.csv(enriched_terms, file.path(output_folder, "enriched_terms.csv"))

enriched_terms <- markers %>% filter(log2FC > 0)
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

write.csv(enriched_terms, file.path(output_folder, "depleted_terms.csv"))

##### Volcano Heatmap #####

p <- volcano_plot(
    genes.de.results = markers,
    pval.threshold = 0.05,
    lfc.threshold = 1,
    add.labels.custom = TRUE,
    add.labels.auto = FALSE,
    genes.to.label = c("ISG15", "FASN", "IFIT1", "SERPINE1", "NLRP2", "TGFB1", "NOTCH1")
) + theme_bw() + theme(aspect.ratio = 1)

pdf(file.path(output_folder, "volcano_plot.pdf"), width = 7, height = 7)
print(p)
dev.off()


##### Gene heatmap #####
load(file.path(bulk_folder, "alex_feb_first_6_protein_coding", "expression_matrix.rda"))
expr_matrix <- expression.matrix[[1]]
g1 <- c("IL13", "ANG", "IL7", "CXCL10", "CCL2", "TIMP2", "SPP1", "CSF3", "IL5", "TNF", "CSF2", "TIMP1", "TGFB2", "PGF", "TNFRSF1B", "IL4", "IGFBP2", "IGFBP3", "IGF1", "IGFBP4", "TGFB1", "CCL5", "PDGFB", "CCL13", "TGFB3", "CSF1", "CUL9", "IL1A", "IL2")
anno <- AnnotationDbi::select(
   org.Hs.eg.db::org.Hs.eg.db,
   keys = rownames(expr_matrix),
   keytype = 'ENSEMBL',
   columns = 'SYMBOL'
 ) %>%
   distinct(ENSEMBL, .keep_all = TRUE) %>%
   mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

all_ens <- sapply(g1, function(x) {
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
g1 <- intersect(g1, names(all_ens))

heatmat <- expr_matrix[all_ens, ]
heatmat <- t(scale(t(heatmat)))
heatmat[heatmat > 3] <- 3
heatmat[heatmat < -3] <- -3
rownames(heatmat) <- names(all_ens)

htmp <- Heatmap(
    heatmat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = c(rep("Ctrl", 3), rep("PMS", 3))
)

pdf(file.path(output_folder, "first_6_gene_heatmap1.pdf"), width = 7, height = 7)
print(htmp)
dev.off()

g2 <- c("ISG15", "HSPG2", "FN1", "CLU", "TNC", "CD44", "A2M", "GDF15", "IL4R", "B2M", "IL4R", "CXCL3")
anno <- AnnotationDbi::select(
   org.Hs.eg.db::org.Hs.eg.db,
   keys = rownames(expr_matrix),
   keytype = 'ENSEMBL',
   columns = 'SYMBOL'
 ) %>%
   distinct(ENSEMBL, .keep_all = TRUE) %>%
   mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

all_ens <- sapply(g2, function(x) {
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
g2 <- intersect(g2, names(all_ens))

heatmat <- expr_matrix[all_ens, ]
heatmat <- t(scale(t(heatmat)))
heatmat[heatmat > 3] <- 3
heatmat[heatmat < -3] <- -3
rownames(heatmat) <- names(all_ens)

htmp <- Heatmap(
    heatmat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = c(rep("Ctrl", 3), rep("PMS", 3))
)

pdf(file.path(output_folder, "first_6_gene_heatmap2.pdf"), width = 7, height = 7)
print(htmp)
dev.off()
