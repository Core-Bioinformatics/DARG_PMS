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

peak_folder <- file.path(project_folder, "output", "peaks")

new_ann_path <- file.path(peak_folder, "annotation_aggregate_peaks.txt")
old_ann_path <- file.path(peak_folder, "annotation_old_aggregate_peaks.txt")

#### NEW SAMPLES ####
so <- qread(file.path(project_folder, "objects", "R", "seurat", "aggregate_norm.qs"), nthreads = ncores)
new_ann <- read.table(new_ann_path, header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
rownames(new_ann) <- paste0(new_ann$Chr, "-", new_ann$Start - 1, "-", new_ann$End)

peak_categs <- list(
    "intergenic" = new_ann %>% filter(Annotation == "Intergenic") %>% rownames,
    "genic_NC" = new_ann %>% filter(Annotation != "Intergenic", Gene.Type != "protein-coding") %>% rownames,
    "genic_PC" = new_ann %>% filter(Annotation != "Intergenic", Gene.Type == "protein-coding") %>% rownames
)

count_mat <- GetAssayData(so, layer = "counts")

uniq_samples <- unique(so$sample_name)
psdb_categs <- matrix(0, nrow = length(peak_categs), ncol = length(uniq_samples))
for (i in seq_along(peak_categs)) {
    for (j in seq_along(uniq_samples)) {
        psdb_categs[i, j] <- sum(count_mat[peak_categs[[i]], so$sample_name == uniq_samples[j]])
    }
}
rownames(psdb_categs) <- names(peak_categs)
colnames(psdb_categs) <- uniq_samples

# for (i in seq_len(ncol(psdb_categs))) {
#     psdb_categs[, i] <- psdb_categs[, i] / sum(psdb_categs[, i]) 
# }

new_psdb <- psdb_categs

#### OLD SAMPLES ####
so_rna <- readRDS("metadata-integrated-2023-scRNAseq-harmony-exclude.rds")
so_atac <- readRDS("snATAC-seq-integrated-filter-Aug-2-2023.rds")
mtd_path <- read.csv("metadata-integrated-harmony-2023-matched-with-ATAC-cells.csv")

mtd_path <- mtd_path %>% filter(barcode %in% colnames(so_rna)) %>% filter(matched_cells %in% colnames(so_atac)) %>% filter(!(matched_cells %in% colnames(so_atac)[so_atac$seurat_clusters == 9]))

so_rna <- so_rna[, mtd_path$barcode]
so_atac <- so_atac[, mtd_path$matched_cells]
old_ann <- read.table(old_ann_path, header = TRUE, comment.char = "#", sep = "\t", fill = TRUE, quote = "")
rownames(old_ann) <- paste0(old_ann$Chr, "-", old_ann$Start - 1, "-", old_ann$End)

peak_categs <- list(
    "intergenic" = old_ann %>% filter(Annotation == "Intergenic") %>% rownames,
    "genic_NC" = old_ann %>% filter(Annotation != "Intergenic", Gene.Type != "protein-coding") %>% rownames,
    "genic_PC" = old_ann %>% filter(Annotation != "Intergenic", Gene.Type == "protein-coding") %>% rownames
)

count_mat <- GetAssayData(so_atac, layer = "counts")

uniq_samples <- unique(so_atac$sample)
psdb_categs <- matrix(0, nrow = length(peak_categs), ncol = length(uniq_samples))
for (i in seq_along(peak_categs)) {
    for (j in seq_along(uniq_samples)) {
        psdb_categs[i, j] <- sum(count_mat[peak_categs[[i]], so_atac$sample == uniq_samples[j]])
    }
}
rownames(psdb_categs) <- names(peak_categs)
colnames(psdb_categs) <- uniq_samples

old_psdb <- psdb_categs


###### Stacked histogram proportion ######
library(reshape2)
library(ggplot2)
figure_path <- file.path(project_folder, "output", "figures")
if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}

current_df <- cbind(
    new_psdb,
    old_psdb
)
for (i in seq_len(ncol(current_df))) {
    current_df[, i] <- current_df[, i] / sum(current_df[, i]) * 100
}

current_df <- melt(current_df)
current_df

pdf(file.path(figure_path, "reads_stacked_histogram.pdf"))
print(ggplot(current_df, aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
        x = "Sample",
        y = "Percentage of fragments incident to peak regions",
        fill = "Region category"
    )
)
dev.off()

##### Statistic tests #####
library(confintr)
####### chi^2 #######
# pairwise sample comparison

chi2_mat <- matrix(NA, nrow = 8, ncol = 8)
rownames(chi2_mat) <- c(colnames(new_psdb), colnames(old_psdb))
colnames(chi2_mat) <- c(colnames(new_psdb), colnames(old_psdb))

current_df <- cbind(
    new_psdb,
    old_psdb
)
for (i in seq(from = 1, to = 7)) {
    for (j in seq(from = i + 1, to = 8)) {
        temp_df <- current_df[, c(i, j)]

        temp_df[, 1] <- temp_df[, 1] / sum(temp_df[, 1]) * 100
        temp_df[, 2] <- temp_df[, 2] / sum(temp_df[, 2]) * 100

        chi2_mat[i, j] <- chisq.test(temp_df)$p.value
        chi2_mat[j, i] <- chi2_mat[i, j]
    }
}

col_fun <- circlize::colorRamp2(c(0, 0.05, 1), c("red", "white", "blue"))
pdf(file.path(figure_path, "atac_reads_chi2_heatmap.pdf"), width = 10, height = 10)
Heatmap(
    chi2_mat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (is.na(chi2_mat[i, j])) {
            return(grid.text("", x, y, gp = gpar(fontsize = 10)))
        }

        if (chi2_mat[i, j] > 0.05) {
            return(grid.text(round(chi2_mat[i, j], 2), x, y, gp = gpar(fontsize = 10)))
        }

        if (chi2_mat[i, j] > 0.01) {
            return(grid.text("*", x, y, gp = gpar(fontsize = 10)))
        }

        if (chi2_mat[i, j] > 0.001) {
            return(grid.text("**", x, y, gp = gpar(fontsize = 10)))
        }

        if (chi2_mat[i, j] > 0.0001) {
            return(grid.text("***", x, y, gp = gpar(fontsize = 10)))
        }

        return(grid.text("****", x, y, gp = gpar(fontsize = 10)))
    },
    col = col_fun
)
dev.off()
