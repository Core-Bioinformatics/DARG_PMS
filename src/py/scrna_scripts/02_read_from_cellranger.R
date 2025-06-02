setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
source(file.path(general_scripts_folder, "RNA", "create_seurat_from_cr_h5.R"))
library(ggplot2)
library(Seurat)
library(qs)
library(tidyverse)

counts_folder <- file.path(project_folder, "preprocessing", "cellranger")
objects_folder <- file.path(project_folder, "objects", "R", "seurat")
if (!dir.exists(objects_folder)) {
    dir.create(objects_folder, recursive = TRUE)
}
metadata_folder <- file.path(project_folder, "metadata")

sample_metadata <- read.csv(file.path(metadata_folder, "2501_samples.csv"), header = TRUE)

n_samples <- nrow(sample_metadata)

so_list <- lapply(seq_len(n_samples), function(i) {
    h5_path <- file.path(counts_folder, sample_metadata$sample_name[i], "outs", "raw_feature_bc_matrix.h5")
    print(h5_path)
    so <- create_so_default(
        h5_path = h5_path,
        sample_name = sample_metadata$sample_name[i],
    )
    so$treatment <- factor(sample_metadata$treatment[i])
    so$sample_type <- factor(sample_metadata$type[i])
    so
})

so_merged <- JoinLayers(merge(
    x = so_list[[1]],
    y = so_list[seq(from = 2, to = n_samples)],
    add.cell.ids = sample_metadata$sample_name
))

so_merged$sample_name <- factor(so_merged$sample_name)
so_merged$treatment <- factor(so_merged$treatment)
so_merged$sample_type <- factor(so_merged$sample_type)
Idents(so_merged) <- "sample_name"

mt_regex <- "^MT-"
rp_regex <- "^RP[SL]\\d+"
mrp_regex <- "^MRP[SL]\\d+"
mt_genes <- grep(mt_regex, rownames(so_merged), value = FALSE)
rp_genes <- grep(rp_regex, rownames(so_merged), value = FALSE)
mrp_genes <- grep(mrp_regex, rownames(so_merged), value = FALSE)

so_merged$percent_mt <- PercentageFeatureSet(so_merged, features = rownames(so_merged)[mt_genes])
so_merged$percent_rp <- PercentageFeatureSet(so_merged, features = rownames(so_merged)[c(rp_genes, mrp_genes)])
so_merged$percent_mt[is.na(so_merged$percent_mt)] <- 100
so_merged$percent_rp[is.na(so_merged$percent_rp)] <- 100
so_merged <- subset(so_merged, features = -c(mt_genes, rp_genes, mrp_genes))

mtd_rna <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_rp")
qc_folder <- file.path(project_folder, "preprocessing", "qc", "RNA")
if (!dir.exists(qc_folder)) {
    dir.create(qc_folder, recursive = TRUE)
}
pdf(file.path(qc_folder, "rna_violin_before_filtering.pdf"))
print(VlnPlot(so_merged, features = mtd_rna, ncol=4, pt.size=0, log = TRUE, raster = TRUE))
dev.off()


# soft filtering
filtered_so <- subset(so_merged, nCount_RNA > 1e3 & nFeature_RNA > 1e3 & percent_mt < 50 & percent_rp < 50)
pdf(file.path(qc_folder, "rna_violin_intermediate_filtering.pdf"))
print(VlnPlot(filtered_so, features = mtd_rna, ncol=4, pt.size=0, log = TRUE, raster = TRUE))
print(VlnPlot(filtered_so, features = mtd_rna, ncol=4, pt.size=0, log = FALSE, raster = TRUE))
dev.off()
print(summary(filtered_so$sample_name))

# final filtering
final_so <- subset(filtered_so, nCount_RNA < 6e4 & nCount_RNA > 3e3 & percent_rp > 15 & percent_rp < 40 & percent_mt < 10 & percent_mt > 2)
pdf(file.path(qc_folder, "rna_violin_final_filtering.pdf"))
print(VlnPlot(final_so, features = mtd_rna, ncol=4, pt.size=0, log = TRUE, raster = TRUE))
print(VlnPlot(final_so, features = mtd_rna, ncol=4, pt.size=0, log = FALSE, raster = TRUE))
dev.off()

summary(final_so$sample_name)

# saturation plots
pdf(file.path(qc_folder, "rna_saturation_plots.pdf"), width = 10, height = 10)
for (x in seq_len(length(mtd_rna) - 1)) {
    for (y in (x+1):length(mtd_rna)) {
        mtd1 <- mtd_rna[x]
        mtd2 <- mtd_rna[y]
        df <- final_so@meta.data[, c(mtd1, mtd2, "sample_name")]

        print(ggplot(df, aes_string(x = mtd1, y = mtd2, color = "sample_name")) +
            geom_point() +
            theme_classic() 
        )

        for (lib_name in unique(df$sample_name)) {
            print(ggplot(df %>% filter(sample_name == lib_name), aes_string(x = mtd1, y = mtd2)) +
                geom_point() +
                theme_classic() +
                ggtitle(lib_name)
            )
        }
    }
}
dev.off()

qs::qsave(final_so, file.path(objects_folder, "aggregate.qs"), nthreads = ncores)
