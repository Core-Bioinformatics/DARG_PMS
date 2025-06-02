setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(foreach)
library(ggpubr)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
qc_folder <- file.path(project_folder, "preprocessing", "qc", "RNA")
options(future.globals.maxSize = 2000 * 1024^3)

suffix_options <- c("")
for (suffix in suffix_options) {
    print(suffix)

    so <- qread(file.path(objects_folder, paste0("aggregate", suffix, ".qs")), nthreads = ncores)

    n_median_genes <- floor(median(so$nFeature_RNA) / 500) * 500
    cc_genes <- intersect(unlist(cc.genes.updated.2019), rownames(so))
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = n_median_genes + length(cc_genes), verbose = TRUE)
    so <- ScaleData(so, features = rownames(so), verbose = TRUE)

    var_features <- setdiff(VariableFeatures(so), cc_genes)[seq_len(n_median_genes)]
    VariableFeatures(so) <- var_features
    print(length(VariableFeatures(so)))
    so <- CellCycleScoring(so, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)

    RhpcBLASctl::blas_set_num_threads(ncores)
    so <- RunPCA(so, npcs = 30, approx = TRUE, verbose = TRUE)
    so <- RunUMAP(so, reduction = "pca", dims = 1:30, verbose = TRUE, umap.method = "uwot")
    qsave(so, file.path(objects_folder, paste0("aggregate", suffix, "_norm.qs")), nthreads = ncores)

    pdf(file.path(qc_folder, paste0("aggregate", suffix, "_umaps.pdf")), width = 25, height = 10)
    discrete_plots <- DimPlot(so, group.by = c("sample_name", "sample_type", "treatment", "Phase"), combine = FALSE)
    cont_plots <- lapply(c("percent_mt", "percent_rp", "nCount_RNA", "nFeature_RNA"), function(feat) {
        FeaturePlot(so, features = feat, pt.size = 0.1, order = TRUE) + scale_colour_viridis_c()
    })
    all_plots <- c(discrete_plots, cont_plots)
    print(ggarrange(
        plotlist = all_plots,
        ncol = 4
    ))
    dev.off()

}

