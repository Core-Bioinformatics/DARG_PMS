setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(foreach)
library(ggpubr)
library(Signac)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
qc_folder <- file.path(project_folder, "preprocessing", "qc", "ATAC")
options(future.globals.maxSize = 2000 * 1024^3)

prefix_options <- c("", "intermediate_")
mtd_list <- c("nFeature_ATAC", "nCount_ATAC", "nucleosome_signal", "TSS.enrichment")
for (prefix in prefix_options) {
    print(prefix)

    so <- qread(file.path(objects_folder, paste0(prefix, "aggregate.qs")), nthreads = ncores)

    so <- RunTFIDF(
        so,
        scale.factor = 1e5
    )

    so <- FindTopFeatures(so)
    so <- RunSVD(so)
    so <- RunUMAP(so, reduction = "lsi", dims = 2:30, verbose = TRUE, umap.method = "uwot")
    
    qsave(so, file.path(objects_folder, paste0(prefix, "aggregate_norm.qs")), nthreads = ncores)

    pdf(file.path(qc_folder, paste0(prefix, "aggregate_umaps.pdf")), width = 25, height = 10)
    discrete_plots <- DimPlot(so, group.by = c("sample_name", "sample_type"), combine = FALSE)
    cont_plots <- lapply(mtd_list, function(feat) {
        FeaturePlot(so, features = feat, pt.size = 0.1, order = TRUE) + scale_colour_viridis_c()
    })
    all_plots <- c(discrete_plots, cont_plots)
    print(ggarrange(
        plotlist = all_plots,
        ncol = 4
    ))
    dev.off()

    so_list <- SplitObject(so, split.by = "sample_type")
    for (i in seq_along(so_list)) {
        so_list[[i]] <- RunTFIDF(
            so_list[[i]],
            scale.factor = 1e5
        )

        so_list[[i]] <- FindTopFeatures(so_list[[i]])
        so_list[[i]] <- RunSVD(so_list[[i]])
    }

    anchors <- FindIntegrationAnchors(
        object.list = so_list,
        anchor.features = rownames(so),
        reduction = "rlsi",
        dims = 2:30
    )

    integrate_embeddings <- IntegrateEmbeddings(
        anchorset = anchors,
        reductions = so@reductions$lsi,
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = 1:30https://www.youtube.com/watch?list=RDiRn9aN9D_c8&v=zntYcrksPCI
    )
    
    integrate_embeddings <- RunUMAP(integrate_embeddings, reduction = "integrated_lsi", dims = 2:30, verbose = TRUE, umap.method = "uwot")

    pdf(file.path(qc_folder, paste0(prefix, "aggregate_integrated_umaps.pdf")), width = 25, height = 10)
    discrete_plots <- DimPlot(integrate_embeddings, group.by = c("sample_name", "sample_type"), combine = FALSE)
    cont_plots <- lapply(mtd_list, function(feat) {
        FeaturePlot(integrate_embeddings, features = feat, pt.size = 0.1, order = TRUE) + scale_colour_viridis_c()
    })
    all_plots <- c(discrete_plots, cont_plots)
    print(ggarrange(
        plotlist = all_plots,
        ncol = 4
    ))
    dev.off()

}


