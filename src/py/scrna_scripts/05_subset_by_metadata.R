setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(ClustAssess)
library(qs)
library(ggplot2)

options(future.globals.maxSize = 2 * 1024^9)

metadata_folder <- file.path(project_folder, "metadata")
seurat_folder <- file.path(project_folder, "objects", "R", "seurat")

subs_list <- rjson::fromJSON(file = file.path(metadata_folder, "subsets.json"))

for (object_name in names(subs_list)) {
    print(object_name)
    seurat_path <- file.path(seurat_folder, paste0(object_name, "_norm.qs"))

    subset_paths <- file.path(seurat_folder, sapply(subs_list[[object_name]], function(sub) paste0(sub$name_subset, "_norm.qs")))

    check_existing <- sapply(subset_paths, function(subset_path) file.exists(subset_path) & file.size(subset_path) > 0)
    names(check_existing) <- sapply(subs_list[[object_name]], function(sub) sub$name_subset)
    if (all(check_existing)) {
        print("All subsets already exist.")
        next
    }

    so <- qread(seurat_path, nthreads = ncores)

    mtd_names <- colnames(so@meta.data)
    for (target_mtd_name in mtd_names[grepl("^ecc_|^stable_clusters_", mtd_names)]) {
        so@meta.data[[target_mtd_name]] <- NULL
    }

    if ("stable_clusters" %in% mtd_names) {
        so@meta.data[[paste0(object_name, "_stable_clusters")]] <- so$stable_clusters
        so$stable_clusters <- NULL
    }

    for (sub in subs_list[[object_name]]) {
        name_subset <- sub$name_subset
        if (file.exists(file.path(seurat_folder, paste0(name_subset, "_norm.qs"))) && file.size(file.path(seurat_folder, paste0(name_subset, "_norm.qs"))) > 0) {
            print(paste0(name_subset, " already exists. Skipping."))
            next
        }
        mtd_name <- sub$name_metadata
        mtd_values <- sub$values

        if (check_existing[name_subset]) {
            print(paste0(name_subset, " already exists. Skipping."))
            next
        }

        so_subset <- subset(so, cells = colnames(so)[as.character(so@meta.data[[mtd_name]]) %in% mtd_values])
        for (mtd_name in colnames(so_subset@meta.data)) {
            if (!is.factor(so_subset@meta.data[[mtd_name]])) {
                next
            }

            so_subset@meta.data[[mtd_name]] <- factor(as.character(so_subset@meta.data[[mtd_name]]))
        }

        n_median_genes <- median(so_subset$nFeature_RNA)
        n_median_genes <- floor(n_median_genes / 500) * 500

        g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(so_subset))
        s_genes <- intersect(cc.genes.updated.2019$s.genes, rownames(so_subset))
        so_subset <- NormalizeData(so_subset, normalization.method = "LogNormalize", scale.factor = 10000)
        so_subset <- FindVariableFeatures(so_subset, selection.method = "vst", nfeatures = n_median_genes + length(g2m_genes) + length(s_genes), verbose = TRUE)
        so_subset <- ScaleData(so_subset, features = rownames(so_subset), verbose = TRUE)

        var_features <- setdiff(VariableFeatures(so_subset), c(g2m_genes, s_genes))[seq_len(n_median_genes)]
        VariableFeatures(so_subset) <- var_features

        RhpcBLASctl::blas_set_num_threads(ncores)
        so_subset <- RunPCA(so_subset, npcs = 30, approx = TRUE, verbose = TRUE)
        so_subset <- RunUMAP(so_subset, reduction = "pca", dims = 1:30, verbose = TRUE, umap.method = "uwot")

        qsave(so_subset, file.path(seurat_folder, paste0(name_subset, "_norm.qs")), nthreads = ncores)
        print(DimPlot(so_subset, group.by = "sample_name") + ggtitle(name_subset))
    }
}
