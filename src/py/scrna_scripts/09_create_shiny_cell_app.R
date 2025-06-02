setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(qs)
library(ShinyCell)

objects_folder <- file.path(project_folder, "objects", "R", "seurat")
ca_folder <- file.path(project_folder, "objects", "R", "clustassess")
mtd_folder <- file.path(project_folder, "metadata")
shiny_folder <- file.path(project_folder, "shiny_apps", "shiny_cell")
if (!dir.exists(shiny_folder)) {
    dir.create(shiny_folder, recursive = TRUE)
}

get_abbv <- function(ftype, fsize) {
    words <- strsplit(ftype, "_")[[1]]
    words <- sapply(words, function(word) stringr::str_to_lower(substr(word, 1, 1)))

    return(paste0(paste(words, collapse = ""), "_", fsize))
}

mtd_configs <- read.csv(file.path(mtd_folder, "configurations.csv"), header = TRUE, comment.char = "#")

for (i in seq_len(nrow(mtd_configs))) {
    if (is.null(mtd_configs$ftype[i])) {
        next
    }

    id <- mtd_configs$id[i]
    print(id)

    if (file.exists(file.path(shiny_folder, id, "sc1gexpr.h5"))) {
        next
    }

    so <- qread(file.path(objects_folder, paste0(id, "_norm.qs")), nthreads = ncores)
    DefaultAssay(so) <- "RNA"

    ca <- qread(file.path(file.path(ca_folder, paste0(id, ".qs"))), nthreads = ncores)
    so@reductions$umap@cell.embeddings <- ca[[mtd_configs$ftype[[i]]]][[as.character(mtd_configs$fsize[[i]])]]$umap
    colnames(so@reductions$umap@cell.embeddings) <- c("umap_1", "umap_2")

    sov4 <- CreateAssayObject(
        data = GetAssayData(object = so, slot = "data")
    )

    so[["RNA_v4"]] <- sov4
    
    sc_conf <- createConfig(so)
    makeShinyApp(
        so,
        sc_conf,
        gene.mapping = TRUE,
        gex.assay = "RNA_v4",
        shiny.dir = file.path(shiny_folder, id),
        shiny.title = paste0(id)
    )
}
