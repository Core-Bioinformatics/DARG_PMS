setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
library(Seurat)
library(ClustAssess)
library(qs)

objects_folder <- file.path(project_folder, "objects", "R")
so_folder <- file.path(objects_folder, "seurat")
ca_folder <- file.path(objects_folder, "clustassess")

mtd_folder <- file.path(project_folder, "metadata")
shiny_folder <- file.path(project_folder, "shiny_apps", "clustassess")
if (!dir.exists(shiny_folder)) {
    dir.create(shiny_folder, recursive = TRUE)
}

mtd_configs <- read.csv(file.path(mtd_folder, "configurations.csv"), header = TRUE, comment.char = "#")

for (i in seq_len(nrow(mtd_configs))) {
    id <- mtd_configs$id[i]
    print(id)

    if (!file.exists(file.path(so_folder, paste0(id, "_norm.qs")))) {
        next
    }
    if (!file.exists(file.path(ca_folder, paste0(id, ".qs")))) {
        next
    }

    if (file.exists(file.path(shiny_folder, id, "app.R"))) {
        next
    }

    so <- qread(file.path(so_folder, paste0(id, "_norm.qs")), nthreads = ncores)
    DefaultAssay(so) <- "RNA"
    ca_object <- qread(file.path(ca_folder, paste0(id, ".qs")), nthreads = ncores)

    write_shiny_app(
        object = so,
        assay_name = "RNA",
        clustassess_object = ca_object,
        shiny_app_title = id,
        project_folder = file.path(shiny_folder, id)
    )

}
