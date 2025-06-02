setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
list2env(rjson::fromJSON(file = "00_configs.json"), envir = .GlobalEnv)
source(file.path(general_scripts_folder, "RNA", "create_seurat_from_cr_h5.R"))
library(ggplot2)
library(Seurat)
library(qs)
library(tidyverse)
library(Signac)
library(GenomicRanges)
library(future)

plan("multicore", workers = ncores)
options(future.globals.maxSize = 50 * 1024^3)

counts_folder <- file.path(project_folder, "preprocessing", "cellranger")
objects_folder <- file.path(project_folder, "objects", "R", "seurat")
if (!dir.exists(objects_folder)) {
    dir.create(objects_folder, recursive = TRUE)
}
metadata_folder <- file.path(project_folder, "metadata")

sample_metadata <- read.csv(file.path(metadata_folder, "samples.csv"), header = TRUE)

n_samples <- nrow(sample_metadata)

# peaks
grange_list <- lapply(seq_len(n_samples), function(i) {
    peak_path <- file.path(counts_folder, sample_metadata$id[i], "outs", "raw_peak_bc_matrix", "peaks.bed")
    peaks <- read.table(
        file = peak_path,
        col.names = c("chr", "start", "end")
    )

    makeGRangesFromDataFrame(peaks)
})

combined.peaks <- reduce(x = do.call(c, grange_list))

peak_widths <- width(combined.peaks)
combined.peaks <- combined.peaks[peak_widths > 20 & peak_widths < 10000]

# fragments
frags_list <- lapply(seq_len(n_samples), function(i) {
    sc_path <- file.path(counts_folder, sample_metadata$id[i], "outs", "singlecell.csv")
    fg_path <- file.path(counts_folder, sample_metadata$id[i], "outs", "fragments.tsv.gz")

    mtd <- read.table(
        file = sc_path,
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ]
    mtd <- mtd[mtd$passed_filters > 500, ]

    CreateFragmentObject(
        path = fg_path,
        cells = rownames(mtd)
    )
})

# seurat object
so_list <- lapply(seq_len(n_samples), function(i) {
    sc_path <- file.path(counts_folder, sample_metadata$id[i], "outs", "singlecell.csv")
    mtd <- read.table(
        file = sc_path,
        stringsAsFactors = FALSE,
        sep = ",",
        header = TRUE,
        row.names = 1
    )[-1, ]
    mtd <- mtd[mtd$passed_filters > 500, ]

    count_mat <- FeatureMatrix(
        fragments = frags_list[[i]],
        features = combined.peaks,
        cells = rownames(mtd)
    )

    so <- CreateSeuratObject(
        CreateChromatinAssay(count_mat, fragments = frags_list[[i]]),
        assay = "ATAC",
        meta.data = mtd
    )
    so$sample_name <- sample_metadata$id[i]
    so$sample_type <- sample_metadata$condition[i]
    return(so)
})

so_merged <- merge(
    x = so_list[[1]],
    y = so_list[seq(from = 2, to = n_samples)],
    add.cell.ids = sample_metadata$id
)

so_merged$sample_name <- factor(so_merged$sample_name)
so_merged$sample_type <- factor(so_merged$sample_type)
Idents(so_merged) <- "sample_name"

so_merged <- NucleosomeSignal(so_merged)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(so_merged) <- annotations

so_merged <- TSSEnrichment(so_merged)
so_merged$pct_reads_in_peaks <- so_merged$peak_region_fragments / so_merged$passed_filters * 100
so_merged$blacklist_ratio <- FractionCountsInRegion(
    object = so_merged,
    assay = 'ATAC',
    regions = blacklist_hg38_unified
)

mtd_list <- c("nFeature_ATAC", "nCount_ATAC", "nucleosome_signal", "TSS.enrichment", "blacklist_ratio")
qc_folder <- file.path(project_folder, "preprocessing", "qc", "ATAC")
if (!dir.exists(qc_folder)) {
    dir.create(qc_folder, recursive = TRUE)
}
pdf(file.path(qc_folder, "atac_violin_before_filtering.pdf"), width = 10, height = 10)
print(VlnPlot(so_merged, features = mtd_list, ncol=length(mtd_list), pt.size=0, log = TRUE, raster = TRUE))
dev.off()


# soft filtering
filtered_so <- subset(so_merged, nCount_ATAC > 1e3 & nCount_ATAC < 4e4 & nFeature_ATAC > 1e3 & nFeature_ATAC < 3e4 & TSS.enrichment > 3 & blacklist_ratio < 0.0025)
pdf(file.path(qc_folder, "atac_violin_intermediate_filtering.pdf"), width = 12, height = 7)
print(VlnPlot(filtered_so, features = mtd_list, ncol=length(mtd_list), pt.size=0, log = TRUE, raster = TRUE))
print(VlnPlot(filtered_so, features = mtd_list, ncol=length(mtd_list), pt.size=0, log = FALSE, raster = TRUE))
dev.off()
print(summary(filtered_so$sample_name))
qs::qsave(filtered_so, file.path(objects_folder, "intermediat_aggregate.qs"), nthreads = ncores)

# final filtering
final_so <- subset(filtered_so, nucleosome_signal < 0.75 & nCount_ATAC > 3e3 & nFeature_ATAC > 3e3)
pdf(file.path(qc_folder, "atac_violin_final_filtering.pdf"), width = 12, height = 7)
print(VlnPlot(final_so, features = mtd_list, ncol=length(mtd_list), pt.size=0, log = TRUE, raster = TRUE))
print(VlnPlot(final_so, features = mtd_list, ncol=length(mtd_list), pt.size=0, log = FALSE, raster = TRUE))
dev.off()

summary(final_so$sample_name)

# saturation plots
pdf(file.path(qc_folder, "atac_saturation_plots.pdf"), width = 10, height = 10)
for (x in seq_len(length(mtd_list) - 1)) {
    for (y in (x+1):length(mtd_list)) {
        mtd1 <- mtd_list[x]
        mtd2 <- mtd_list[y]
        df <- final_so@meta.data[, c(mtd1, mtd2, "sample_name")]

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

peak_path <- file.path(project_folder, "output", "peaks")
if (!dir.exists(peak_path)) {
    dir.create(peak_path, recursive = TRUE)
}
peak_list <- rownames(final_so)
peak_list <- strsplit(peak_list, "-")

peak_df <- data.frame(
    chr = sapply(peak_list, function(x) x[1]),
    start = as.numeric(sapply(peak_list, function(x) x[2])),
    end = as.numeric(sapply(peak_list, function(x) x[3]))
)

peak_df$peak_id <- paste0("peak_", seq_len(nrow(peak_df)))
peak_df$empty <- ""
peak_df$strand <- "+"

write.table(
    peak_df,
    file = file.path(peak_path, "aggregate_peaks.bed"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
)

so_atac <- readRDS("snATAC-seq-integrated-filter-Aug-2-2023.rds")
peak_list <- rownames(so_atac)
peak_list <- strsplit(peak_list, "-")

peak_df <- data.frame(
    chr = sapply(peak_list, function(x) x[1]),
    start = as.numeric(sapply(peak_list, function(x) x[2])),
    end = as.numeric(sapply(peak_list, function(x) x[3]))
)

peak_df$peak_id <- paste0("peak_", seq_len(nrow(peak_df)))
peak_df$empty <- ""
peak_df$strand <- "+"

write.table(
    peak_df,
    file = file.path(peak_path, "old_aggregate_peaks.bed"),
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
)
