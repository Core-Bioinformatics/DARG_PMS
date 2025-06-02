library(bulkAnalyseR)
library(ggplot2)
library(dplyr)

source_folder <- 
expr_matrix <- file.path(source_folder, "07.Expr_matrix/expr_matrix.csv")
metadata <- file.path(source_folder ,"metadata.csv")
output_dir <- file.path(source_folder, "08.Bulkanalyser")
gtf_path <- "GRch38_p113.gtf"


gtf <- read.table(gtf_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#") 
gtf <- gtf %>% filter(V3 == "gene")
gtf$gene_id <- sapply(regmatches(gtf$V9, regexpr("gene_id [a-zA-Z0-9_\\-]+;", gtf$V9)), function(x) gsub("gene_id ", "", gsub(";", "", x)))
gtf$gene_type <- sapply(regmatches(gtf$V9, regexpr("gene_biotype [a-zA-Z0-9_\\-]+;", gtf$V9)), function(x) gsub("gene_biotype ", "", gsub(";", "", x)))

ens_ids <- list(
  "all_genes" = gtf$gene_id,
  "all_genes_chr_1_22" = gtf %>% filter(V1 %in% paste0("", 1:22)) %>% pull(gene_id),
  "protein_coding" = gtf %>% filter(gene_type == "protein_coding") %>% pull(gene_id)
)

ens_ids[["protein_coding_chr_1_22"]] <- intersect(ens_ids[["protein_coding"]], ens_ids[["all_genes_chr_1_22"]])

if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

expr.matrix = read.csv(file = expr_matrix, row.names = 1)

metadata = read.csv(file = metadata, comment.char = "#")

app.title = "BulkAnalyseR App for Alex data - February 2025 samples"
organism = "hsapiens"
organism.db = "org.Hs.eg.db"


for(i in 2:ncol(metadata)) {
    metadata[,i] <- as.factor(metadata[, i])  
}
rownames(metadata) <- metadata$id

metadata <- metadata %>% filter(id %in% paste0("C", seq_len(6)))

expr.matrix <- expr.matrix[, metadata$id]
prefix <- "alex_feb_first_6_"

for (gene_types in names(ens_ids)) {
  shiny_dir <- file.path(output_dir, paste0(prefix, gene_types))
  if (file.exists(file.path(shiny_dir, "app.R")) && file.size(file.path(shiny_dir, "app.R")) > 0) {
    next
  }
  print(gene_types)
  sub_matrix <- expr.matrix[ens_ids[[gene_types]],]
  sub_matrix <- preprocessExpressionMatrix(sub_matrix, output.plot = TRUE)
  ggsave(file.path(output_dir, paste0("plot_", gene_types, ".pdf"))) + scale_color_brewer()
  write.table(sub_matrix, file.path(dirname(expr_matrix), paste0('normalised_expression_', gene_types, '.csv')))

  saveRDS(sub_matrix, file.path(dirname(expr_matrix), paste0("matrix_processed_", gene_types, ".rds")))

  generateShinyApp(
    shiny.dir = file.path(output_dir, paste0(prefix, gene_types)),
    app.title = app.title,
    modality = "RNA",
    expression.matrix = sub_matrix,
    metadata = metadata,#[1:22,],
    organism = organism,
    org.db = organism.db 
  )



}

