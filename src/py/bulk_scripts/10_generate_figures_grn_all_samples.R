#library(magrittr)
library(tidyr)
library(data.table)
library(biomaRt)
library(org.Hs.eg.db)

source_folder <- "/servers/iss-corescratch/am3019/250220_SP_alex_bulk/X204SC25013974-Z01-F001"
expr_folder <- file.path(source_folder, "07.Expr_matrix")
bulk_folder <- file.path(source_folder, "08.Bulkanalyser")
output_folder <- file.path(source_folder, "09.Figures_all_samples")

load(file.path(bulk_folder, "alex_feb_protein_coding", "expression_matrix.rda"))
load(file.path(bulk_folder, "alex_feb_protein_coding", "metadata.rda"))

gene_list <- list(
    # "NSCs" = c("SOX2", "PAX6", "NES", "FABP7", "YBX1", "PTPRZ1"),
    # "Glial progenitors" = c("HES1", "SLC1A3", "VIM", "CKB", "MEG3", "PTN", "CLU", "MT1X", "CST3", "ID4"),
    # "Senescence" = c("CDKN2A", "CDKN1A", "CDKN1C", "CDK5", "TP53", "RB1", "ATM", "PARP1"),
    # "Inflammatory" = c("NFKB1", "IFNG", "CSF2", "IL1R1", "IL4", "CSF3", "TNF", "IL1B", "NLRP3", "NOS2", "B2M", "IL6", "TGFB1"),
    # "Interferon associated" = c("STAT1", "STAT6", "ISG15", "IRF1", "IRF2", "IRF4", "CGAS"),
    "Larger Category" = c("CDKN2A", "CDKN1A", "CDKN1C", "TP53", "NFKB1", "IFNG", "TNF", "B2M", "STAT1", "STAT6", "ISG15", "CGAS", "IRF1", "IRF2", "IRF4", "IL6")
    # "Senescence Inflammatory" = c("B2M", "CGAS", "TBX2", "PLK2", "NUAK1", "BCL6", "YPEL3", "CDKN1A", "BMAL1", "SRF", "IGF1R", "ZNF277", "NEK6", "CALR", "WRN", "SIRT1", "MAPK10", "SMC6", "CDKN2A", "CDKN1C", "CDK5", "TP53", "RB1", "ATM", "PARP1", "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3", "NFKB", "IFNG", "CSF2", "IL1R1", "IL4", "CSF3", "TNF", "IL1B", "NLRP3", "NOS2", "B2M", "IL6", "TGFB1", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CSF1", "CSF2", "HMGB1", "TIMP2", "SPP1", "SERPINE1", "SERPINE2", "MMP1", "MMP10", "MMP12"),
    # "Interferon" =  c("STAT1", "STAT6", "ISG15", "IRF1", "IRF2", "IRF4", "IFITM3", "OASL", "BST2", "USP18", "IFIT1", "IFIT2", "IFIT3", "CD44", "IRF7", "TRIM45", "HLA-DPB1", "GBP2", "TRIM22", "TRIM62", "TRIM5", "IFNGR1", "OAS1", "OAS3", "OAS2", "IFNG", "TRIM38", "IRF6", "MT2A", "RAF1", "IRF8", "SOCS1", "SOCS3", "IRF9"),
    # "panel E" = c("CDKN2A", "CDKN1A", "CDKN1C", "TP53", "NFKB", "STAT1", "ISG15", "IFNG", "B2M", "TNF", "IL1R1", "IRF1", "IL1B", "IL6", "TGFB1", "STAT6", "IRF2", "IRF4", "CGAS")
)
gene_mapping <- mapIds(org.Hs.eg.db, keys = unique(unlist(gene_list)), column = "ENSEMBL", keytype = "SYMBOL")
gene_mapping_list <- lapply(gene_list, function(x) { gene_mapping[x]} )
names(gene_mapping_list) <- names(gene_list)

gene_mapping_list$negative <- c("ENSG00000259527", "ENSG00000176715", "ENSG00000123352", "ENSG00000166169", "ENSG00000122042", "ENSG00000103275", "ENSG00000153944", "ENSG00000155093", "ENSG00000171984", "ENSG00000250033", "ENSG00000137824", "ENSG00000104044", "ENSG00000167772", "ENSG00000100271", "ENSG00000179168", "ENSG00000224470", "ENSG00000225595", "ENSG00000124155", "ENSG00000152484", "ENSG00000005238", "ENSG00000158457", "ENSG00000111215", "ENSG00000231887", "ENSG00000065882", "ENSG00000185453", "ENSG00000146909", "ENSG00000100473", "ENSG00000091073", "ENSG00000129691", "ENSG00000135525", "ENSG00000172785", "ENSG00000277399", "ENSG00000132740", "ENSG00000147234", "ENSG00000137409", "ENSG00000115241", "ENSG00000134371", "ENSG00000135069", "ENSG00000109771", "ENSG00000146067", "ENSG00000178449", "ENSG00000186994", "ENSG00000151883", "ENSG00000085465", "ENSG00000142449", "ENSG00000259826", "ENSG00000178860", "ENSG00000109920", "ENSG00000219395", "ENSG00000244040", "ENSG00000201183", "ENSG00000111832", "ENSG00000133250", "ENSG00000100994", "ENSG00000109572", "ENSG00000083857", "ENSG00000216819", "ENSG00000204899", "ENSG00000062582", "ENSG00000110851", "ENSG00000160783", "ENSG00000258884", "ENSG00000102384", "ENSG00000154832", "ENSG00000271824", "ENSG00000116120", "ENSG00000139112", "ENSG00000158169", "ENSG00000048028", "ENSG00000185650", "ENSG00000196335", "ENSG00000106133", "ENSG00000290831", "ENSG00000140948", "ENSG00000169174", "ENSG00000129353", "ENSG00000133612", "ENSG00000162222", "ENSG00000171574")
mapping <- mapIds(org.Hs.eg.db, keys = unlist(gene_mapping_list), column = "SYMBOL", keytype = "ENSEMBL")

get_link_list_rename <- function(weightMat, plotConnections){
  GENIE3::getLinkList(weightMat, plotConnections) %>%
    dplyr::mutate(from = as.character(.data$regulatoryGene), 
                  to = as.character(.data$targetGene), 
                  value = .data$weight, 
                  regulatoryGene = NULL, 
                  targetGene = NULL,
                  weight = NULL)
}

plot_GRN_multiomics <- function(weightMat, plotConnections){
  
  edges <- get_link_list_rename(weightMat, plotConnections)
  edges$id = paste0(pmin(edges$from,edges$to),pmax(edges$from,edges$to))
  edges = edges[order(-edges$value),]
  edges = edges[!duplicated(edges$id),]
  edges = edges[,1:3]
  edges$color = 'black'

  # add nodes for each end of the edges and colour them by modality
  nodes <- tibble::tibble(
    id = c(edges$to, edges$from),
    label = c(edges$to, edges$from),
    value=1,
    font.color='black',
    font.size=50,
    font.bold=T
  ) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)

  nodes$value = 100
  nodes = nodes %>% dplyr::filter(id %in% names(mapping))
  new_labels = c()
  for (label in nodes$label) {
      new_labels = append(new_labels, mapping[[label]])
  }
  nodes$label = new_labels
  

  edges = edges[edges$from != edges$to,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  edges = edges[edges$to %in% nodes$id,]
  edges = edges[edges$from %in% nodes$id,]
  nodes = nodes[nodes$id %in% c(edges$from,edges$to),]
  nodes$label = gsub('_','\n',nodes$label)
  nodes$label = gsub('/','/\n',nodes$label)
  nodes = nodes[nodes$id!='',]
  
  network <- visNetwork::visNetwork(nodes, edges,width=992,height=744) %>%
    visNetwork::visGroups(groupname = "NSCs", color = '#bada55', shape = "dot", size=50) %>%
    visNetwork::visGroups(groupname = "Inflammation", color = '#ff7373', shape = "dot", size=50) %>%
    visNetwork::visGroups(groupname = "Lipids", color = '#b0e0e6', shape = "dot", size=50) %>%
    visNetwork::visGroups(groupname = "Hub", color = '#BFE4F7', shape = "box") %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based",
                           forceAtlas2Based = list(gravitationalConstant = -150,springConstant=0.1,centralGravity=0.005)) %>%
    visNetwork::visLegend()%>%
    visNetwork::visLayout(randomSeed = 2023)
  nodes$color = '#8bc0e0'
  return(list('nodes'=nodes,'edges'=edges,'network'=network))
}

exp.mat <- expression.matrix[[1]]


conditions = list(
  'pms cm'=paste0("C", 7:15),
  'senolytic cm'=paste0("C", 16:24)
)

edges = c(50)

for (n_edges in edges) {
    for (condition in names(conditions)) {
        columns = conditions[[condition]]
        for (set in names(gene_mapping_list)) {
            rows <- gene_mapping_list[[set]]
            exp.mat.temp = exp.mat[rownames(exp.mat)%in%c(rows),]
            exp.mat.temp = exp.mat.temp[,colnames(exp.mat.temp)%in%columns]
            set.seed(2023)
            res <- GENIE3::GENIE3(as.matrix(exp.mat.temp), targets = rownames(exp.mat.temp),
                  nCores = 16)
            vis_result <- plot_GRN_multiomics(res, n_edges)$network
            visNetwork::visSave(
                vis_result,
                file.path(output_folder, paste(condition, '_', set, '_', n_edges, '.html', sep='')),
                selfcontained = F, background = "white"
            )
        } 
    }
}
