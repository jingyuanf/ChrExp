### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()
setwd("./")
source(file.path("src", "utils.R"))
pacman::p_load(reshape2, GenomicRnages)

### DATA 
# data <- read.table(gzfile("./data/epimap-gene-exp/merged_log2fpkm.mtx.gz"), sep = "\t", header = TRUE)
data <- read.table(gzfile("./data/epimap-gene-exp/merged_qn_log2fpkm.mtx.gz"), sep = "\t", header = TRUE)

### CONVERT TO MATRIX 
data_reshape <- reshape2::dcast(gene ~ id, value.var = "log2fpkm", data = data)
row.names(data_reshape) <- data_reshape$gene
data_reshape$gene <- NULL

### CONVERT ENSG* TO GENE_ID
# metadata <- readRDS("./data/263-HG19-CELLTYPES-RNA-METADATA.rds")
metadata <- readRDS("./data/epimap-metadata/gencode_gene_metadata_tss.rds")

### Convert metadata gene ID "ENSG*.*" to "ENSG*"
metadata$gene_id_ss = lapply(metadata$gene_id, function(x){
    strsplit(x, split='\\.')[[1]][1]
})

metadata$gene_id = metadata$gene_id_ss

common <- intersect(metadata$gene_id, row.names(data_reshape))
data_reshape <- data_reshape[common, ]
metadata <- metadata[metadata$gene_id %in% common,]

stopifnot(identical(row.names(data_reshape), unlist(metadata$gene_id)))
row.names(data_reshape) <- unlist(metadata$gene_id)
stopifnot(identical(row.names(data_reshape), unlist(metadata$gene_id)))


### EXPORT
# saveRDS(data_reshape, compress = TRUE, "./data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds")
saveRDS(data_reshape, compress = TRUE, "./data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM-QN.rds")
