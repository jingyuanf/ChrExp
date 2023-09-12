### THIS CODE IS FOR FORMATTING ATLAS_GENE AND ATLAS_GENE_METADATA
### ROW NAMES WILL BE GENE NAMES
### ROWS WILL BE SORTED BY GENE NAMES


### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()


### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, ggplot2, GenomicRanges) # parallel


if (interactive()){
  wd <- "/u/home/f/fujy2038/project-ernst/Project/"
  setwd(wd)
  source(file.path("batch", "workflow", "utils", "utils.R"))

  ATLAS_GENE_train <- readRDS("./Correct_data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds")
  ATLAS_GENE_test <- readRDS("./Correct_data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds")
  ATLAS_GENE_METADATA <- readRDS("./Correct_data/epimap-metadata/gencode_gene_metadata.rds")

  new_file_location <- "epimap-gene-exp"
}


if (!interactive()) {
  ### CLI w/ conda env 
  # <USER>:~/ source activate beast
  # <USER>:~/ Rscript --no-save --no-restore --verbose model-across-celltype-signal-mu-h3k27ac.R all 2000 > output.Rout 2>&1
  # see example bash file
  args <- commandArgs(trailingOnly = TRUE)
  wd <- args[1]
  setwd(wd)
  source(file.path("src", "utils2.R"))


  ATLAS_GENE_train <- try(readRDS(args[2]))
  if (class(ATLAS_GENE_train) == "try-error"){
    message("Gene expression training file not found: Please input a valid gene expression training file location! ")
    ATLAS_GENE_train
  }
  message("Data: ATLAS_GENE_train")

  ATLAS_GENE_test <- try(readRDS(args[3]))
  if (class(ATLAS_GENE_test) == "try-error"){
    message("Gene expression testing file not found: Please input a valid gene expression testing file location! ")
    ATLAS_GENE_test
  }
  message("Data: ATLAS_GENE_test")


  ATLAS_GENE_METADATA <- try(readRDS(args[4]))
  if (class(ATLAS_GENE_METADATA) == "try-error"){
    message("Gene expression meta file not found: Please input a valid gene expression meta file location! ")
    ATLAS_GENE_METADATA
  }
  message("Data: ATLAS_GENE_METADATA")
  
  new_file_location <- args[5]
}


ATLAS_GENE_METADATA$gene_id_ss = lapply(ATLAS_GENE_METADATA$gene_id, function(x){
    strsplit(x, split='\\.')[[1]][1]
})
ATLAS_GENE_METADATA$gene_id <- unlist(ATLAS_GENE_METADATA$gene_id_ss)


if(length(grep("ENSG", row.names(ATLAS_GENE_train))) > 0){
  common_gene_id <- intersect(row.names(ATLAS_GENE_train), ATLAS_GENE_METADATA$gene_id)
  common_rows <- which(ATLAS_GENE_METADATA$gene_id %in% common_gene_id)

  ATLAS_GENE_METADATA <- ATLAS_GENE_METADATA[common_rows,]
  ATLAS_GENE_train <- ATLAS_GENE_train[ATLAS_GENE_METADATA$gene_id,]

  check_dup <- which(duplicated(ATLAS_GENE_METADATA$gene_name))
  if (length(check_dup) > 1){
    ATLAS_GENE_METADATA <- ATLAS_GENE_METADATA[-check_dup,]
    ATLAS_GENE_train <- ATLAS_GENE_train[-check_dup,]
  }
  row.names(ATLAS_GENE_train) <- ATLAS_GENE_METADATA$gene_name

} else {
  common_gene_name <- intersect(row.names(ATLAS_GENE_train), ATLAS_GENE_METADATA$gene_name)
  common_rows <- which(ATLAS_GENE_METADATA$gene_name %in% common_gene_name)

  ATLAS_GENE_METADATA <- ATLAS_GENE_METADATA[common_rows,]
  ATLAS_GENE_train <- ATLAS_GENE_train[ATLAS_GENE_METADATA$gene_name,]
}



## Then, align ATLAS_GENE_test with ATLAS_GENE_METADATA

if(length(grep("ENSG", row.names(ATLAS_GENE_test))) > 0){
  common_gene_id <- intersect(row.names(ATLAS_GENE_test), ATLAS_GENE_METADATA$gene_id)
  common_rows <- which(ATLAS_GENE_METADATA$gene_id %in% common_gene_id)

  ATLAS_GENE_METADATA <- ATLAS_GENE_METADATA[common_rows,]
  ATLAS_GENE_train <- ATLAS_GENE_train[common_rows,]
  ATLAS_GENE_test <- ATLAS_GENE_test[ATLAS_GENE_METADATA$gene_id,]
  row.names(ATLAS_GENE_test) <- ATLAS_GENE_METADATA$gene_name

} else {
  common_gene_name <- intersect(row.names(ATLAS_GENE_test), ATLAS_GENE_METADATA$gene_name)
  common_rows <- which(ATLAS_GENE_METADATA$gene_name %in% common_gene_name)

  ATLAS_GENE_METADATA <- ATLAS_GENE_METADATA[common_rows,]
  ATLAS_GENE_train <- ATLAS_GENE_train[common_rows,]
  ATLAS_GENE_test <- ATLAS_GENE_test[ATLAS_GENE_METADATA$gene_name,]
}

dir.create(file.path(wd, new_file_location))
saveRDS(ATLAS_GENE_METADATA, file.path(wd, "Correct_data", new_file_location, "ATLAS_GENE_METADATA.rds"))
saveRDS(ATLAS_GENE_train, file.path(wd, "Correct_data", new_file_location, "ATLAS_GENE_train.rds"))
saveRDS(ATLAS_GENE_test, file.path(wd, "Correct_data", new_file_location, "ATLAS_GENE_test.rds"))

### OPTIONAL: CHECK GENE EXPRESSION DISTRIBUTION
mean_l <- list()
sd_l <- list()
for (ct in colnames(ATLAS_GENE_train)){
  g_mean <- mean(ATLAS_GENE_train[,ct])
  g_sd <- sd(ATLAS_GENE_train[,ct])

  mean_l <- append(mean_l, g_mean)
  sd_l <- append(sd_l, g_sd)
}

mean_l <- unlist(mean_l)
sd_l <- unlist(sd_l)

summary(mean_l)
summary(sd_l)

q(save = "no")
