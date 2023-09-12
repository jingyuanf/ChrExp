### Environment ###
rm(list=ls(all=TRUE))
options(warn=-1)
options(scipen=999)
setwd("./")
source("./utils/utils.R")

pacman::p_load(tidyverse, GenomicRanges, plyr, dplyr) # caret, doParallel
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

source(file.path("batch", "workflow", "utils", "utils.R"))

parser <- OptionParser()
option_list <- list(
    make_option(c("-m", "--mark"), default="./markers/list_of_markers.txt",
        help="Total list of markers to work with")
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))

mark_f <- opt$mark
mark_f <- "./markers/list_of_markers.txt"
marks <- readLines(mark_f)

mark_table <- data.frame(matrix(nrow=0, ncol=3))
colnames(mark_table) = c("sample", "mark", "file")


## Append histone modification files to mark table
for (mark in marks){
  filedir = sprintf("./data/epimap-%s-hg19/raw-tracks-wig/", mark)
  # filedir = "./data/mouse-chromimpute/tracks"
  # dir.create(filedir, recursive=TRUE)

  file_names = list.files(filedir, pattern = "*-FULL-*")
  celltypes = lapply(file_names, function(x){
    strsplit(x, split="-")[[1]][1]
  })
  celltypes = unlist(celltypes)

  # ## Add this section for adding mark-dependent directories ## TODO: Figure this out
  # file_names_dir = lapply(file_names, function(x){
  #   paste0(sprintf("/epimap-%s-hg19/raw-tracks-wig/", mark), x)
  # })
  # file_names_dir = unlist(file_names_dir)

  new_table = data.frame(sample = celltypes, mark = mark, file = file_names)
  mark_table = rbind(mark_table, new_table)
}

## Append RNA files to mark table (Make two separate mark tables for TSS annotation and Gene Body annotation)
filedir_tss = "./data/epimap-gene-exp/raw-tracks-tss/"
file_names = list.files(filedir_tss)
celltypes = lapply(file_names, function(x){
  strsplit(x, split="_")[[1]][1]
})
celltypes = unlist(celltypes)

# file_names_dir = lapply(file_names, function(x){
#   paste0(sprintf("epimap-gene-exp/raw-tracks-tss/"), x)
# })
# file_names_dir = unlist(file_names_dir)

new_table = data.frame(sample = celltypes, mark = "RNA", file = file_names)
mark_table_full_tss = rbind(mark_table, new_table)


filedir_gb = "./data/epimap-gene-exp/raw-tracks-genebody/"
file_names = list.files(filedir_gb)
celltypes = lapply(file_names, function(x){
  strsplit(x, split="_")[[1]][1]
})
celltypes = unlist(celltypes)

new_table = data.frame(sample = celltypes, mark = "RNA", file = file_names)
mark_table_full_gb = rbind(mark_table, new_table)

## Save mark tables
write.table(mark_table_full_tss, "./bss-predictions/chromimpute/epimap/mark_table_full_tss_w_dir.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(mark_table_full_gb, "./bss-predictions/chromimpute/epimap/mark_table_full_gb_w_dir.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

