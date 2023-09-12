### Environment ###
rm(list=ls(all=TRUE))
options(warn=-1)
options(scipen=999)
setwd("/u/home/f/fujy2038/project-ernst/Project/")
source(file.path("batch", "workflow", "utils", "utils.R"))

pacman::p_load(tidyverse, GenomicRanges, plyr, dplyr) # caret, doParallel

chip <- readRDS("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")
seqlength_df = as.data.frame(seqlengths(chip))
write.table(seqlength_df, "/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-hg19/chrom_info.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
write.table(seqlength_df, "/u/home/f/fujy2038/project-ernst/Project/bss-predictions/chromimpute/epimap/chrom_info.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
