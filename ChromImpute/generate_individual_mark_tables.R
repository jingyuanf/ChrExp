### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()
setwd("/u/home/f/fujy2038/project-ernst/Project/")
source(file.path("src", "utils2.R"))

### LIBS
pacman::p_load(tidyverse, GenomicRanges)

### DATA


## This is using gene body
mark_table = "/u/home/f/fujy2038/project-ernst/Project/bss-predictions/chromimpute/epimap/mark_table_full_gb.txt"

## The other one is using TSS

#diff_cov <- 100
if (!interactive()) {
  ### CLI w/ conda env 
  # <USER>:~/ source activate beast
  # <USER>:~/ Rscript --no-save --no-restore --verbose model-across-celltype-signal-mu-h3k27ac.R all 2000 > output.Rout 2>&1
  args <- commandArgs(trailingOnly = TRUE)
  mark_table <- as.numeric(args[1])
}
mark_table_full = read.table(mark_table)
names(mark_table_full) <- c("celltype", "mark", "file")

marks = unique(mark_table_full$mark)
dir.create("/u/home/f/fujy2038/project-ernst/Project/bss-predictions/chromimpute/epimap/mark_tables/")

for (mark in marks[1:length(marks)-1]){
  mark_tgt_comb = c(mark, "RNA")
  mark_table_ss = mark_table_full[mark_table_full$mark %in% mark_tgt_comb,]
  write.table(mark_table_ss, sprintf("/u/home/f/fujy2038/project-ernst/Project/bss-predictions/chromimpute/epimap/mark_tables/mark_table_%s.txt", mark), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

mark = "RNA"
mark_table_ss = mark_table_full[mark_table_full$mark == mark,]
write.table(mark_table_ss, sprintf("/u/home/f/fujy2038/project-ernst/Project/bss-predictions/chromimpute/epimap/mark_tables/mark_table_%s.txt", mark), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



## Old code for producing a subset of mark-sample combination based on training, testing and validating cell types
# mark_table_ss = mark_table_full[mark_table_full$V1 %in% ct_list,]
# # write.table(mark_table_ss, sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-h3k27ac-hg19/mark_tables/mark_table_%s_ct_%s_%s.txt", num_train_ct, num_group, ct), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# write.table(mark_table_ss, sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-chromimpute-tss/mark_tables/mark_table_%s_ct_%s_%s.txt", num_train_ct, num_group, ct), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# # mark_table_ss = read.table("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-h3k27ac-hg19/mark_tables/mark_table_20_ct_0_BSS00046.txt")
