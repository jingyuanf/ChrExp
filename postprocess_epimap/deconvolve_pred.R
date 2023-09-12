rm(list=ls(all=TRUE))
options(warn=-1)
options(scipen=999)
setwd("/u/home/f/fujy2038/project-ernst/Project/")


pacman::p_load(tidyverse, ggplot2, GenomicRanges, plyr, dplyr) # caret, doParallel


### OLD OPTION HANDLES FOR TRAIN-VALID-TEST SPLIT
# num_cell_to_convolve = 62
# num_train_ct = 20
# num_group = 0
# chr = "chr1"
# #diff_cov <- 100
# if (!interactive()) {
#   ### CLI w/ conda env 
#   # <USER>:~/ source activate beast
#   # <USER>:~/ Rscript --no-save --no-restore --verbose model-across-celltype-signal-mu-h3k27ac.R all 2000 > output.Rout 2>&1
#   args <- commandArgs(trailingOnly = TRUE)
#   num_cell_to_convolve <- as.numeric(args[1])
#   num_train_ct <- as.numeric(args[2]) # filename
#   num_group <- as.numeric(args[3])
#   chr <- args[4]
# }

### NEW OPTION HANDLE FOR DIFFERENT MARKS AND LEAVE-ONE-OUT EXPERIMENTS
mark = "H3K36me3"
chr = "chr1"
proj_name = sprintf("epimap-%s-hg19-non-dup", mark)
pred_dir = "/u/home/f/fujy2038/project-ernst/Project/bss-predictions/"
testing_ct = sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-%s-hg19/all-celltypes-non-dup-%s.rds", mark, mark)
training_ct = sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-%s-hg19/all-celltypes-non-dup-%s.rds", mark, mark)
bulk_tracks_f = sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-%s-hg19/FULL-ATLAS-EPIMAP-HG19-%s-RAW.rds", mark, mark)
#diff_cov <- 100
if (!interactive()) {
  ### CLI w/ conda env 
  # <USER>:~/ source activate beast
  # <USER>:~/ Rscript --no-save --no-restore --verbose model-across-celltype-signal-mu-h3k27ac.R all 2000 > output.Rout 2>&1
  args <- commandArgs(trailingOnly = TRUE)
  mark <- args[1]
  chr <- args[2]
  proj_name <- args[3]
  dir_name_chr <- args[4]
  pred_dir <- args[5]
  testing_ct <- args[6]
  training_ct <- args[7]
  bulk_tracks_f <- args[8]
}

source("/u/home/f/fujy2038/project-ernst/Project/batch/workflow/utils/utils.R")
source("/u/home/f/fujy2038/project-ernst/Project/batch/workflow/utils/deconvolve_utils.R")


across_dir = sprintf("%s/across-model/", pred_dir)
chr_dir = sprintf("%s/chromimpute/%s/", pred_dir, dir_name_chr)

proj_name_across = sprintf("%s_test", proj_name)
proj_name_chr = sprintf("output_bw")


### This is specifically for leave-one-out experiments
total_testing_celltypes <- readRDS(testing_ct)
total_training_celltypes <- readRDS(training_ct)

bulk_tracks <- readRDS(bulk_tracks_f)
# ATLAS_GENE_TEST <- readRDS("/u/home/f/fujy2038/project-ernst/Project/Correct_data/preprocessed_files/ATLAS_GENE_test.rds")
cells_to_convolve <- as.character(total_testing_celltypes)

convolved <- granges(bulk_tracks)
convolved$score <- rowSums(as.data.frame(mcols(bulk_tracks)[,cells_to_convolve]))
convolved <- subset_to_chr(convolved, chr)

## ACROSS CELLTYPE PREDICTION ##
if(!file.exists(file.path(across_dir, proj_name_across, "final-pred", chr, "across_gr.rds"))){
  across_gr <- process_across(across_dir, proj_name_across, chr)
} else {
  across_gr <- readRDS(file.path(across_dir, proj_name_across, "final-pred", chr, "across_gr.rds"))
}

## WITHIN CELLTYPE PREDICTION ##
# if(!file.exists(file.path("bss-predictions", "within-model", proj_name_within, "final_pred", "2", "with_distance", num_train_ct, "bw_files", sprintf("within_gr_log2_%s.rds", chr)))){
#   within_gr <- process_within(pred_dir, proj_name_within, chr)
# } else {
#   within_gr <- readRDS(file.path("bss-predictions", "within-model", proj_name_within, "final_pred", "2", "with_distance", num_train_ct, "bw_files", sprintf("within_gr_log2_%s.rds", chr)))
# }

## CHROMIMPUTE CELLTYPE PREDICTION ##
if(!file.exists(file.path(chr_dir, proj_name_chr, sprintf("chrom_gr_%s.rds", chr)))){
  chrom_gr <- process_chrom(chr_dir, proj_name_chr, chr)
} else {
  chrom_gr <- readRDS(file.path(chr_dir, proj_name_chr, sprintf("chrom_gr_%s.rds", chr)))
}

# thrd=5
# int_gr <- process_int(pred_dir, proj_name_int, chr, thrd)
# cells_to_convolve <- intersect(intersect(names(mcols(across_gr)), names(mcols(chrom_gr))), cells_to_convolve)

across_gr_decon <- deconvolve(across_gr, cells_to_convolve, convolved)
chrom_gr_decon <- deconvolve(chrom_gr, cells_to_convolve, convolved)
# within_gr_decon <- deconvolve(within_gr, cells_to_convolve, convolved)
# int_gr_decon <- deconvolve(int_gr, cells_to_convolve, convolved)

proj_name <- sprintf("%s_%s", proj_name, chr)

dir.create(file.path("deconvolved_tracks", proj_name, "leave_one_out"), recursive=TRUE)
saveRDS(across_gr_decon, sprintf("./deconvolved_tracks/%s/leave_one_out/across_decon.rds", proj_name))
print("Across deconvolved saved!")
saveRDS(chrom_gr_decon, sprintf("./deconvolved_tracks/%s/leave_one_out/chrom_decon.rds", proj_name))
print("ChromImpute deconvolved saved!")

# saveRDS(within_gr_decon, sprintf("./deconvolved_tracks/%s/%s/within_decon.rds", proj_name, num_train_ct))
# saveRDS(int_gr_decon, sprintf("./deconvolved_tracks/%s/%s/int_decon.rds", proj_name, num_train_ct))


##### OLD CELLTYPES INFO IGNORE #####
# # total_testing_celltypes = readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/holdout-celltypes-rm-dup/BSS_training_validation_celltypes/testing_ct_16_celltypes_%s.rds", num_group))
# total_testing_celltypes = readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/holdout-celltypes-rm-dup/BSS_training_validation_celltypes/testing_ct_62_celltypes_%s.rds", num_group))
# # total_validing_celltypes = readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/holdout-celltypes-rm-dup/BSS_training_validation_celltypes/validing_ct_16_celltypes_%s_%s.rds", num_train_ct, num_group))
# total_training_celltypes <- readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/holdout-celltypes-rm-dup/BSS_training_validation_celltypes/training_ct_%s_celltypes_%s.rds", num_train_ct, num_group))
# bulk_tracks <- readRDS("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-h3k27ac-hg19/QN-FULL-ATLAS-EPIMAP-HG19-H3K27AC-RAW.rds")
# # bulk_tracks <- readRDS("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-h3k27ac-hg19/epimap_110822_10kb/FULL-ATLAS-EPIMAP-HG19-H3K27AC-RAW-10kb.rds")
