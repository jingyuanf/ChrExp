#!/usr/bin/env Rscript
### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })
getwd()
gc()
memory.limit(size=Inf)

#.libPaths(c("/Library/Frameworks/R.framework/Versions/3.6/Resources/library", .libPaths()))
#print(.libPaths())
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))


### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, ggplot2, GenomicRanges, caret, foreach, doParallel, data.table) # parallel
source("/u/home/f/fujy2038/project-ernst/Project/batch/workflow/utils/deconvolve_utils.R")
source("/u/home/f/fujy2038/project-ernst/Project/batch/workflow/utils/evaluation_utils.R")
source("/u/home/f/fujy2038/project-ernst/Project/batch/workflow/utils/plot_utils.R")

parser <- OptionParser()
option_list <- list(
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c("-c", "--chip"), default="./data/atlas_chip.rds",
              help = "Set file path for ChIP-seq data [default \"%default\"]",
              metavar="chip-seq data"),
  make_option(c("--mark"), default="H3K27ac", help = "The mark that is used for analysis", metavar="histone mark"),
  make_option(c("-p", "--proj_name"), default="new_proj",
              help = "Define a project name. Results from this experiment will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--testing_ct"), default="./path_to_testing_ct/testing_ct.txt",
              help = "A file listing celltypes to use for testing [default \"%default\"]",
              metavar="testing celltypes"),
  make_option(c("--leave_one_out"), action='store_true', default=FALSE,
              help = "Use --leave-one-out when you are working on leave one out experiments (train with all other cell types and test with one cell type). Default is False."),
  make_option(c("--chunk_size", default=100000, help = "Chunk size for data reading into the model", metavar="chunk size")),
  make_option(c("--num_chunks", default=13, help = "Number of chunks in training", metavar="num chunk")),
  make_option(c("--model_name", default="ensembl-model", help = "Subdirectory path in prediction folder. For example, cluster-model", metavar="model name"))
  ####### TODO #########
  ### make option knn ###
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/", opt$proj_name), recursive=TRUE)
saveRDS(opt, sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/opt_handles_combine.rds", opt$proj_name))

opt$mark <- "H3K4me1"
opt$proj_name_file <- sprintf("epimap-%s-hg19-non-dup-cluster_test/epimap-%s-hg19-non-dup-cluster_test_BSS00007", opt$mark, opt$mark)
opt <- readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/opt_handles_combine.rds", opt$proj_name_file))

opt$proj_name <- sprintf("epimap-%s-hg19-non-dup-cluster_test", opt$mark)
opt$model_name <- "cluster-model"
## EXAMPLE: READ THIS OPTION HANDLE
# opt <- readRDS("/u/home/f/fujy2038/project-ernst/Project/temp_file/epimap_h3k27ac_0_20_1118_w_dist/opt_handles_combine.rds")
# opt$leave_one_out = TRUE
# opt$mark = "H3K4me1"
# opt$chip = sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-%s-hg19/FULL-ATLAS-EPIMAP-HG19-%s-RAW.rds", opt$mark, opt$mark)
# opt$proj_name = sprintf("epimap-%s-hg19-non-dup_test", opt$mark)
# opt$testing_ct = sprintf("/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-%s-hg19/all-celltypes-non-dup-%s.rds", opt$mark, opt$mark)
# opt$chunk_size = "1e+05"
# opt$num_chunks = 13

## TODO: Implement the code for combining the results
wd <- opt$wd
setwd(wd)

loo <- opt$leave_one_out
across_dir <- opt$outdir
proj_name <- opt$proj_name
chunk_size <- opt$chunk_size
num_chunks <- opt$num_chunks
model_name <- opt$model_name

testing_ct <- try(readRDS(opt$testing_ct))
# testing_ct <- "BSS00007"
if (class(testing_ct) == "try-error"){
  message("Testing cell-types file not found: Please input a valid testing cell-types file location! (Needs to be a GRANGES object)")
  testing_ct
}
message("Data: Testing cell tyeps")

ATLAS_CHIP <- try(readRDS(opt$chip))
if (class(ATLAS_CHIP) == "try-error"){
  message("ChIP-seq file not found: Please input a valid ChIP-seq file location! (Needs to be a GRANGES object)")
  ATLAS_CHIP
}
message("Data: ATLAS_CHIP")
names(mcols(ATLAS_CHIP)) <- unlist(lapply(names(mcols(ATLAS_CHIP)), function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x
}))


### COMBINE RESULTS FROM DIFFERENT CHUNKS

if (!loo) {
  pred_path = file.path(across_dir, proj_name, "predictions", model_name)
  all_chr = list.files(pred_path)
  for (chr in all_chr){
    testing_pred_df_reg_full_l <- lapply(seq(1,num_chunks), function(chunk){
      testing_pred_df_reg_full <- readRDS(file.path(pred_path, chr, sprintf("testing_pred_df_reg_full_%s_%s.rds", chunk, chunk_size)))
      testing_pred_df_reg_full
    })

    testing_pred_df_reg_full <- do.call(rbind, testing_pred_df_reg_full_l)

    testing_pred_df_reg_global_l <- lapply(seq(1,num_chunks), function(chunk){
      testing_pred_df_reg_global <- readRDS(file.path(pred_path, chr, sprintf("testing_pred_df_reg_global_%s_%s.rds", chunk, chunk_size)))
      testing_pred_df_reg_global
    })

    testing_pred_df_reg_global <- do.call(rbind, testing_pred_df_reg_global_l)

    saveRDS(testing_pred_df_reg_full, file.path(pred_path, chr, "testing_pred_df_reg_full.rds"))
    saveRDS(testing_pred_df_reg_global, file.path(pred_path, chr, "testing_pred_df_reg_global.rds"))
  }
} else {
  for (ct in testing_ct){
    proj_name_ct = paste0(proj_name, "/", proj_name, "_", ct)
    pred_path = file.path(across_dir, proj_name_ct, "predictions", model_name)
    all_chr = list.files(pred_path)

    for (chr in all_chr){
      testing_pred_df_reg_full_l <- lapply(seq(1,num_chunks), function(chunk){
        testing_pred_df_reg_full <- readRDS(file.path(across_dir, proj_name_ct, "predictions", model_name, chr, sprintf("testing_pred_df_reg_full_%s_%s.rds", chunk, chunk_size)))
        testing_pred_df_reg_full
      })

      testing_pred_df_reg_full <- do.call(rbind, testing_pred_df_reg_full_l)

      testing_pred_df_reg_global_l <- lapply(seq(1,num_chunks), function(chunk){
        testing_pred_df_reg_global <- readRDS(file.path(across_dir, proj_name_ct, "predictions", model_name, chr, sprintf("testing_pred_df_reg_global_%s_%s.rds", chunk, chunk_size)))
        testing_pred_df_reg_global
      })

      testing_pred_df_reg_global <- do.call(rbind, testing_pred_df_reg_global_l)

      saveRDS(testing_pred_df_reg_full, file.path(across_dir, proj_name_ct, "predictions", model_name, chr, "testing_pred_df_reg_full.rds"))
      saveRDS(testing_pred_df_reg_global, file.path(across_dir, proj_name_ct, "predictions", model_name, chr, "testing_pred_df_reg_global.rds"))
    }
  }
}

### COMBINE RESULTS FROM DIFFERENT TESTING CELL TYPES
# testing_ct="BSS00007"

if (!loo) {
  # chr=all_chr[1]
  pred_path = file.path(across_dir, proj_name, "predictions", model_name)
  all_chr = list.files(pred_path)
  pred_full <- readRDS(file.path(pred_path, chr, "testing_pred_df_reg_full.rds"))
  pred_global <- readRDS(file.path(pred_path, chr, "testing_pred_df_reg_global.rds"))

  chip_gr_full <- granges(ATLAS_CHIP)
  mcols(chip_gr_full) <- data.frame(matrix(nrow=nrow(mcols(ATLAS_CHIP)), ncol=ncol(pred_full)))
  names(mcols(chip_gr_full)) <- colnames(pred_full)

  chip_gr_global <- granges(ATLAS_CHIP)
  mcols(chip_gr_global) <- data.frame(matrix(nrow=nrow(mcols(ATLAS_CHIP)), ncol=ncol(pred_global)))
  names(mcols(chip_gr_global)) <- colnames(pred_global)

  for (chr in all_chr){
    pred_full <- try(readRDS(file.path(pred_path, chr, "testing_pred_df_reg_full.rds")))
    pred_global <- try(readRDS(file.path(pred_path, chr, "testing_pred_df_reg_global.rds")))
    if (class(pred_full) == "try-error"){
      message(sprintf("Prediction file of %s not found!", chr))
      pred_full
    } else {
      mcols(chip_gr_full)[seqnames(chip_gr_full) == chr] <- pred_full
    }

    if (class(pred_global) == "try-error"){
      message(sprintf("Prediction file of %s not found!", chr))
      pred_global
    } else {
      mcols(chip_gr_global)[seqnames(chip_gr_global) == chr] <- pred_global
    }
  }

  dir.create(file.path(across_dir, proj_name, "final-pred"), recursive=TRUE)
  saveRDS(chip_gr_full, file.path(across_dir, proj_name, "final-pred", "testing_pred_reg_full.rds"))
  saveRDS(chip_gr_global, file.path(across_dir, proj_name, "final-pred", "testing_pred_reg_global.rds"))

} else {
  chip_gr_full <- granges(ATLAS_CHIP)
  mcols(chip_gr_full) <- data.frame(matrix(nrow=nrow(mcols(ATLAS_CHIP)), ncol=length(testing_ct)))
  names(mcols(chip_gr_full)) <- testing_ct

  chip_gr_global <- granges(ATLAS_CHIP)
  mcols(chip_gr_global) <- data.frame(matrix(nrow=nrow(mcols(ATLAS_CHIP)), ncol=length(testing_ct)))
  names(mcols(chip_gr_global)) <- testing_ct

  for (ct in testing_ct){
    proj_name_ct = paste0(proj_name, "/", proj_name, "_", ct)
    pred_path = file.path(across_dir, proj_name_ct, "predictions", model_name)
    all_chr = list.files(pred_path)

    for (chr in all_chr){
      pred_full <- try(readRDS(file.path(pred_path, chr, "testing_pred_df_reg_full.rds")))
      pred_global <- try(readRDS(file.path(pred_path, chr, "testing_pred_df_reg_global.rds")))
      if (class(pred_full) == "try-error"){
        message(sprintf("Prediction file of %s not found!", chr))
        pred_full
      } else {
        mcols(chip_gr_full)[seqnames(chip_gr_full) == chr,ct] <- pred_full
      }

      if (class(pred_global) == "try-error"){
        message(sprintf("Prediction file of %s not found!", chr))
        pred_global
      } else {
        mcols(chip_gr_global)[seqnames(chip_gr_global) == chr,ct] <- pred_global
      }
    
      # output_file = file.path(across_dir, proj_name_combined, "final-pred", sprintf("%s_%s_prediction.bw", chr, ct))
      # chip_gr_ct = chip_gr[seqnames(chip_gr) == chr,ct]
      # names(mcols(chip_gr_ct)) = "score"
      # export.granges(chip_gr, file = output_file, format = "bigwig")
    }

  }
  proj_name_combined = proj_name
  dir.create(file.path(across_dir, proj_name_combined, "final-pred"), recursive=TRUE)
  saveRDS(chip_gr_full, file.path(across_dir, proj_name, "final-pred", "testing_pred_reg_full.rds"))
  saveRDS(chip_gr_global, file.path(across_dir, proj_name, "final-pred", "testing_pred_reg_global.rds"))
}


