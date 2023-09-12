#!/usr/bin/env Rscript

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
pacman::p_load(tidyverse, ggplot2, GenomicRanges, caret, foreach, doParallel, data.table, rtracklayer) # parallel

parser <- OptionParser()
option_list <- list(
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c("--proj_name_test"), default="new_proj_test",
              help = "Define a project name. Prediction results for testing cell types will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name test"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--save_intermediate"), action='store_true', default=FALSE,
              help = "Use this flag if you want to save all intermediate files. Otherwise only the models will be saved and other intermediate files will be deleted.")

  ####### TODO #########
  ### make option knn ###
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))

wd <- opt$wd
setwd(wd)
source(file.path("batch", "workflow", "utils", "utils.R"))
source(file.path("batch", "workflow", "ChromExp", "model_helper_functions.R"))

save_intermediate <- opt$save_intermediate

proj_name_test <- opt$proj_name_test
across_dir <- opt$outdir


if (!save_intermediate){
  unlink(file.path(across_dir, proj_name_test, "intermediate-files"), recursive = TRUE, force = TRUE)
}