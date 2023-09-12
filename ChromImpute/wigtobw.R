### Environment ###
rm(list=ls(all=TRUE))
options(warn=-1)
options(scipen=999)
setwd("./")
source("./utils/utils.R")

library(rtracklayer)
library(tools)

mark="H3K36me3"
proj_name=sprintf("epimap-%s-chromimpute-rna-all-pred-chr1", mark)
chr="chr1"
chr_pred="chr1"
testing_ct="./data/epimap-H3K36me3-hg19/all-celltypes-non-dup-H3K36me3.rds"
chrimp_dir <- "./bss-predictions/chromimpute"

### NEW OPTION HANDLE FOR LEAVE-ONE-OUT EXPERIMENTS ON MULTIPLE MARKERS
if (!interactive()) {
  ### CLI w/ conda env 
  # <USER>:~/ source activate beast
  # <USER>:~/ Rscript --no-save --no-restore --verbose model-across-celltype-signal-mu-h3k27ac.R all 2000 > output.Rout 2>&1
  args <- commandArgs(trailingOnly = TRUE)
  # sc_k27_file_path <- args[1]
  # file <- args[2]
  mark = args[1]
  proj_name = args[2]
  chr = args[3]
  testing_ct = args[4]
  chrimp_dir = args[5]
}



ATLAS_CHIP <- readRDS("./data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")
seqinfo_track <- seqinfo(ATLAS_CHIP)
testing_celltypes <- readRDS(testing_ct)

list_files <- list.files(sprintf("%s/%s/output/", chrimp_dir, proj_name), pattern = ".wig.gz")

## Only get the predictions for the chromosome that we need
if (chr != "all"){
  files_ss <- list_files[grepl(paste0(chr, "_impute"), list_files)]
} else {
  files_ss <- list_files
}

dir.create(sprintf("%s/%s/output_bw/tracks/log2/", chrimp_dir, proj_name), recursive = TRUE)
dir.create(sprintf("%s/%s/output_bw/tracks/raw/", chrimp_dir, proj_name), recursive = TRUE)


for(ct in testing_celltypes){
  print(ct)

  files_ct <- files_ss[grepl(ct, files_ss)]
  all_chr <- lapply(files_ct, function(x){
    strsplit(x, split="_")[[1]][1]
  })
  all_chr = unique(unlist(all_chr))

  all_wig = lapply(all_chr, function(chr_pred){
      wig_file <- import.wig(sprintf("%s/%s/output/%s_impute_%s_%s.wig.gz", chrimp_dir, proj_name, chr_pred, ct, mark))
      keep_rows = dim(mcols(ATLAS_CHIP[seqnames(ATLAS_CHIP) == chr_pred, ]))[1]
      return(GRanges(wig_file[1:keep_rows,]))
  })

  all_gr = unlist(as(all_wig, "GRangesList"))
  all_gr$score[is.na(all_gr$score)] = 0
  
  all_gr_log2 = all_gr
  mcols(all_gr_log2) = log2(as.matrix(mcols(all_gr)) + 1)

  seqinfo(all_gr) <- seqinfo_track
  seqinfo(all_gr_log2) <- seqinfo_track

  rtracklayer::export.bw(all_gr, sprintf("%s/%s/output_bw/tracks/raw/%s_ChromImpute.bw", chrimp_dir, proj_name, ct))
  rtracklayer::export.bw(all_gr_log2, sprintf("%s/%s/output_bw/tracks/log2/%s_ChromImpute.bw", chrimp_dir, proj_name, ct))
}

