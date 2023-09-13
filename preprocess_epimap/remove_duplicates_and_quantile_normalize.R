#!/usr/bin/env Rscript

### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()
setwd("./mypath")
source(file.path("batch", "workflow", "utils", "utils.R"))
suppressPackageStartupMessages(library("optparse"))

pacman::p_load(GenomicRanges, rtracklayer, tidyverse, preprocessCore, readxl)

### GET PARAMETERS
parser <- OptionParser()
option_list <- list(
  make_option(c("-m", "--mark"), default="H3K27ac",
              help="Specify the chromatin mark that you want to download from EpiMap. Default is H3K27ac.", 
              metavar="histone mark"),
  make_option(c("--wig"), action='store_true', default=FALSE,
              help="Use --wig if you want your results to be output in wig format. If you want '.bw' format then don't use this flag.", 
              metavar="Use wig or bw"),
  make_option(c("--rm_dup"), action='store_true', default=FALSE,
              help="Use --rm_dup if you want to remove duplicates. If you want all cell types then don't use this flag.")
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))
mark <- opt$mark
use_wig <- opt$wig
rm_dup <- opt$rm_dup

# mark = "H3K36me3"

### READ IN GENE EXP AND CHROMATIN SIGNAL INFORMATION
ATLAS_GENE <- readRDS("./data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds")
ATLAS_CHIP <- readRDS(file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("FULL-ATLAS-EPIMAP-HG19-%s-RAW.rds", mark)))

### GET CELL TYPE NAMES FROM GENE EXPRESSION AND CHROMATIN SIGNAL 
celltypes_exp <- colnames(ATLAS_GENE)
celltypes_chip <- colnames(mcols(ATLAS_CHIP))
celltypes_overlap <- intersect(celltypes_exp, celltypes_chip)

saveRDS(celltypes_overlap, file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("all-celltypes-with-dup-%s.rds", mark)))
writeLines(celltypes_overlap, file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("all-celltypes-with-dup-%s.txt", mark)))

### GET ANNOTATION FOR EPIMAP CELLTYPES
if (rm_dup){
  epimap_meta <- as.data.frame(read_excel("./data/Imputation_Metadata.xlsx"))

  ### FROM OVERLAPPING CELL TYPES, REMOVE DUPLICATED ONES
  ### PS: Duplicated cell types are defined by "Extended Info", "Lifestage", "Age", "AgeUnits", "Sex", "Project"
  epimap_meta_overlap <- epimap_meta[epimap_meta$BSSID %in% celltypes_overlap, ]
  epimap_meta_overlap$unique_info <- paste0(epimap_meta_overlap$`Extended Info`, "_", epimap_meta_overlap$Lifestage, "_", epimap_meta_overlap$Age, "_", epimap_meta_overlap$AgeUnits, "_", epimap_meta_overlap$Sex, "_", epimap_meta_overlap$Project)
  unique_status <- unique(epimap_meta_overlap$unique_info)


  celltype_ind <- split(seq_along(epimap_meta_overlap$unique_info), epimap_meta_overlap$unique_info)
  unique_ind <- lapply(unique_status, function(x){
    return(celltype_ind[[x]][1])
  })
  unique_ind <- unlist(unique_ind)

  celltypes_rm_dup <- celltypes_overlap[unique_ind]

  saveRDS(celltypes_rm_dup, file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("all-celltypes-non-dup-%s.rds", mark)))
  writeLines(celltypes_rm_dup, file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("all-celltypes-non-dup-%s.txt", mark)))
}

### SUBSET CHROMATIN SIGNAL DATA TO NON-DUPLICATING CELL TYPES
if (rm_dup){
  ATLAS_CHIP <- ATLAS_CHIP[,celltypes_rm_dup] # Cell types after removing duplicates
} else {
  ATLAS_CHIP <- ATLAS_CHIP[,celltypes_overlap] # Overlapping cell types, without removing duplicates
}

### NORMALIZE QUANTILES
ATLAS_CHIP_qn <- granges(ATLAS_CHIP)
ATLAS_CHIP_df <- as.matrix(mcols(ATLAS_CHIP))

ATLAS_CHIP_df_qn <- normalize.quantiles(ATLAS_CHIP_df)
colnames(ATLAS_CHIP_df_qn) <- colnames(ATLAS_CHIP_df)

ATLAS_CHIP_df_qn[is.na(ATLAS_CHIP_df_qn)] <- 0 # Remove all NA's make them 0
ATLAS_CHIP_df[is.na(ATLAS_CHIP_df)] <- 0

mcols(ATLAS_CHIP_qn) <- ATLAS_CHIP_df_qn
mcols(ATLAS_CHIP) <- ATLAS_CHIP_df

### SAVE QUANTILE NORMALIZED ATLAS CHIP
saveRDS(ATLAS_CHIP_qn, file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("FULL-EPIMAP-HG19-%s-QN.rds", mark)))

### OUTPUT QUANTILE NORMALIZED TRACKS 
ATLAS_CHIP_qn <- readRDS(file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("FULL-EPIMAP-HG19-%s-QN.rds", mark)))
chr_list = paste0("chr", c(seq(1,22), "X"))
dir.create(file.path("data", sprintf("epimap-%s-hg19", mark), "qn-tracks"), recursive=TRUE)
dir.create(file.path("data", sprintf("epimap-%s-hg19", mark), "raw-tracks"), recursive=TRUE)
dir.create(file.path("data", sprintf("epimap-%s-hg19", mark), "qn-tracks-wig"), recursive=TRUE)
dir.create(file.path("data", sprintf("epimap-%s-hg19", mark), "raw-tracks-wig"), recursive=TRUE)


for (ct in celltypes_rm_dup){
  print(ct)
  track_qn <- ATLAS_CHIP_qn[,ct]
  track_raw <- ATLAS_CHIP[,ct]

  names(mcols(track_qn)) <- "score"
  names(mcols(track_raw)) <- "score"

  for (chr in chr_list){
    track_ss_qn <- track_qn[seqnames(track_qn) == chr,]
    track_ss_raw <- track_raw[seqnames(track_raw) == chr,]

    if (use_wig){
      rtracklayer::export.wig(track_ss_qn, file.path("data", sprintf("epimap-%s-hg19", mark), "qn-tracks-wig", sprintf("%s-%s-HG19-%s-QN.wig", ct, chr, mark)))
      rtracklayer::export.wig(track_ss_raw, file.path("data", sprintf("epimap-%s-hg19", mark), "raw-tracks-wig", sprintf("%s-%s-HG19-%s-RAW.wig", ct, chr, mark)))
    } else {
      rtracklayer::export.bw(track_ss_qn, file.path("data", sprintf("epimap-%s-hg19", mark), "qn-tracks", sprintf("%s-%s-HG19-%s-QN.bw", ct, chr, mark)))
      rtracklayer::export.bw(track_ss_raw, file.path("data", sprintf("epimap-%s-hg19", mark), "raw-tracks", sprintf("%s-%s-HG19-%s-RAW.bw", ct, chr, mark)))
    }
  }

  if (use_wig){
    rtracklayer::export.wig(track_qn, file.path("data", sprintf("epimap-%s-hg19", mark), "qn-tracks-wig", sprintf("%s-FULL-HG19-%s-QN.wig", ct, mark)))
    rtracklayer::export.wig(track_raw, file.path("data", sprintf("epimap-%s-hg19", mark), "raw-tracks-wig", sprintf("%s-FULL-HG19-%s-RAW.wig", ct, mark)))
  } else {
    rtracklayer::export.bw(track_qn, file.path("data", sprintf("epimap-%s-hg19", mark), "qn-tracks", sprintf("%s-FULL-HG19-%s-QN.bw", ct, mark)))
    rtracklayer::export.bw(track_raw, file.path("data", sprintf("epimap-%s-hg19", mark), "raw-tracks", sprintf("%s-FULL-HG19-%s-RAW.bw", ct, mark)))
  }

}
