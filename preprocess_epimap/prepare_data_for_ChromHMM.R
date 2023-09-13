#!/usr/bin/env Rscript
rm(list=ls(all=TRUE))
options(warn=-1)
options(scipen=999)
source("./utils/utils.R")
source("./utils/deconvolve_utils.R")

library(rtracklayer)
library(tools)
library(GenomicRanges)

# chrom_info <- read.table("./data/epimap-hg19/chrom_info.txt")
# names(chrom_info) <- c("chr", "length")

chr = "chr1"
mark = "H3K4me1"
ct = "BSS00007"
wd = "./"

setwd(wd)

## GET BED FILES FROM MACS2 FOR PEAK CALLING
proj_name = sprintf("epimap-%s-hg19", mark)
bfile = sprintf("./data/%s/raw-bed/%s-%s-HG19-%s-RAW.bed", proj_name, ct, chr, mark)
bdata = read.table(bfile, skip = 1, header = FALSE, sep ='\t')

bdata_region <- bdata[,1:3]
names(bdata_region) <- c("chr", "start", "end")

bdata_gr <- makeGRangesFromDataFrame(bdata_region)

## GET 200BP CHROMATIN GRANGES AND FIND OVERLAPS
region_200bp <- readRDS("./data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")
region_chr = subset_to_chr(region_200bp,chr)

region_len = length(region_chr)
overlaps <- findOverlaps(region_chr, bdata_gr)

## ANNOTATE PEAK REGIONS IN THE GRANGES WITH 1'S, THE OTHER REGION AS ZEROS
peak_df = data.frame(matrix(0, nrow=region_len, ncol=1))
colnames(peak_df) <- mark
peak_df[queryHits(overlaps),1] = 1 

## WRITE INTO BINARY FILE
dir.create(sprintf("./data/%s/raw-binarized/", proj_name))
ofile = sprintf("./data/%s/raw-binarized/%s-%s-HG19-%s_binary.txt", proj_name, ct, chr, mark)
writeLines(paste(c(ct, chr), collapse='\t'), con=ofile)
write.table(peak_df, file=ofile, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)
