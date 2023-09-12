#!/usr/bin/env Rscript

### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()
setwd("./")
source(file.path("batch", "workflow", "utils", "utils.R"))
suppressPackageStartupMessages(library("optparse"))

pacman::p_load(reshape2, GenomicRanges, rtracklayer, tidyverse, XML)

### GET PARAMETERS
parser <- OptionParser()
option_list <- list(
  make_option(c("-m", "--mark"), default="H3K27ac",
              help="Specify the chromatin mark that you want to download from EpiMap. Default is H3K27ac.", 
              metavar="histone mark")
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))
mark <- opt$mark

# mark = "H3K36me3"
### EXAMPLE BINNED DATA 
#si <- readRDS(file.path("data", "encode-atac-hg19", "hg19-seqinfo.rds"))
tmp <- readRDS(file.path("./data/GR-FULL-ATLAS-EPIMAP-HG19-H3K27AC-LOG2.rds")) # GRanges object with 15181508 ranges and 206 metadata columns, width = 200
mcols(tmp) <- data.frame(matrix(nrow = 15181508, ncol=206))
bins <- granges(tmp)

### GET TOTAL CELLTYPES FROM RNA
rna <- readRDS(file.path("./data/epimap-gene-exp/FULL-ATLAS-EPIMAP-HG19-RNA-GENEID-LOG2FPKM.rds"))
celltypes <- names(rna)


### GET EPIGENOME DATA FROM WEBSITE
url <- "https://epigenome.wustl.edu/epimap/data/observed/"
webpage <- RCurl::getURL(url)
tc <- textConnection(webpage)
contents <- readLines(tc)
files <- getHTMLLinks(contents)
files <- grep(mark, files, value = TRUE)
files_celltypes <- gsub(".sub", "", sapply(strsplit(files, "_", fixed=TRUE), function(x) (x[3])))
shared <- intersect(files_celltypes, celltypes)
files <- files[files_celltypes %in% shared]

si <- seqinfo(tmp)

### DOWNLOAD DATA AND READ TO LIST
dir.create(file.path("data", sprintf("epimap-%s-hg19", mark)), recursive=TRUE)

list_of_grs <- lapply(files, function(x){
  ID <- gsub(".sub", "", sapply(strsplit(x, "_", fixed=TRUE), function(x) (x[3])))
  message(ID)
  
  if(sprintf("RAW-%s.bw", ID) %in% list.files(file.path("data", sprintf("epimap-%s-hg19", mark)), pattern = "*.bw")){
    message("Already processed")
    binned_data <- rtracklayer::import.bw(file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("RAW-%s.bw", ID)))
  } else {
    message("Processing")
    # download
    bw <- paste(url, x, sep = "/")
    td <- tempdir(); message(paste(td, x))
    tf <- tempfile(tmpdir = td, fileext = ".bigWig")
    download.file(bw, tf, method='curl')
    
    # to gr
    gr <- rtracklayer::import.bw(tf)
    common_lvls <- gtools::mixedsort(intersect(seqlevels(gr), seqlevels(si)))
    seqlevels(gr, pruning.mode = "coarse") <- common_lvls
    seqinfo(gr) <- si[common_lvls]
    gr <- sort(gr[seqnames(gr) %in% paste0("chr", c(seq(1,22), "X"))])
    
    # bin
    score <- GenomicRanges::coverage(gr, weight = "score")
    binned_data <- GenomicRanges::binnedAverage(bins, score, "score")
    # binned_data$score <- log2(binned_data$score + 1)
    # binned_data
    stopifnot(identical(length(bins), length(binned_data)))
    
    # export and clean up
    out <- file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("RAW-%s.bw", gsub(".sub", "", sapply(strsplit(x, "_", fixed=TRUE), function(x) (x[3])))))
    rtracklayer::export.bw(binned_data, con = out) 
    clean <- sprintf("rm %s", tf); message(clean)
    system(clean)
  }
  
  return(binned_data)
})


### AGGREGRATE 
bw <- list.files(file.path("data", sprintf("epimap-%s-hg19", mark)), pattern = "*.bw", full.names = TRUE)
list_of_grs <- lapply(bw, rtracklayer::import.bw)
names(list_of_grs) <- gsub(".sub", "", sapply(strsplit(tools::file_path_sans_ext(basename(bw)), "-", fixed=TRUE), function(x) (x[2])))

atlas <- granges(list_of_grs[[1]])
atlas_mcols <- do.call(cbind.data.frame, lapply(list_of_grs, mcols))
mcols(atlas) <- atlas_mcols
names(mcols(atlas)) <- names(list_of_grs)
atlas


### SAVE
saveRDS(atlas, compress = TRUE, file.path("data", sprintf("epimap-%s-hg19", mark), sprintf("FULL-ATLAS-EPIMAP-HG19-%s-RAW.rds", mark)))
