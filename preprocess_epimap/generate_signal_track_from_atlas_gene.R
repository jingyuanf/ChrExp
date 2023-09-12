#!/usr/bin/env Rscript

### Environment ###
rm(list=ls(all=TRUE))
options(warn=-1)
options(scipen=999)
setwd("./")
source(file.path("batch", "workflow", "utils", "utils.R"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

pacman::p_load(tidyverse, GenomicRanges, plyr, dplyr) # caret, doParallel
#library(biganalytics)

parser <- OptionParser()
option_list <- list(
  make_option(c("--tss"), action="store_true", default=FALSE,
              help="Generate tracks annotated in only TSS region (Otherwise generate tracks annotated in all gene bodies)"),
  make_option(c("--ct"), default="BSS00007", 
              help="Save track for which cell type")
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))

tss <- opt$tss
ct <- opt$ct

### HUMAN EPIMAP DATA
## First get gene expression data

expr <- readRDS("./data/epimap-gene-exp/ATLAS_GENE_train.rds")
celltypes <- colnames(expr)

if (!file.exists("./data/epimap-gene-exp/gene_exp_celltypes.txt")){
    saveRDS(celltypes, "./data/epimap-gene-exp/gene_exp_celltypes.rds")
    writeLines(celltypes, "./data/epimap-gene-exp/gene_exp_celltypes.txt")
}

expr_meta <- readRDS("./data/epimap-gene-exp/ATLAS_GENE_METADATA.rds")

if (!file.exists("./data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")){
    # Then get ChIP-seq data to get track format (Here we use H3K27ac as an example)
    chip <- readRDS("./data/epimap-H3K9ac-hg19/FULL-ATLAS-EPIMAP-HG19-H3K9ac-RAW.rds")
    chip <- granges(chip)

    # Save GRANGES ChIP-seq data 
    dir.create("./data/epimap-hg19/", recursive=TRUE)
    saveRDS(chip_gr, "./data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")
}

## Directly read GRANGES ChIP-seq data
if (tss){
    if (!file.exists("./data/epimap-gene-exp/ATLAS_GENE_EXP_TRACK_TSS.rds")){
        chip <- readRDS("./data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")

        expr_meta_tss <- reduce_to_tss(expr_meta)
        overlaps <- findOverlaps(chip, expr_meta_tss)
        identical(expr_meta$gene_name, rownames(expr))

        df = data.frame(matrix(nrow=15181508, ncol=ncol(expr)))
        colnames(df) <- colnames(expr)
        df[is.na(df)] = 0
        df[queryHits(overlaps), ] <- expr[subjectHits(overlaps), ]

        mcols(chip) = df
        saveRDS(chip, "./data/epimap-gene-exp/ATLAS_GENE_EXP_TRACK_TSS.rds")
    }
}

if (!tss){
    if(!file.exists("./data/epimap-gene-exp/ATLAS_GENE_EXP_TRACK_GENEBODY.rds")){
        chip <- readRDS("./data/epimap-hg19/FULL-ATLAS-EPIMAP-HG19-GRANGES.rds")
        
        overlaps <- findOverlaps(chip, expr_meta)
        identical(expr_meta$gene_name, rownames(expr))

        df = data.frame(matrix(nrow=15181508, ncol=ncol(expr)))
        colnames(df) <- colnames(expr)
        df[is.na(df)] = 0
        df[queryHits(overlaps), ] <- expr[subjectHits(overlaps), ]

        mcols(chip) = df
        saveRDS(chip, "./data/epimap-gene-exp/ATLAS_GENE_EXP_TRACK_GENEBODY.rds")
    }
}



## Generate gene exp track annotated in TSS
if (tss){
    chip <- readRDS("./data/epimap-gene-exp/ATLAS_GENE_EXP_TRACK_TSS.rds")
    print(ct)
    filedir = "./data/epimap-gene-exp/raw-tracks-tss/"
    dir.create(filedir, recursive=TRUE)

    con_expr = sprintf("%s_rna_epimap.bedgraph", ct)

    expr_track = chip[,ct]
    mcols(expr_track)[is.na(mcols(expr_track))] = 0
    colnames(mcols(expr_track)) = "score"
    rtracklayer::export(expr_track, con=file.path(filedir, con_expr), format="bedGraph")

} else {
    chip <- readRDS("./data/epimap-gene-exp/ATLAS_GENE_EXP_TRACK_GENEBODY.rds")
    print(ct)
    filedir = "./data/epimap-gene-exp/raw-tracks-genebody/"
    dir.create(filedir, recursive=TRUE)

    con_expr = sprintf("%s_rna_epimap.bedgraph", ct)

    expr_track = chip[,ct]
    mcols(expr_track)[is.na(mcols(expr_track))] = 0
    colnames(mcols(expr_track)) = "score"
    rtracklayer::export(expr_track, con=file.path(filedir, con_expr), format="bedGraph")
}

