#!/usr/bin/env Rscript

### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
getwd()
setwd("./")
source(file.path("batch", "workflow", "utils", "utils.R"))

### LIBS (GET PARAMETERS)
pacman::p_load(tidyverse, ggplot2, GenomicRanges) # parallel

ATLAS_GENE_META <- rtracklayer::import("./data/epimap-metadata/gencode.v19.annotation.gtf.gz")

## Subset Meta data to only genes
ATLAS_GENE_META_genes <- ATLAS_GENE_META[ATLAS_GENE_META$type == "gene",]

gene_type_list = c("protein_coding")
chr_list = paste0("chr", c(seq(1,22),"X"))

## Subset Meta data to protein coding genes and chr1-22 and chrX
ATLAS_GENE_META_pc <- ATLAS_GENE_META_genes[ATLAS_GENE_META_genes$gene_type %in% gene_type_list,]
ATLAS_GENE_META_pc <- ATLAS_GENE_META_pc[seqnames(ATLAS_GENE_META_pc) %in% chr_list,]

## Subset Meta data to gene status of KNOWN (instead of NOVEL or PUTATIVE)
ATLAS_GENE_META_known <- ATLAS_GENE_META_pc[ATLAS_GENE_META_pc$gene_status == "KNOWN", ]
ATLAS_GENE_META_tss <- reduce_to_tss(ATLAS_GENE_META_known)

## Save Meta data
saveRDS(ATLAS_GENE_META_known, "./data/epimap-metadata/gencode_gene_metadata.rds")
saveRDS(ATLAS_GENE_META_tss, "./data/epimap-metadata/gencode_gene_metadata_tss.rds")


# ATLAS_GENE_META_sel <- readRDS("./data/Gene_annotation/gencode_v19_annotation_meta_protein_coding.rds")
# old_ATLAS_GENE_META <- readRDS("/u/home/f/fujy2038/project-ernst/Project/data/preprocessed_files/ATLAS_GENE_METADATA.rds")
# unique(old_ATLAS_GENE_META$gene_type)
