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
pacman::p_load(tidyverse, ggplot2, GenomicRanges, caret, foreach, doParallel, data.table, rtracklayer) # parallel

parser <- OptionParser()
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("--train_chr"), default="all",
              help = "Which chromosome used for training (\"chr1\"-\"chrX\", or \"all\" if use all chromosomes) [default \"%default\"]",
              metavar="train chromosome"),
  make_option(c("--test_chr"), default="all",
              help = "Which chromosome to use for testing (\"chr1\"-\"chrX\", or \"all\" if use all chromosomes) [default \"%default\"]",
              metavar="test chromosome"),
  make_option(c("--num_nearest_celltypes"), type="integer", default=5, help="Number of k nearest training cell types to look for [default %default]", metavar="num nearest celltypes"),
  make_option(c("--num_nearest_gene"), type="integer", default=5, help="Number of nearest genes to build model [default %default]", metavar="num nearest gene"),
  make_option(c("--training_ct"), default="./path_to_training_ct/training_ct.txt",
              help = "A file listing celltypes to use for training [default \"%default\"]",
              metavar="training celltypes"),
  make_option(c("--testing_ct"), default="./path_to_testing_ct/testing_ct.txt",
              help = "A file listing celltypes to use for testing [default \"%default\"]",
              metavar="testing celltypes"),
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c( "--gene_train"), default="./data/atlas_gene.rds",
              help = "Set file path for Gene Expression data used for training [default \"%default\"]",
              metavar="training gene expression"),
  make_option(c("--gene_test"), default="./data/atlas_gene.rds",
              help = "Set file path for Gene Expression data used for testing [default \"%default\"]",
              metavar="testing gene expression"),
  make_option(c("-c", "--chip"), default="./data/chip_dir/",
              help = "Set directory for ChIP-seq data [default \"%default\"]",
              metavar="chip-seq data directory"),
  make_option(c("-m", "--meta"), default="./data/atlas_gene_meta.rds",
              help = "Set file path for ChIP-seq data [default \"%default\"]",
              metavar="gene metadata"),
  make_option(c("--mark", default="H3K27ac", help = "The mark that is used for analysis", metavar="histone mark")),
  make_option(c("--proj_name_train"), default="new_proj_train",
              help = "Define a project name. Trained models and training features will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name train"),
  make_option(c("--proj_name_test"), default="new_proj_test",
              help = "Define a project name. Prediction results for testing cell types will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name test"),
  make_option(c("-r", "--outdir"), default="./predictions/across_model/",
              help = "Define the output directory (where most across-model predictions will be in)",
              metavar="output directory"),
  make_option(c("--seed"), default=42,
              help = "Set seed for this experiment (such as training-validation split)",
              metavar="seed"),
  make_option(c("--dist"), default="corr",
              help = "Method to use to calculate nearest training cell-types. 'corr' for Pearson Correlation (default) and 'euc' for Euclidean distance)",
              metavar="dist"),
  make_option(c("--global"), action='store_true', default=FALSE,
              help = "Use --global when you only want to use global features (will be fast). Otherwise both global and local features will be used.", 
              metavar="global"), ## NOT YET IMPLEMENTED
  make_option(c("--leave_one_out"), action='store_true', default=FALSE,
              help = "Use --leave-one-out when you are working on leave one out experiments (train with all other cell types and test with one cell type). Default is False."),
  make_option(c("--task_id"), default=1, help = "Task ID for this experiment", metavar = "task id"),
  make_option(c("--sample_size", default=100000, help = "Sample size for data used for model training", metavar="sample size")),
  make_option(c("--which_chunk", default=1, help = "Which chunk to predict in testing", metavar="which chunk")),
  make_option(c("--chunk_size", default=100000, help = "Chunk size to predict in testing", metavar="chunk size")),
  make_option(c("--save_intermediate"), action='store_true', default=FALSE,
              help = "Use this flag if you want to save all intermediate files. Otherwise only the models will be saved and other intermediate files will be deleted.")

  ####### TODO #########
  ### make option knn ###
)
# opt$proj_name <- "090822_test"

opt <- parse_args(OptionParser(option_list=option_list))
if(opt$verbose){
  print(opt)
}

# opt$proj_name_train <- "epimap-H3K4me1-hg19-non-dup_train"
# opt$proj_name_test <- "epimap-H3K4me1-hg19-non-dup_test"
# opt$task_id <- 2
# opt <- readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s/opt_handles.rds", opt$proj_name_test, opt$task_id))

dir.create(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s/", opt$proj_name_test, opt$task_id), recursive=TRUE)
saveRDS(opt, sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s/opt_handles.rds", opt$proj_name_test, opt$task_id))

opt$proj_name_train <- "epimap-H3K27ac-hg19-non-dup_train/epimap-H3K27ac-hg19-non-dup_test_BSS00007"
opt$proj_name_test <- "epimap-H3K27ac-hg19-non-dup_test/epimap-H3K27ac-hg19-non-dup_test_BSS00007"
opt$task_id <- 1
opt <- readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s/opt_handles.rds", opt$proj_name_test, opt$task_id))

# opt$num_chunks <- 13

k_num_train <- as.numeric(opt$num_nearest_celltypes)
knn <- as.numeric(opt$num_nearest_gene)

wd <- opt$wd
setwd(wd)
source(file.path("batch", "workflow", "utils", "utils.R"))
source(file.path("batch", "workflow", "ChromExp", "model_helper_functions.R"))

train_chr <- opt$train_chr
test_chr <- opt$test_chr
mark <- opt$mark
save_intermediate <- opt$save_intermediate

k <- opt$k

ATLAS_GENE_METADATA <- try(readRDS(opt$meta))
if (class(ATLAS_GENE_METADATA) == "try-error"){
  message("Metadata file not found: Please input a valid gene metadata file location! ")
  ATLAS_GENE_METADATA
}
message("Data: ATLAS_GENE_METADATA")

ATLAS_GENE_train <- try(readRDS(opt$gene_train))
if (class(ATLAS_GENE_train) == "try-error"){
  message("Training gene expression file not found: Please input a valid training gene expression file location! ")
  ATLAS_GENE_train
}
message("Data: ATLAS_GENE_train")
colnames(ATLAS_GENE_train) <- unlist(lapply(colnames(ATLAS_GENE_train), function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x
}))

ATLAS_GENE_test <- try(readRDS(opt$gene_test))
if (class(ATLAS_GENE_test) == "try-error"){
  message("Testing Gene expression file not found: Please input a valid testing gene expression file location! ")
  ATLAS_GENE_test
}
message("Data: ATLAS_GENE_test")
colnames(ATLAS_GENE_test) <- unlist(lapply(colnames(ATLAS_GENE_test), function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x
}))

training_celltypes <- try(readRDS(opt$training_ct))
if (class(training_celltypes) == "try-error"){
  message("Training celltypes file not found!!")
  training_celltypes
}
training_celltypes <- as.character(training_celltypes)
training_celltypes <- unlist(lapply(training_celltypes, function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x
}))

loo <- opt$leave_one_out

if(loo){
  testing_celltypes <- opt$testing_ct
  testing_celltypes = str_replace_all(testing_celltypes, " ", ".")
  testing_celltypes = str_replace_all(testing_celltypes, "-", ".")
} else {
  testing_celltypes <- try(readRDS(opt$testing_ct))
  if (class(testing_celltypes) == "try-error"){
    message("Held-out celltypes file not found!!")
    testing_celltypes
  }
  testing_celltypes <- as.character(testing_celltypes)
  testing_celltypes <- unlist(lapply(testing_celltypes, function(x){
    x = str_replace_all(x, " ", ".")
    x = str_replace_all(x, "-", ".")
    x
  }))
}

if(loo){
  all_celltypes <- training_celltypes
  training_celltypes <- setdiff(training_celltypes, testing_celltypes)
}

which_chunk = as.numeric(opt$which_chunk)
sample_size = as.numeric(opt$sample_size)
chunk_size = as.numeric(opt$chunk_size)

test_ct_label = testing_celltypes[1]
num_test_ct = length(testing_celltypes)
num_train_ct = length(training_celltypes)

# chip_path = "/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-chromimpute/converted_tracks/"

proj_name_train <- opt$proj_name_train
proj_name_test <- opt$proj_name_test
across_dir <- opt$outdir
seed <- opt$seed
dist_method <- opt$dist
global <- opt$global

stopifnot(identical(row.names(ATLAS_GENE_test), ATLAS_GENE_METADATA$gene_id) || identical(row.names(ATLAS_GENE_test), ATLAS_GENE_METADATA$gene_name))
stopifnot(identical(row.names(ATLAS_GENE_train), ATLAS_GENE_METADATA$gene_id) || identical(row.names(ATLAS_GENE_test), ATLAS_GENE_METADATA$gene_name))

row.names(ATLAS_GENE_test) <- ATLAS_GENE_METADATA$gene_id
row.names(ATLAS_GENE_train) <- ATLAS_GENE_METADATA$gene_id

common_genes <- intersect(intersect(mcols(ATLAS_GENE_METADATA)$gene_id, row.names(ATLAS_GENE_test)), row.names(ATLAS_GENE_train))

ATLAS_GENE_test <- ATLAS_GENE_test[common_genes,]
ATLAS_GENE_train <- ATLAS_GENE_train[common_genes,]

X_train <- ATLAS_GENE_train[, training_celltypes, drop = FALSE]
X_test <- ATLAS_GENE_test[, testing_celltypes, drop = FALSE]

X_train_st = scale(X_train)
X_test_st = scale(X_test)
dist_test = compute_euc_distance(X_train_st, X_test_st)
corr_test = compute_corr(X_train_st, X_test_st)

# Y_train = readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("ATLAS_CHIP_train_%s.rds", sample_size)))

### READ IN CHIP SEQ DATA FROM FOLDER
## Generate GR first
ct = training_celltypes[1]
chip_ct = import.bw(paste0(opt$chip,sprintf('%s-%s-HG19-%s-RAW.bw', ct, test_chr, mark)))

## Now generate chip seq data for a set of sampled positions
set.seed(seed)
if (which_chunk*chunk_size > nrow(mcols(chip_ct))){
  read_pos = ((which_chunk-1)*chunk_size+1):(nrow(mcols(chip_ct)))
} else {
  read_pos = ((which_chunk-1)*chunk_size+1):((which_chunk)*chunk_size)
}

# dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "sampled_pos"))
# saveRDS(sampled_pos, file.path(across_dir, proj_name_train, "intermediate-files", "sampled_pos", sprintf("%s_sampled_pos.rds", sample_size)))

chip_ct_chunk = chip_ct[read_pos,]
chip_ct_gr = granges(chip_ct_chunk)
chip_ct_gr_train = chip_ct_gr

# chip_ct_gr = granges(chip_ct)

### TODO
for (ct in testing_celltypes){
 chip_ct = import.bw(paste0(opt$chip,sprintf('%s-%s-HG19-%s-RAW.bw', ct, test_chr, mark)))
 names(mcols(chip_ct)) <- ct
 chip_ct_chunk <- chip_ct[read_pos,]
 mcols(chip_ct_gr) = cbind(as.data.frame(mcols(chip_ct_gr)), as.data.frame(mcols(chip_ct_chunk[,ct])))
}

ATLAS_CHIP <- chip_ct_gr
k_list <- seq(1,k_num_train)

###
for (ct in training_celltypes){
 chip_ct = import.bw(paste0(opt$chip,sprintf('%s-%s-HG19-%s-RAW.bw', ct, test_chr, mark)))
 names(mcols(chip_ct)) <- ct
 chip_ct_chunk <- chip_ct[read_pos,]
 mcols(chip_ct_gr_train) = cbind(as.data.frame(mcols(chip_ct_gr_train)), as.data.frame(mcols(chip_ct_chunk[,ct])))
}

Y_train <- chip_ct_gr_train

#### GENERATE GLOBAL KNN FEATURES FROM EUCLIDEAN DIST
if (dist_method == "euc"){
  for (test_ct in testing_celltypes){
    first_k_train_ct = names(sort(dist_test[,test_ct]))[1:k_num_train]
    knn_features_all = generate_global_avg(chip=ATLAS_CHIP, k_list=k_list, first_k_train_ct=first_k_train_ct, y_train=Y_train)
    dir.create(file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr), recursive=TRUE)
    saveRDS(knn_features_all, compress=TRUE, file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr, sprintf("%s-%s-%s-avg-knn-tracks.rds", test_ct, which_chunk, chunk_size)))
  }
}

#### GENERATE GLOBAL KNN FEATURES FROM PEARSON CORRELATION
if (dist_method == "corr"){
  for (test_ct in testing_celltypes){
    if(!file.exists(file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr, sprintf("%s-%s-%s-avg-knn-tracks.rds", test_ct, which_chunk, chunk_size)))){
      first_k_train_ct = names(sort(corr_test[,test_ct], decreasing = TRUE, na.last = TRUE))[1:k_num_train]
      knn_features_all = generate_global_avg(chip=ATLAS_CHIP, k_list=k_list, first_k_train_ct=first_k_train_ct, y_train=Y_train)
      dir.create(file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr), recursive=TRUE)
      saveRDS(knn_features_all, compress=TRUE, file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr, sprintf("%s-%s-%s-avg-knn-tracks.rds", test_ct, which_chunk, chunk_size)))
    }
  }
}

### COMBINE GLOBAL FEATURES FROM DIFFERENT CELLTYPES INTO A LIST
knn_global_test_all = lapply(testing_celltypes, function(test_ct){
  knn_global_feature = readRDS(file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr, sprintf("%s-%s-%s-avg-knn-tracks.rds", test_ct, which_chunk, chunk_size)))
  knn_global_feature
})
names(knn_global_test_all) <- testing_celltypes


if(!file.exists(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", sprintf("nearest_celltypes_test_l_%s_%s.rds", test_ct_label, num_test_ct)))){
  nearest_celltypes_test_l = list()

  # knn_features_all_test = readRDS(file.path(across_dir, proj_name, "intermediate-files", "avg-tracks-knn-test", chr, sprintf("%s-avg-knn-tracks.rds", test_ct)))
  for (gene in row.names(ATLAS_GENE_train)){
    test_df = data.frame(matrix(nrow=k_num_train, ncol=length(testing_celltypes)))
    colnames(test_df) = testing_celltypes

    if (sd(ATLAS_GENE_train[gene,]) > 0){  ## Could tune this parameter
      for (ct in testing_celltypes){
        test_df[,ct] <- names(sort(abs(X_train_st[gene,] - X_test_st[gene,ct])))[1:k_num_train]
      }
      nearest_celltypes_test_l[[gene]] <- test_df
    }
  }
  dir.create(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results"), recursive=TRUE)
  saveRDS(nearest_celltypes_test_l, file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", "nearest_celltypes_test_l.rds"))
}

message("Nearest cell types calculated!")
nearest_celltypes_test_l <- readRDS(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", "nearest_celltypes_test_l.rds"))

ct <- training_celltypes[1]
if(!file.exists(file.path(across_dir, proj_name_test, "intermediate-files", sprintf("features_%s_%s_%s", test_chr, which_chunk, chunk_size), test_chr, "nearest_gene_features.rds"))){
  feature_ct <- feature_eng_train(signal = Y_train, 
                    expr = X_train_st, 
                    expr_metadata = ATLAS_GENE_METADATA, 
                    celltype = ct, 
                    k_nearest = 2, 
                    normalize=FALSE)
  dir.create(file.path(across_dir, proj_name_test, "intermediate-files", sprintf("features_%s_%s_%s", test_chr, which_chunk, chunk_size), test_chr), recursive=TRUE)
  saveRDS(feature_ct, compress=TRUE, file.path(across_dir, proj_name_test, "intermediate-files", sprintf("features_%s_%s_%s", test_chr, which_chunk, chunk_size), test_chr, "nearest_gene_features.rds"))
}

feature_ct <- readRDS(file.path(across_dir, proj_name_test, "intermediate-files", sprintf("features_%s_%s_%s", test_chr, which_chunk, chunk_size), test_chr, "nearest_gene_features.rds"))
feature_ct

Y_train_df = as.data.frame(mcols(Y_train))

if(!file.exists(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr, sprintf("knn_local_features_test_ct_list_%s_%s.rds", which_chunk, chunk_size)))){
  knn_local_features_test_ct_list <- lapply(testing_celltypes, function(test_ct){
    # which_feature_uniq <- 1
    temp_list=list()
    print(test_ct)

    local1_all_features <- lapply(unique(feature_ct$Gene1), function(gene1){
      local1_features <- generate_local_features(gene=gene1,
                          gene_col="Gene1",
                          feature=feature_ct,
                          nearest_celltype_list=nearest_celltypes_test_l,
                          train_ct=test_ct,
                          k_list=k_list,
                          y_train_df=Y_train_df,
                          new_col_name="local1_")
      local1_features
    })
    local1_all_features_df = do.call(rbind, local1_all_features)
    local1_all_features_df <- local1_all_features_df[order(local1_all_features_df$rows),]
    
    local2_all_features <- lapply(unique(feature_ct$Gene2), function(gene2){
      local2_features <- generate_local_features(gene=gene2,
                          gene_col="Gene2",
                          feature=feature_ct,
                          nearest_celltype_list=nearest_celltypes_test_l,
                          train_ct=test_ct,
                          k_list=k_list,
                          y_train_df=Y_train_df,
                          new_col_name="local2_")
      local2_features
    })
    local2_all_features_df = do.call(rbind, local2_all_features)
    local2_all_features_df <- local2_all_features_df[order(local2_all_features_df$rows),]
    
    local_all_features <- merge(local1_all_features_df, local2_all_features_df, by="rows")
    local_all_features
  })
  dir.create(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr), recursive=TRUE)
  saveRDS(knn_local_features_test_ct_list, file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr, sprintf("knn_local_features_test_ct_list_%s_%s.rds", which_chunk, chunk_size)))
}

print("Local Features generated!!")

### COMBINE LOCAL AND GLOBAL FEATURES
knn_local_features_test_ct_list <- readRDS(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr, sprintf("knn_local_features_test_ct_list_%s_%s.rds", which_chunk, chunk_size)))
names(knn_local_features_test_ct_list) <- testing_celltypes

k_nearest_dist = 1 ## Hardcoded parameter value. Needs to be changed.

knn_features_all_ct_test_df <- data.frame(matrix(nrow=0,ncol=2*length(k_list)+2*k_nearest_dist+2))
colnames(knn_features_all_ct_test_df) <- c("Cell_type", c(paste0("local1_", k_list), paste0("local2_", k_list)), paste0("Expr", seq(1, k_nearest_dist, 1)), paste0("Dist", seq(1, k_nearest_dist, 1)), "locus")

# # ### Generate testing cell types' predictions
ptm <- proc.time()

k_list <- seq(1,k_num_train)
k_nearest_dist = 1 ## Hardcoded parameter value. Needs to be changed.

# knn_features_all_ct_test_df <- data.frame(matrix(nrow=0,ncol=2*length(k_list)+2*k_nearest_dist+2))
# colnames(knn_features_all_ct_test_df) <- c("Cell_type", c(paste0("local1_", k_list), paste0("local2_", k_list)), paste0("Expr", seq(1, k_nearest_dist, 1)), paste0("Dist", seq(1, k_nearest_dist, 1)), "locus")

if(!file.exists(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr, sprintf("knn_all_features_test_%s_%s.rds", which_chunk, chunk_size)))){
  for (test_ct in testing_celltypes){
    # test_ct <- validation_celltypes[[1]]
    # test_ct=testing_celltypes[1]
    knn_features_test_ss <- knn_local_features_test_ct_list[[test_ct]]
    knn_features_test_ss_df <- knn_features_test_ss[, c(paste0("local1_", k_list), paste0("local2_", k_list))]
    knn_features_test_ss_df$locus <- knn_features_test_ss$rows
    dist_exp_feature_test <- feature_eng_dist_exp_test(signal=ATLAS_CHIP, expr=ATLAS_GENE_test, expr_metadata=ATLAS_GENE_METADATA, celltype=test_ct, k_nearest = k_nearest_dist, normalize=FALSE)
    dist_exp_feature_df_test <- as.data.frame(mcols(dist_exp_feature_test)[,-1])
    names(dist_exp_feature_df_test) <- c(paste0("Expr", seq(1, k_nearest_dist, 1)), paste0("Dist", seq(1, k_nearest_dist, 1)))

    knn_features_test_ss_df <- cbind(dist_exp_feature_df_test, knn_features_test_ss_df)
    knn_features_test_ss_df$Cell_type = test_ct
    knn_features_test_ss_df <- knn_features_test_ss_df[,c("Cell_type", c(paste0("local1_", k_list), paste0("local2_", k_list)), paste0("Expr",seq(1,k_nearest_dist)), paste0("Dist",seq(1,k_nearest_dist)), "locus")]

    knn_features_all_ct_test_df <- rbind(knn_features_all_ct_test_df, knn_features_test_ss_df)
  }

  for(test_ct in testing_celltypes){
    knn_features_all = readRDS(file.path(across_dir, proj_name_test, "intermediate-files", "avg-tracks-knn-test", test_chr, sprintf("%s-%s-%s-avg-knn-tracks.rds", test_ct, which_chunk, chunk_size)))
    knn_features_all_df = as.data.frame(mcols(knn_features_all))
    colnames(knn_features_all_df) = paste0("avg_", k_list)
    knn_features_all_ct_test_df[which(knn_features_all_ct_test_df$Cell_type == test_ct),paste0("global_", k_list)] <- knn_features_all_df
  }
  dir.create(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr), recursive=TRUE)
  saveRDS(knn_features_all_ct_test_df, file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr, sprintf("knn_all_features_test_%s_%s.rds", which_chunk, chunk_size)))
}

knn_features_all_ct_test_df <- readRDS(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", test_chr, sprintf("knn_all_features_test_%s_%s.rds", which_chunk, chunk_size)))
# knn_features_all_ct_test_df <- readRDS(file.path(across_dir, proj_name_test, "intermediate-files", "intermediate-results", chr, sprintf("knn_all_features_test_%s_%s.rds", which_chunk, chunk_size)))

# ### DEBUG CODE: CHECK FOR DIFFERENT CELL TYPES WHETHER THEY ARE SHARING SAME FEATURES
# check_ct = unique(knn_features_all_ct_test_df$Cell_type)
# head(knn_features_all_ct_test_df[knn_features_all_ct_test_df$Cell_type == "BSS00007",])
# head(knn_features_all_ct_test_df[knn_features_all_ct_test_df$Cell_type == "BSS00043",])

training_pred_lm_full_df <- data.frame(matrix(nrow=nrow(knn_features_all_ct_test_df), ncol=length(training_celltypes)))
colnames(training_pred_lm_full_df) <- training_celltypes
training_pred_lm_global_df <- training_pred_lm_full_df
training_pred_reg_full_df <- training_pred_lm_global_df
training_pred_reg_global_df <- training_pred_reg_full_df

for (train_ct in training_celltypes){
  if (sample_size == 0){
    pruned_regression_tree_full <- readRDS(file.path(across_dir, proj_name_train, "models", "ensmbl-model", train_chr, sprintf("reg_model_full_val_%s.rds", train_ct)))
    pruned_regression_tree_global <- readRDS(file.path(across_dir, proj_name_train, "models", "ensmbl-model", train_chr, sprintf("reg_model_global_val_%s.rds", train_ct))) 
  } else {
    pruned_regression_tree_full <- readRDS(file.path(across_dir, proj_name_train, "models", "ensmbl-model", train_chr, sprintf("reg_model_full_val_%s_%s.rds", train_ct, sample_size)))
    pruned_regression_tree_global <- readRDS(file.path(across_dir, proj_name_train, "models", "ensmbl-model", train_chr, sprintf("reg_model_global_val_%s_%s.rds", train_ct, sample_size)))
  }

  training_pred_reg_full <- predict(pruned_regression_tree_full, knn_features_all_ct_test_df)
  training_pred_reg_global <- predict(pruned_regression_tree_global, knn_features_all_ct_test_df)

  training_pred_reg_full_df[,train_ct] <- training_pred_reg_full
  training_pred_reg_global_df[,train_ct] <- training_pred_reg_global
}

dir.create(file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr), recursive=TRUE)
saveRDS(training_pred_reg_full_df, file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr, sprintf("training_pred_df_reg_full_%s_%s.rds", which_chunk, chunk_size)))
saveRDS(training_pred_reg_global_df, file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr, sprintf("training_pred_df_reg_global_%s_%s.rds", which_chunk, chunk_size)))

training_pred_reg_full_df <- readRDS(file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr, sprintf("training_pred_df_reg_full_%s_%s.rds", which_chunk, chunk_size)))
training_pred_reg_global_df <- readRDS(file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr, sprintf("training_pred_df_reg_global_%s_%s.rds", which_chunk, chunk_size)))

# training_pred_reg_full_df$average = rowSums(training_pred_reg_full_df)/length(training_celltypes)
# training_pred_reg_global_df$average = rowSums(training_pred_reg_global_df)/length(training_celltypes)

training_pred_reg_full_df$average = rowMeans(training_pred_reg_full_df)
training_pred_reg_global_df$average = rowMeans(training_pred_reg_global_df)

knn_features_all_ct_test_df$reg_pred_full <- training_pred_reg_full_df$average
knn_features_all_ct_test_df$reg_pred_global <- training_pred_reg_global_df$average

nrow_test_pred = nrow(knn_features_all_ct_test_df[knn_features_all_ct_test_df$Cell_type == testing_celltypes[1],])
# knn_features_all_ct_test_df_ss = knn_features_all_ct_test_df[knn_features_all_ct_test_df$Cell_type == test_ct & knn_features_all_ct_test_df$locus == 1,]

testing_pred_df_reg_full <- data.frame(matrix(nrow=nrow_test_pred, ncol=length(testing_celltypes)))
colnames(testing_pred_df_reg_full) <- testing_celltypes

for(test_ct in testing_celltypes){
  print(test_ct)
  pred_signal = knn_features_all_ct_test_df$reg_pred_full[knn_features_all_ct_test_df$Cell_type == test_ct]
  testing_pred_df_reg_full[,test_ct] = pred_signal
}

testing_pred_df_reg_global <- data.frame(matrix(nrow=nrow_test_pred, ncol=length(testing_celltypes)))
colnames(testing_pred_df_reg_global) <- testing_celltypes

for(test_ct in testing_celltypes){
  pred_signal = knn_features_all_ct_test_df$reg_pred_global[knn_features_all_ct_test_df$Cell_type == test_ct]
  testing_pred_df_reg_global[,test_ct] = pred_signal
}

dir.create(file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr), recursive = TRUE)
saveRDS(testing_pred_df_reg_full, file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr, sprintf("testing_pred_df_reg_full_%s_%s.rds", which_chunk, chunk_size)))
saveRDS(testing_pred_df_reg_global, file.path(across_dir, proj_name_test, "predictions", "ensmbl-model", test_chr, sprintf("testing_pred_df_reg_global_%s_%s.rds", which_chunk, chunk_size)))
print("Predictions generated and saved!")
print(proc.time()-ptm)


### Remove intermediate files after all regression trees are built
if (!save_intermediate){
  unlink(file.path(across_dir, proj_name_test, "intermediate-files", sprintf("features_%s_%s_%s", test_chr, which_chunk, chunk_size)), recursive = TRUE, force = TRUE)
}

