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
  make_option(c("--chr"), default="all",
              help = "Which chromosome to use (\"chr1\"-\"chrX\", or \"all\" if use all chromosomes) [default \"%default\"]",
              metavar="chromosome"),
  make_option(c("--num_nearest_celltypes"), type="integer", default=5, help="Number of k nearest training cell types to look for [default %default]", metavar="num nearest celltypes"),
  make_option(c("--num_nearest_gene"), type="integer", default=5, help="Number of nearest genes to build model [default %default]", metavar="num nearest gene"),
  make_option(c("--training_ct"), default="./path_to_training_ct/training_ct.txt",
              help = "A file listing celltypes to use for training [default \"%default\"]",
              metavar="training celltypes"),
  make_option(c("--holdout"), default="None",
              help = "Held out cell type (or an rds file specifying a list of held out cell types). If held out cell type in the list of training celltypes, remove it from training.",
              metavar="holdout"),
  make_option(c("--wd"), default="./path_to_working_directory",
              help = "Set working directory [default \"%default\"]",
              metavar="working directory"),
  make_option(c( "--gene_train"), default="./data/atlas_gene.rds",
              help = "Set file path for Gene Expression data used for training [default \"%default\"]",
              metavar="training gene expression"),
  make_option(c("-c", "--chip"), default="./data/chip_dir/",
              help = "Set directory for ChIP-seq data [default \"%default\"]",
              metavar="chip-seq data directory"),
  make_option(c("-m", "--meta"), default="./data/atlas_gene_meta.rds",
              help = "Set file path for ChIP-seq data [default \"%default\"]",
              metavar="gene metadata"),
  make_option(c("--mark"), default="H3K27ac", help = "The mark that is used for analysis", metavar="histone mark"),
  make_option(c("--proj_name_train"), default="new_proj_train",
              help = "Define a project name. Trained models and training features will be stored in a folder named by this project name, inside a larger folder defined by OUTPUT DIRECTORY.",
              metavar="project name train"),
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
  make_option(c("--sample_size"), default=100000, help = "Sample size for data reading into the model", metavar="sample size"),
  make_option(c("--task_id"), default=1, help = "Task ID for this experiment", metavar = "task id"),
  make_option(c("--qn"), action='store_true', default=FALSE,
              help = "Use --qn when you are working on quantile normalized data. Otherwise don't use this flag. Default is False (raw data)."),
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

dir.create(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s", opt$proj_name_train, opt$task_id), recursive=TRUE)
saveRDS(opt, sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s/opt_handles.rds", opt$proj_name_train, opt$task_id))

# opt$proj_name_train <- "epimap-H3K9me3-hg19-non-dup_train"
# opt$task_id <- 1
# opt <- readRDS(sprintf("/u/home/f/fujy2038/project-ernst/Project/temp_file/%s/%s/opt_handles.rds", opt$proj_name_train, opt$task_id))

### TAKE INPUT FROM ARGUMENT HANDLES

k_num_train <- as.numeric(opt$num_nearest_celltypes)
knn <- as.numeric(opt$num_nearest_gene)

wd <- opt$wd
setwd(wd)
source(file.path("batch", "workflow", "utils", "utils.R"))
source(file.path("batch", "workflow", "ChromExp", "model_helper_functions.R"))

num_train_ct <- opt$num_train_ct
chr <- opt$chr
k <- opt$k
mark <- opt$mark
holdout_f <- opt$holdout
save_intermediate <- opt$save_intermediate
# loo <- opt$leave_one_out


training_ct <- try(readRDS(opt$training_ct))
if (class(training_ct) == "try-error"){
  message("Training celltypes file not found!!")
  training_ct
}
training_ct <- as.character(training_ct)
training_ct <- unlist(lapply(training_ct, function(x){
  x = str_replace_all(x, " ", ".")
  x = str_replace_all(x, "-", ".")
  x
}))

if (holdout_f != "None"){
  if(!endsWith(holdout_f, ".rds")){
    holdout <- holdout_f
    holdout = str_replace_all(holdout, " ", ".")
    holdout = str_replace_all(holdout, "-", ".")
  } else {
    holdout <- try(readRDS(holdout_f))
    if (class(holdout) == "try-error"){
      message("Held-out celltypes file not found!!")
      holdout
    }
    holdout <- as.character(holdout)
    holdout <- unlist(lapply(holdout, function(x){
      x = str_replace_all(x, " ", ".")
      x = str_replace_all(x, "-", ".")
      x
    }))
  }
}


if(holdout != "None"){
  all_ct <- training_ct
  training_ct <- setdiff(training_ct, holdout)
}

num_train_ct <- length(training_ct)

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

# chip_path = opt$chip

ATLAS_GENE_METADATA <- try(readRDS(opt$meta))
if (class(ATLAS_GENE_METADATA) == "try-error"){
  message("Metadata file not found: Please input a valid gene metadata file location! ")
  ATLAS_GENE_METADATA
}
message("Data: ATLAS_GENE_METADATA")

proj_name_train <- opt$proj_name_train
across_dir <- opt$outdir
seed <- opt$seed
dist_method <- opt$dist
global <- opt$global


### SET WORKING DIRECTORY AND FORMAT INPUT DATA
chip_path = opt$chip

# training_celltypes <- intersect(training_ct, intersect(colnames(ATLAS_GENE_train), names(mcols(ATLAS_CHIP))))
training_celltypes <- intersect(training_ct, colnames(ATLAS_GENE_train))
names(ATLAS_GENE_METADATA) <- ATLAS_GENE_METADATA$gene_id
ATLAS_GENE_train <- ATLAS_GENE_train[,training_celltypes]
stopifnot(identical(row.names(ATLAS_GENE_train), ATLAS_GENE_METADATA$gene_id) || identical(row.names(ATLAS_GENE_train), ATLAS_GENE_METADATA$gene_name))

row.names(ATLAS_GENE_train) <- ATLAS_GENE_METADATA$gene_id

common_genes <- intersect(names(ATLAS_GENE_METADATA), row.names(ATLAS_GENE_train))

ATLAS_GENE_train <- ATLAS_GENE_train[common_genes,]
ATLAS_GENE_METADATA <- ATLAS_GENE_METADATA[common_genes,]

ATLAS_GENE_METADATA <- reduce_to_tss(ATLAS_GENE_METADATA)
# stopifnot(identical(colnames(ATLAS_GENE_train), names(mcols(ATLAS_CHIP))))
# sprintf("Atlas contains %s cell types", ncol(mcols(ATLAS_CHIP)))

sample_size = as.numeric(opt$sample_size)

### READ IN CHIP SEQ DATA FROM FOLDER

## Generate GR first
ct = training_ct[1]
chip_ct = import.bw(paste0(opt$chip,sprintf('%s-%s-HG19-%s-RAW.bw', ct, chr, mark)))
# chip_ct = import.wig(gzfile(sprintf('/u/home/f/fujy2038/project-ernst/Project/Correct_data/epimap-chromimpute/converted_tracks/%s_%s_h3k27ac_epimap.wig.wig.gz',chr, ct)))

## Now generate chip seq data for a set of sampled positions
if (!file.exists(file.path(across_dir, proj_name_train, "intermediate-files", "sampled_pos", sprintf("%s_sampled_pos.rds", sample_size)))){
  set.seed(seed)
  sampled_pos = sample(1:nrow(mcols(chip_ct)), sample_size)
  dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "sampled_pos"), recursive=TRUE)
  saveRDS(sampled_pos, file.path(across_dir, proj_name_train, "intermediate-files", "sampled_pos", sprintf("%s_sampled_pos.rds", sample_size)))
}

sampled_pos <- readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "sampled_pos", sprintf("%s_sampled_pos.rds", sample_size)))

chip_ct_chunk = chip_ct[sampled_pos,]
chip_ct_chunk_gr = granges(chip_ct_chunk)

### TODO
for (ct in training_ct){
 chip_ct = import.bw(paste0(opt$chip,sprintf('%s-%s-HG19-%s-RAW.bw', ct, chr, mark)))
 names(mcols(chip_ct)) <- ct
 chip_ct_chunk = chip_ct[sampled_pos,]
 mcols(chip_ct_chunk_gr) = cbind(as.data.frame(mcols(chip_ct_chunk_gr)), as.data.frame(mcols(chip_ct_chunk[,ct])))
}

ATLAS_CHIP <- chip_ct_chunk_gr

### GENERATE TESTING DATA
# X_test <- ATLAS_GENE_test[, holdout, drop = FALSE]

### PARTITION DATA FOR MODEL
set.seed(seed)

X_train <- ATLAS_GENE_train[, training_celltypes, drop = FALSE]
Y_train <- ATLAS_CHIP[, training_celltypes, drop = FALSE]

dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr), recursive=TRUE)
saveRDS(Y_train, file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("ATLAS_CHIP_train_%s.rds", sample_size)))

k_list <- seq(1,k_num_train)

######################################
########## GLOBAL FEATURES ###########
######################################

# ### GENERATE FEATURES FOR TRAINING, VALIDATION AND TESTING CELLTYPES (IN ALL LOCI)
### First generate a distance matrix (distance of each testing cell type to each training cell type)
### Then pick the k nearest nearest cell types for each testing cell type, and average the track of those cell types and use those as features

#### STANDARDIZE GENE EXPRESSION OF TRAINING AND TESTING DATA
X_train_st = scale(X_train)

#### COMPUTE EUCLIDEAN DIST
dist_train = compute_euc_distance(X_train_st, X_train_st)

#### COMPUTE PEARSON CORR
corr_train = compute_corr(X_train_st, X_train_st)

# dir.create(file.path(across_dir, proj_name, "intermediate-files", "avg-tracks-knn-val", sprintf("%s-avg-knn-tracks.rds", valid_ct)))

#### GENERATE GLOBAL KNN FEATURES FROM EUCLIDEAN DIST
if (dist_method == "euc"){
  for (train_ct in training_celltypes){
    first_k_train_ct = names(sort(dist_train[,train_ct]))[1:k_num_train]
    knn_features_all = generate_global_avg(chip=ATLAS_CHIP, k_list=k_list, first_k_train_ct=first_k_train_ct, y_train=Y_train)
    dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr), recursive=TRUE)
    saveRDS(knn_features_all, compress=TRUE, file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("%s-%s-avg-knn-tracks.rds", train_ct, sample_size)))
  }
}

#### GENERATE GLOBAL KNN FEATURES FROM PEARSON CORRELATION
if (dist_method == "corr"){
  for (train_ct in training_celltypes){
      if(!file.exists(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("%s-%s-avg-knn-tracks.rds", train_ct, sample_size)))){
        first_k_train_ct = names(sort(corr_train[,train_ct], decreasing = TRUE, na.last = TRUE))[2:(k_num_train+1)]
        knn_features_all = generate_global_avg(chip=ATLAS_CHIP, k_list=k_list, first_k_train_ct=first_k_train_ct, y_train=Y_train)
        dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr), recursive=TRUE)
        saveRDS(knn_features_all, compress=TRUE, file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("%s-%s-avg-knn-tracks.rds", train_ct, sample_size)))
      }
  }
}

print("Features saved!")

### COMBINE GLOBAL FEATURES FROM DIFFERENT CELLTYPES INTO A LIST
knn_global_train_all = lapply(training_celltypes, function(train_ct){
  knn_global_feature = readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("%s-%s-avg-knn-tracks.rds", train_ct, sample_size)))
  knn_global_feature
})
names(knn_global_train_all) <- training_celltypes

print("Global Features generated!!")

#####################################
########## ADD LOCAL FEATURES ###########
#####################################

#### GENERATE NEAREST GENE IDS FOR EACH LOCI #####
ct <- training_celltypes[1]
if(!file.exists(file.path(across_dir, proj_name_train, "intermediate-files", sprintf("features_%s_%s", chr, sample_size), chr, sprintf("nearest-gene-features.rds", ct)))){
  feature_ct <- feature_eng_train(signal = Y_train, 
                    expr = X_train_st, 
                    expr_metadata = ATLAS_GENE_METADATA, 
                    celltype = ct, 
                    k_nearest = 2, 
                    normalize=FALSE)
  dir.create(file.path(across_dir, proj_name_train, "intermediate-files", sprintf("features_%s_%s", chr, sample_size), chr), recursive=TRUE)
  saveRDS(feature_ct, compress=TRUE, file.path(across_dir, proj_name_train, "intermediate-files", sprintf("features_%s_%s", chr, sample_size), chr, sprintf("nearest-gene-features.rds", ct)))
}

feature_ct <- readRDS(file.path(across_dir, proj_name_train, "intermediate-files", sprintf("features_%s_%s", chr, sample_size), chr, "nearest-gene-features.rds"))
feature_ct

##### GENERATE LIST OF NEAREST CELLTYPES ######

if(!file.exists(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", sprintf("nearest_celltypes_train_l_%s.rds", sample_size)))){
  nearest_celltypes_train_l = list()
  for (gene in row.names(ATLAS_GENE_train)){
    train_df = data.frame(matrix(nrow=k_num_train, ncol=length(training_celltypes)))
    colnames(train_df) = training_celltypes

    if (sd(ATLAS_GENE_train[gene,]) > 0){  ## Could tune this parameter
      for (ct in training_celltypes){
        train_df[,ct] <- names(sort(abs(X_train_st[gene,] - X_train_st[gene,ct])))[2:(k_num_train+1)]
      }
      nearest_celltypes_train_l[[gene]] <- train_df
    }
  }
  dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results"), recursive=TRUE)
  saveRDS(nearest_celltypes_train_l, file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", sprintf("nearest_celltypes_train_l_%s.rds", sample_size)))
}

nearest_celltypes_train_l <- readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", sprintf("nearest_celltypes_train_l_%s.rds", sample_size)))

#### GENERATE LOCAL FEATURES FOR VARIATING GENES ####
Y_train_df = as.data.frame(mcols(Y_train))

### USE PARALLELIZATION TO GENERATE LOCAL FEATURES
#setup parallel backend to use many processors

if(!file.exists(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_local_features_train_ct_list_%s.rds", sample_size)))){
  knn_local_features_train_ct_list <- lapply(training_celltypes, function(train_ct){
    # which_feature_uniq <- 1
    temp_list=list()
    print(train_ct)

    local1_all_features <- lapply(unique(feature_ct$Gene1), function(gene1){
      local1_features <- generate_local_features(gene=gene1,
                          gene_col="Gene1",
                          feature=feature_ct,
                          nearest_celltype_list=nearest_celltypes_train_l,
                          train_ct=train_ct,
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
                          nearest_celltype_list=nearest_celltypes_train_l,
                          train_ct=train_ct,
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
  dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr), recursive=TRUE)
  saveRDS(knn_local_features_train_ct_list, file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_local_features_train_ct_list_%s.rds", sample_size)))
}


### TODO: Simplify the code here.
print("Local Features generated!!")

### COMBINE LOCAL AND GLOBAL FEATURES
k_nearest_dist = 1 ## Hardcoded parameter value. Needs to be changed.

knn_local_features_train_ct_list <- readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_local_features_train_ct_list_%s.rds", sample_size)))

names(knn_local_features_train_ct_list) <- training_celltypes

knn_features_all_ct_train_df <- data.frame(matrix(nrow=0,ncol=2*length(k_list)+2*k_nearest_dist+3))
colnames(knn_features_all_ct_train_df) <- c("Cell_type", "Signal", c(paste0("local1_", k_list), paste0("local2_", k_list)), paste0("Expr", seq(1, k_nearest_dist, 1)), paste0("Dist", seq(1, k_nearest_dist, 1)), "locus")


##################### 22/01/13 TODO: Update the code below to separate proj_name_train and proj_name_test ################

## GENERATE DISTANCE FEATURES FOR REGRESSION TREE
if(!file.exists(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_all_features_train_%s.rds", sample_size)))){
  for (train_ct in training_celltypes){
    # train_ct <- validation_celltypes[[1]]
    # train_ct=training_celltypes[1]
    knn_features_train_ss <- knn_local_features_train_ct_list[[train_ct]]
    knn_features_train_ss_df <- knn_features_train_ss[, c(paste0("local1_", k_list), paste0("local2_", k_list))]
    knn_features_train_ss_df$locus <- knn_features_train_ss$rows
    dist_exp_feature_train <- feature_eng_dist_exp(signal=ATLAS_CHIP, expr=ATLAS_GENE_train, expr_metadata=ATLAS_GENE_METADATA, celltype=train_ct, k_nearest = k_nearest_dist, normalize=FALSE)
    dist_exp_feature_df_train <- as.data.frame(mcols(dist_exp_feature_train)[,-1])
    names(dist_exp_feature_df_train) <- c(paste0("Expr", seq(1, k_nearest_dist, 1)), paste0("Dist", seq(1, k_nearest_dist, 1)))

    knn_features_train_ss_df <- cbind(dist_exp_feature_df_train, knn_features_train_ss_df)
    knn_features_train_ss_df$Cell_type = train_ct
    knn_features_train_ss_df$Signal = unlist(unname(data.frame(mcols(ATLAS_CHIP[,train_ct]))))
    knn_features_train_ss_df <- knn_features_train_ss_df[,c("Cell_type", "Signal", c(paste0("local1_", k_list), paste0("local2_", k_list)), paste0("Expr",seq(1,k_nearest_dist)), paste0("Dist",seq(1,k_nearest_dist)), "locus")]

    knn_features_all_ct_train_df <- rbind(knn_features_all_ct_train_df, knn_features_train_ss_df)
  }
  for(train_ct in training_celltypes){
    knn_features_all = readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "avg-tracks-knn-train", chr, sprintf("%s-%s-avg-knn-tracks.rds", train_ct, sample_size)))
    knn_features_all_df = as.data.frame(mcols(knn_features_all))
    colnames(knn_features_all_df) = paste0("avg_", k_list)
    knn_features_all_ct_train_df[which(knn_features_all_ct_train_df$Cell_type == train_ct),paste0("global_", k_list)] <- knn_features_all_df
  }
  dir.create(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr), recursive=TRUE)
  saveRDS(knn_features_all_ct_train_df, file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_all_features_train_%s.rds", sample_size)))
}

knn_features_all_ct_train_df <- readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_all_features_train_%s.rds", sample_size)))

print("All Features Saved!")

## BUILD REGRESSION TREE 
library(rpart) #for fitting decision trees
library(rpart.plot) #for plotting decision trees

knn_features_all_ct_train_df <- readRDS(file.path(across_dir, proj_name_train, "intermediate-files", "intermediate-results", chr, sprintf("knn_all_features_train_%s.rds", sample_size)))
ptm <- proc.time()

if (!file.exists(file.path(across_dir, proj_name_train, "models", "ensmbl-model", chr, sprintf("reg_model_full_val_%s_%s.rds", train_ct, sample_size)))){
  for(train_ct in training_celltypes){
    knn_features_train_ss_df = knn_features_all_ct_train_df[knn_features_all_ct_train_df$Cell_type == train_ct,]
    local_features <- c(paste0("local1_", k_list), paste0("local2_", k_list))
    global_features <- paste0("global_", k_list)
    dist_features <- paste0("Dist", seq(1,k_nearest_dist))
    exp_features <- paste0("Expr", seq(1,k_nearest_dist))
    full_formula <- as.formula(paste0("Signal ~ ", paste(local_features, collapse = ' + '), "+", paste(global_features, collapse = ' + '), "+", paste(dist_features, collapse = ' + '), "+", paste(exp_features, collapse = ' + ')))
    global_formula <- as.formula(paste0("Signal ~ ", paste(global_features, collapse = ' + '), "+", paste(dist_features, collapse = ' + ')))

    regression_tree_full <- rpart(full_formula, data=knn_features_train_ss_df, method = "anova", control=rpart.control(minbucket=20, cp=0))
    regression_tree_global <- rpart(global_formula, data=knn_features_train_ss_df, method = "anova", control=rpart.control(minbucket=20, cp=0))
    
    best_full <- regression_tree_full$cptable[which.min(regression_tree_full$cptable[,"xerror"]),"CP"]
    best_global <- regression_tree_global$cptable[which.min(regression_tree_global$cptable[,"xerror"]),"CP"]

    pruned_regression_tree_full <- prune(regression_tree_full, cp=best_full)
    pruned_regression_tree_global <- prune(regression_tree_global, cp=best_global)

    #plot the pruned tree
    dir.create(file.path(across_dir, proj_name_train, "models", "ensmbl-model", chr), recursive=TRUE)
    png(filename=file.path(across_dir, proj_name_train, "models", "ensmbl-model", chr, sprintf("Pruned_regression_tree_full_%s_%s.png", sample_size, train_ct)), width=4, height=3, unit="in", res=600)
    prp(pruned_regression_tree_full,
        faclen=0, #use full names for factor labels
        extra=1, #display number of obs. for each terminal node
        roundint=F, #don't round to integers in output
        digits=5) #display 5 decimal places in output
    dev.off()

    png(filename=file.path(across_dir, proj_name_train, "models", "ensmbl-model", chr, sprintf("Pruned_regression_tree_global_%s_%s.png", sample_size, train_ct)), width=4, height=3, unit="in", res=600)
    prp(pruned_regression_tree_global,
        faclen=0, #use full names for factor labels
        extra=1, #display number of obs. for each terminal node
        roundint=F, #don't round to integers in output
        digits=5) #display 5 decimal places in output
    dev.off()
    
    saveRDS(pruned_regression_tree_global, file.path(across_dir, proj_name_train, "models", "ensmbl-model", chr, sprintf("reg_model_global_val_%s_%s.rds", train_ct, sample_size)))
    saveRDS(pruned_regression_tree_full, file.path(across_dir, proj_name_train, "models", "ensmbl-model", chr, sprintf("reg_model_full_val_%s_%s.rds", train_ct, sample_size)))
  }
  print("Regression Tree built!")
  print(proc.time()-ptm)
}


### Remove intermediate files after all regression trees are built
if (!save_intermediate){
  unlink(file.path(across_dir, proj_name_train, "intermediate-files"), recursive = TRUE, force = TRUE)
}



### Remove unused features
# rm(nearest_celltypes_train_l)



