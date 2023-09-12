### GENERAL FUNCTIONS 

add_locus_dist_info <- function(df_out, df_src){
    if(nrow(df_out) == nrow(df_src)){
       df_out$locus <- df_src$locus
       df_out$Dist1 <- df_src$Dist1
       return(df_out)
    } else {
        stop("The number of rows of two data frames are not the same!")
    }
}

make_long <- function(df, ct, value){
  df_long <- pivot_longer(df, cols=all_of(ct), names_to="Cell_type", values_to=value)
  df_long <- as.data.frame(df_long)
  return(df_long)
}

read_wig <- function(x, format='wig', genome='mm9') {
  # keepStandardChromosomes()
  suppressMessages(library(rtracklayer))  
  merged_wig <- import.wig(x, format=format, genome=genome)
  merged_wig <- keepSeqlevels(merged_wig, paste0('chr', c(seq(1,19), 'X', 'Y')), pruning.mode="coarse")
  return(merged_wig)
}

extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*")){ warning("'*' ranges were treated as '+'") }
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  # trim(x)
  x
}

reduce_to_tss <- function(x) {
  return( resize(x, 1, fix="start", use.names=TRUE, ignore.strand=FALSE) ) 
}

mround <- function(x,base){ 
  base*round(x/base) 
} 

### FUNCTIONS FOR READING IN DATA

read_across <- function(across_dir, proj_name_across, chr){
    across_full_pred_gr <- try(readRDS(file.path(across_dir, proj_name_across, "final-pred", chr, "across_gr.rds")), silent = TRUE)
    if (class(across_full_pred_gr) == "try-error"){
        stop("Across prediction read error!")
    }
    across_full_pred_gr <- subset_to_chr(across_full_pred_gr, chr)
    return(across_full_pred_gr)
}

read_across_old <- function(across_dir, proj_name_across_old, chr){
    across_full_pred_gr_old <- try(readRDS(file.path(across_dir, proj_name_across_old, "intermediate-files", "ensmbl-model", chr, "across_gr_log2.rds")), silent = TRUE)
    if (class(across_full_pred_gr_old) == "try-error"){
        stop("Old across prediction read error!")
    }
    across_full_pred_gr_old <- subset_to_chr(across_full_pred_gr_old, chr)
    return(across_full_pred_gr_old)
}

read_across_global <- function(across_dir, proj_name_across_global, chr){
    across_global_pred <- try(readRDS(file.path(across_dir, proj_name_across_global, "intermediate-files", "ensmbl-model", chr, "testing_pred_df_reg_global.rds")), silent = TRUE)
    if (class(across_global_pred) == "try-error"){
        stop("Across prediction (only global) read error!")
    }
    across_global_pred <- subset_to_chr(across_global_pred, chr)
    return(across_global_pred)
}

read_within <- function(within_dir, proj_name_within, chr){
    within_pred <- try(readRDS(file.path(within_dir, proj_name_within, "final_pred", "2", "with_distance", num_train_ct, "bw_files", sprintf("within_gr_log2_%s.rds", chr))), silent = TRUE)
    if (class(within_pred) == "try-error"){
        stop("Within prediction read error!")
    }  
    within_pred <- subset_to_chr(within_pred, chr)
    return(within_pred)
}

read_integration <- function(int_dir, proj_name_int, chr, thrd){
    int_pred = try(readRDS(sprintf("./bss-predictions/integration/%s/integration-acr-within-thrd-%s-new.rds", proj_name_int, thrd)))
    if (class(int_pred) == "try-error"){
        stop("Integration prediction read error!")
    }  
    int_pred <- subset_to_chr(int_pred, chr)
    return(int_pred)
}

read_chromimpute <- function(chrom_dir, proj_name_chrom, chr){
    chrom_pred = try(readRDS(file.path(chrom_dir, proj_name_chrom, sprintf("chrom_gr_%s.rds", chr))))
    if (class(chrom_pred) == "try-error"){
        stop("ChromImpute prediction read error!")
    }  
    chrom_pred <- subset_to_chr(chrom_pred, chr)
    return(chrom_pred)
}

### FUNCTIONS FOR METRIC COMPUTATION

### PART 1: FUNCTIONS FOR GENERATING TOTAL CROSS-LOCI CORRELATIONS AND RELATED PLOTTING
compute_total_cross_loci_mcols <- function(pred_df, truth_df, mark, corr_method, pred_method){
  if (pred_method != "Mean_baseline"){
    celltypes = intersect(names(pred_df), names(truth_df))
    corr_list <- lapply(celltypes, function(ct){
      corr <- cor(pred_df[,ct], truth_df[,ct], method=corr_method)
      corr
  })
  } else {
    celltypes = names(truth_df)
    corr_list <- lapply(celltypes, function(ct){
      corr <- cor(pred_df, truth_df[,ct], method=corr_method)
      corr
    })
  }

  corr_list <- unlist(corr_list)
  corr_df <- data.frame(corr = corr_list, mark = mark, method = pred_method, corr_method = corr_method)
  return(corr_df)
}

# wilcox.test(corr_df)

compute_ct_specific_cross_loci_mcols <- function(pred_df, truth_df, train_df, mark, corr_method, pred_method, proj_name, thrd){
  num_train_ct = ncol(train_df)
  num_other_ct_have_peak <- seq(0, num_train_ct-1)

  if(!file.exists(sprintf("./evaluations/%s/rows_to_test_l_%s_%s.rds", proj_name, thrd, mark))){
    dir.create(sprintf("./evaluations/%s/", proj_name), recursive=TRUE)
    rows_to_test_l <- generate_rows_to_test_l(train_df, truth_df, num_train_ct, thrd)
    saveRDS(rows_to_test_l, sprintf("./evaluations/%s/rows_to_test_l_%s_%s.rds", proj_name, thrd, mark))
  }
  rows_to_test_l <- readRDS(sprintf("./evaluations/%s/rows_to_test_l_%s_%s.rds", proj_name, thrd, mark))
  # rows_to_test_l <- readRDS(sprintf("./evaluations/%s/rows_to_test_l_%s.rds", proj_name, thrd))

  corr_df <- data.frame(matrix(nrow=0, ncol=6))
  names(corr_df) <- c("Correlation", "Celltype", "Method", "Corr_method", "Mark", "Num_ct_hold_peak")

  if (pred_method != "Mean_baseline"){
    celltypes = intersect(names(pred_df), names(truth_df))

    for (num_other_ct in num_other_ct_have_peak){
      rows_to_test <- rows_to_test_l[[toString(num_other_ct)]]
      if(length(rows_to_test) > 50){
        for (ct in celltypes){
          corr <- cor(pred_df[rows_to_test,ct], truth_df[rows_to_test,ct], method = corr_method)
          corr_add <- data.frame(Correlation = corr, Celltype = ct, Method = pred_method, Corr_method = corr_method, Mark = mark, Num_ct_hold_peak = num_other_ct)
          corr_df <- rbind(corr_df, corr_add)
        }
      }
    }
  } else {
    celltypes = names(truth_df)

    for (num_other_ct in num_other_ct_have_peak){
      rows_to_test <- rows_to_test_l[[toString(num_other_ct)]]
      if(length(rows_to_test) > 50){
        for (ct in celltypes){
          corr <- cor(pred_df[rows_to_test], truth_df[rows_to_test,ct], method = corr_method)
          corr_add <- data.frame(Correlation = corr, Celltype = ct, Method = pred_method, Corr_method = corr_method, Mark = mark, Num_ct_hold_peak = num_other_ct)
          corr_df <- rbind(corr_df, corr_add)
        }
      }
    }    
  }
  return(corr_df)
}

generate_rows_to_test_l <- function(train_df, test_df, num_train_ct, thrd){
  rows_to_test_l <- list()
  num_other_ct_have_peak <- seq(0,num_train_ct-1)
  training_celltypes <- colnames(train_df)
  testing_celltypes <- colnames(test_df)
  for (num_other_ct in num_other_ct_have_peak){
    if (num_other_ct < max(num_other_ct_have_peak)){
      rows_to_test <- which(rowSums(train_df[,training_celltypes] >= thrd)  == num_other_ct & rowSums(test_df[,testing_celltypes] >= thrd)  > 0)
    } else {
      rows_to_test <- which(rowSums(train_df[,training_celltypes] >= thrd)  >= num_other_ct & rowSums(test_df[,testing_celltypes] >= thrd)  > 0)
    }
    rows_to_test_l[[toString(num_other_ct)]] <- rows_to_test
  }
  return(rows_to_test_l)
}

############### OTHER EVALUATION FUNCTIONS ####################
compute_cross_ct_corr_dist <- function(rows_to_test, df_long, pred_list, method, dist_categories){
    across_ct_corr_df_pred = data.frame(matrix(nrow=0, ncol=length(pred_list)+1))
    names(across_ct_corr_df_pred) = c(pred_list, "dist")
    for (num_dist in 1:length(dist_categories)){
        if (num_dist == 1){
            k_feature_ss <- df_long[df_long$Dist1 < dist_categories[num_dist], ]
        } else {
            k_feature_ss <- df_long[df_long$Dist1 < dist_categories[num_dist] & df_long$Dist1 >= dist_categories[num_dist-1], ]
        }
        message(num_dist)
        rows_to_eval = intersect(k_feature_ss$locus, rows_to_test)
        for (locus in rows_to_eval){
            locus_ss = which(k_feature_ss$locus == locus)
            k_feature_locus_ss = k_feature_ss[locus_ss,]
            temp_df <- data.frame(matrix(nrow=1, ncol=5))
            names(temp_df) = c(pred_list, "dist")
            for(pred in names(temp_df)[1:length(pred_list)]){
                temp_df[,pred] = cor(k_feature_locus_ss[,pred], k_feature_locus_ss[,"Signal"], method=method)
            }
            temp_df$dist = dist_categories[num_dist]
            across_ct_corr_df_pred <- rbind(across_ct_corr_df_pred, temp_df)
        }
    }
    return(across_ct_corr_df_pred)
}

compute_cross_loci_corr_dist <- function(df_long, pred_list, method, dist_categories){
    across_loci_corr_df_pred = data.frame(matrix(nrow=0, ncol=length(pred_list)+1))
    names(across_loci_corr_df_pred) = c(pred_list, "dist")
    for (num_dist in 1:length(dist_categories)){
        if (num_dist == 1){
            k_feature_ss <- df_long[df_long$Dist1 < dist_categories[num_dist], ]
        } else {
            k_feature_ss <- df_long[df_long$Dist1 < dist_categories[num_dist] & df_long$Dist1 >= dist_categories[num_dist-1], ]
        }
        message(num_dist)
        for (celltype in unique(k_feature_ss$Cell_type)){
            celltype_ss = which(k_feature_ss$Cell_type == celltype)
            k_feature_celltype_ss = k_feature_ss[celltype_ss,]
            temp_df <- data.frame(matrix(nrow=1, ncol=5))
            names(temp_df) = c(pred_list, "dist")
            for(pred in names(temp_df)[1:length(pred_list)]){
                temp_df[,pred] = cor(k_feature_celltype_ss[,pred], k_feature_celltype_ss[,"Signal"], method=method)
            }
            temp_df$dist = dist_categories[num_dist]
            across_loci_corr_df_pred <- rbind(across_loci_corr_df_pred, temp_df)
        }
    }
    return(across_loci_corr_df_pred)
}

compute_cross_loci_corr_strat_ct <- function(df_long, num_other_ct_have_peak, rows_to_test_l,corr_method="pearson", pred_list, true_mat, train_mat){
  corr_df <- data.frame(matrix(nrow=0, ncol=5))
  corr_method="pearson"
  names(corr_df) = c("num_train_ct_w_signal", corr_method, "method", "value", "celltype")
  train_tri = rowSums(train_mat)/ncol(train_mat)
  tri = rowSums(true_mat)/ncol(true_mat)

  for (num_other_ct in num_other_ct_have_peak){
    rows_to_test <- rows_to_test_l[[toString(num_other_ct)]]
    df_ss_locus = df_long[df_long$locus %in% rows_to_test, ]
    train_tri_ss = train_tri[rows_to_test]
    tri_ss = tri[rows_to_test]
    for (ct in testing_celltypes){
      df_ss_locus_ct = df_ss_locus[df_ss_locus$Cell_type == ct,]
      for (pred_m in pred_list){
          corr_val = cor(df_ss_locus_ct[,pred_m], true_mat[rows_to_test,ct], method = corr_method)
          corr_df = rbind(corr_df, data.frame(num_train_ct_w_signal=num_other_ct, corr_method=corr_method, method=pred_m, value=corr_val, celltype=ct))
      }
      train_tri_val = cor(train_tri_ss, true_mat[rows_to_test,ct], method = corr_method)
      corr_df = rbind(corr_df, data.frame(num_train_ct_w_signal=num_other_ct, corr_method=corr_method, method="train_average", value=train_tri_val, celltype=ct))

      tri_val = cor(tri_ss, true_mat[rows_to_test,ct], method = corr_method)
      corr_df = rbind(corr_df, data.frame(num_train_ct_w_signal=num_other_ct, corr_method=corr_method, method="truth_average", value=tri_val, celltype=ct))
    }
  }
  return(corr_df)
}

compute_cross_loci_corr_strat_state <- function(df_long, corr_method="pearson", pred_list, true_mat, train_mat){
  corr_df <- data.frame(matrix(nrow=0, ncol=5))
  corr_method="pearson"
  names(corr_df) = c("state", corr_method, "method", "value", "celltype")
  train_tri = rowSums(train_mat)/ncol(train_mat)
  tri = rowSums(true_mat)/ncol(true_mat)
  state_l <- sort(unique(df_long$state))

  for (state in state_l){
    for (ct in testing_celltypes){
      df_ss_ct = df_long[df_long$Cell_type == ct,]
      df_ss_ct = df_ss_ct[order(df_ss_ct$locus),]
      rows_to_test = which(df_ss_ct$state == state)
      train_tri_ss = train_tri[rows_to_test]
      tri_ss = tri[rows_to_test]
      df_ss_locus = df_ss_ct[rows_to_test, ]

      for (pred_m in pred_list){
          corr_val = cor(df_ss_locus[,pred_m], true_mat[rows_to_test,ct], method = corr_method)
          corr_df = rbind(corr_df, data.frame(state=state, corr_method=corr_method, method=pred_m, value=corr_val, celltype=ct))
      }
      train_tri_val = cor(train_tri_ss, true_mat[rows_to_test,ct], method = corr_method)
      corr_df = rbind(corr_df, data.frame(state=state, corr_method=corr_method, method="train_average", value=train_tri_val, celltype=ct))

      tri_val = cor(tri_ss, true_mat[rows_to_test,ct], method = corr_method)
      corr_df = rbind(corr_df, data.frame(state=state, corr_method=corr_method, method="truth_average", value=tri_val, celltype=ct))
    }
  }
  return(corr_df)
}

compute_cross_ct_corr_strat_state <- function(df_long, corr_method="pearson", pred_list, true_mat, num_smp=100){
  corr_df <- data.frame(matrix(nrow=0, ncol=5))
  corr_method="pearson"
  names(corr_df) = c("state", corr_method, "method", "value", "locus")
  state_l <- sort(unique(df_long$state))

  for (state in state_l){
    df_ss_state = df_long[df_long$state %in% state, ]
    locus_state = unique(df_ss_state$locus)
    set.seed(10)
    if (length(locus_state) > num_smp){
        locus_state_smp <- sample(locus_state, num_smp)
    } else {
        locus_state_smp <- locus_state
    }    
    
    for (locus in locus_state_smp){
      row_ss = which(df_ss_state$locus == locus)
      for (pred_m in pred_list){
          ct = df_ss_state[row_ss, "Cell_type"]
          corr_val = cor(df_ss_state[row_ss, pred_m], unlist(unname(true_mat[locus,ct])), method = corr_method)
          corr_df = rbind(corr_df, data.frame(state=state, corr_method=corr_method, method=pred_m, value=corr_val, locus=locus))
      }
    }
  }
  return(corr_df)
}

compute_cross_ct_corr_strat_ct <- function(df_long, num_other_ct_have_peak, rows_to_test_l,corr_method="pearson", pred_list, true_mat, num_smp=100){
  corr_df <- data.frame(matrix(nrow=0, ncol=5))
  corr_method="pearson"
  names(corr_df) = c("num_train_ct_w_signal", corr_method, "method", "value", "locus")
  for (num_other_ct in num_other_ct_have_peak){
    rows_to_test <- rows_to_test_l[[toString(num_other_ct)]]
    df_ss_locus = df_long[df_long$locus %in% rows_to_test, ]
    message(num_other_ct)
    set.seed(10)
    if (length(rows_to_test) > num_smp){
        rows_to_test_smp <- sample(rows_to_test, num_smp)
    } else {
        rows_to_test_smp <- rows_to_test
    }    
    
    for (row in rows_to_test_smp){
      row_ss = which(df_long$locus == row)
      for (pred_m in pred_list){
          corr_val = cor(df_long[row_ss, pred_m], unlist(unname(true_mat[row,])), method = corr_method)
          corr_df = rbind(corr_df, data.frame(num_train_ct_w_signal=num_other_ct, corr_method=corr_method, method=pred_m, value=corr_val, locus=row))
      }
    }
  }
  return(corr_df)
}

compute_cross_loci_corr_general <- function(df_long, corr_method="pearson", pred_list){
  corr_df <- data.frame(matrix(nrow=0, ncol=4))
  corr_method="pearson"
  names(corr_df) = c(corr_method, "method", "value", "Cell_type")
  celltypes = unique(df_long$Cell_type)
  for (ct in celltypes){
      row_ss = which(df_long$Cell_type == ct)
      for (pred_m in pred_list){
          corr_val = cor(df_long[row_ss, pred_m], df_long[row_ss,"Signal"], method = corr_method)
          corr_df = rbind(corr_df, data.frame(corr_method=corr_method, method=pred_m, value=corr_val, Cell_type=ct))
      }
    }
  return(corr_df)
}

compute_cross_ct_corr_strat_ct_acr_chr <- function(df_long, num_other_ct_have_peak, rows_to_test_l,corr_method="pearson", pred_list,num_smp=100){
  corr_df <- data.frame(matrix(nrow=0, ncol=5))
  corr_method="pearson"
  names(corr_df) = c("num_train_ct_w_signal", corr_method, "method", "value", "locus")
  for (num_other_ct in num_other_ct_have_peak){
    rows_to_test <- rows_to_test_l[[toString(num_other_ct)]]
    df_ss_locus = df_long[df_long$locus %in% rows_to_test, ]
    message(num_other_ct)
    set.seed(10)
    if (length(rows_to_test) > num_smp){
        rows_to_test_smp <- sample(rows_to_test, num_smp)
    } else {
        rows_to_test_smp <- rows_to_test
    }
    
    for (row in rows_to_test_smp){
      row_ss = which(df_long$locus == row)
      corr_val = cor(df_long[row_ss, "across_pred"], df_long[row_ss, "chrom_pred"], method = corr_method)
      corr_df = rbind(corr_df, data.frame(num_train_ct_w_signal=num_other_ct, corr_method=corr_method, value=corr_val, locus=row))
    }
  }
  return(corr_df)
}

compute_cross_ct_corr_strat_ct_acr_acr_old <- function(df_long, num_other_ct_have_peak, rows_to_test_l,corr_method="pearson", pred_list,num_smp=100){
  corr_df <- data.frame(matrix(nrow=0, ncol=5))
  corr_method="pearson"
  names(corr_df) = c("num_train_ct_w_signal", corr_method, "method", "value", "locus")
  for (num_other_ct in num_other_ct_have_peak){
    rows_to_test <- rows_to_test_l[[toString(num_other_ct)]]
    df_ss_locus = df_long[df_long$locus %in% rows_to_test, ]
    message(num_other_ct)
    set.seed(10)
    if (length(rows_to_test) > num_smp){
        rows_to_test_smp <- sample(rows_to_test, num_smp)
    } else {
        rows_to_test_smp <- rows_to_test
    }
    
    for (row in rows_to_test_smp){
      row_ss = which(df_long$locus == row)
      corr_val = cor(df_long[row_ss, "across_pred"], df_long[row_ss, "across_pred_old"], method = corr_method)
      corr_df = rbind(corr_df, data.frame(num_train_ct_w_signal=num_other_ct, corr_method=corr_method, value=corr_val, locus=row))
    }
  }
  return(corr_df)
}

corr_acr_loci_general_two_cols <- function(df_long, col1, col2, method="pearson"){
    celltypes <- unique(df_long$Cell_type)
    corr_val_v <- lapply(celltypes, function(ct){
        col1_ct <- df_long[df_long$Cell_type == ct, col1]
        col2_ct <- df_long[df_long$Cell_type == ct, col2]
        corr_val <- cor(col1_ct, col2_ct, method=method)
        return(corr_val)
    })
    return(unlist(corr_val_v))
}

