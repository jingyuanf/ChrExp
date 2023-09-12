feature_eng_dist_exp <- function(signal, expr, expr_metadata, celltype, k_nearest = 1, normalize=FALSE) {
  message(sprintf("%s\tknn=%s", celltype, k_nearest))
  ct_gene <- expr[, celltype, drop = FALSE]
  ct_signal <- signal[, celltype, drop = FALSE]
  mcols(expr_metadata)$Expr <- expr[, celltype]
  knn <- findKNN(ct_signal, expr_metadata[seqnames(expr_metadata) %in% seqlevels(ct_signal),], k = k_nearest, ignore.strand = TRUE)
  knn_expr <- extractList(unname(expr_metadata$Expr), knn)

  # expr_metadata_tss = reduce_to_tss(expr_metadata)
  knn_dists <- abs(extractList(start(expr_metadata), knn) - start(ct_signal))
  
  d <- as.matrix(knn_dists)
  e <- as.matrix(knn_expr)

  if(normalize == TRUE){
    d <- (1/d)/rowSums(1/d)
    e <- e/rowSums(e)
    e <- t(scale(t(e), center = TRUE, scale = TRUE))
  }

  knn_feats <- cbind.data.frame(e,d)
  names(knn_feats) <- c(paste0("Expr", seq(1, k_nearest, 1)), paste0("Dist", seq(1, k_nearest, 1)))
  mcols(ct_signal) <- cbind(mcols(ct_signal), knn_feats)
  return(ct_signal)
}

feature_eng_dist_exp_test <- function(signal, expr, expr_metadata, celltype, k_nearest = 1, normalize=FALSE) {
  message(sprintf("%s\tknn=%s", celltype, k_nearest))
  ct_gene <- expr[, celltype, drop = FALSE]
  ct_signal <- signal[, 1, drop = FALSE]
  mcols(expr_metadata)$Expr <- expr[, celltype]
  knn <- findKNN(ct_signal, expr_metadata[seqnames(expr_metadata) %in% seqlevels(ct_signal),], k = k_nearest, ignore.strand = TRUE)
  knn_expr <- extractList(unname(expr_metadata$Expr), knn)

  # expr_metadata_tss = reduce_to_tss(expr_metadata)
  knn_dists <- abs(extractList(start(expr_metadata), knn) - start(ct_signal))
  
  d <- as.matrix(knn_dists)
  e <- as.matrix(knn_expr)

  if(normalize == TRUE){
    d <- (1/d)/rowSums(1/d)
    e <- e/rowSums(e)
    e <- t(scale(t(e), center = TRUE, scale = TRUE))
  }

  knn_feats <- cbind.data.frame(e,d)
  names(knn_feats) <- c(paste0("Expr", seq(1, k_nearest, 1)), paste0("Dist", seq(1, k_nearest, 1)))
  mcols(ct_signal) <- cbind(mcols(ct_signal), knn_feats)
  return(ct_signal)
}

feature_eng_train <- function(signal, expr, expr_metadata, celltype, k_nearest = 1, normalize=FALSE) {
  # message(sprintf("%s\tknn=%s", celltype, k_nearest))
  ct_gene <- expr[, celltype, drop = FALSE]
  ct_signal <- signal[, celltype, drop = FALSE]
  mcols(expr_metadata)$Expr <- expr[, celltype]
  knn <- findKNN(ct_signal, expr_metadata[seqnames(expr_metadata) %in% seqlevels(ct_signal),], k = k_nearest, ignore.strand = TRUE)
  knn_genes <- extractList(expr_metadata$gene_id, knn)

  g <- as.data.frame(as.matrix(knn_genes))

  knn_feats <- g
  colnames(knn_feats) <- paste0("Gene", seq(1, k_nearest, 1))

  mcols(ct_signal) <- cbind(mcols(ct_signal), knn_feats)
  return(ct_signal)
}

feature_eng_test <- function(signal, expr, expr_metadata, celltype, k_nearest = 1, normalize) {
  message(sprintf("%s\tknn=%s", celltype, k_nearest))
  ct_gene <- expr[, celltype, drop = FALSE]
  ct_signal <- granges(signal)
  mcols(expr_metadata)$Expr <- expr[, celltype]
  knn <- findKNN(ct_signal, expr_metadata[seqnames(expr_metadata) %in% seqlevels(ct_signal),], k = k_nearest, ignore.strand = TRUE)
  knn_genes <- extractList(expr_metadata$gene_name, knn)
   g <- as.data.frame(as.matrix(knn_genes))
  
  knn_feats <- g
  colnames(knn_feats) <- paste0("Gene", seq(1, k_nearest, 1))
  mcols(ct_signal) <- cbind(mcols(ct_signal), knn_feats)
  return(ct_signal)
}

compute_euc_distance <- function(train, test){
  dist_mat = matrix(nrow=ncol(train), ncol=ncol(test))
  colnames(dist_mat) <- colnames(test)
  row.names(dist_mat) <- colnames(train)

  for (train_ct in colnames(train)){
    for (test_ct in colnames(test)){
      dist = sqrt(sum((train[,train_ct] - test[,test_ct])^2))
      dist_mat[train_ct, test_ct] = dist
    }
  }
  return(dist_mat)
}

compute_corr <- function(train, test){
  dist_mat = matrix(nrow=ncol(train), ncol=ncol(test))
  colnames(dist_mat) <- colnames(test)
  row.names(dist_mat) <- colnames(train)

  for (train_ct in colnames(train)){
    for (test_ct in colnames(test)){
      dist = cor(train[,train_ct], test[,test_ct])
      dist_mat[train_ct, test_ct] = dist
    }
  }
  return(dist_mat)
}

generate_global_avg <- function(chip, k_list, first_k_train_ct, y_train){
  knn_features_all = granges(chip)
  mcols(knn_features_all) <- data.frame(matrix(nrow=length(knn_features_all), ncol=length(k_list)))
  names(mcols(knn_features_all)) <- k_list
  for (ind in 1:length(k_list)){
    k_gene = k_list[ind]
    mcols(knn_features_all)[,ind] <- rowSums(as.data.frame(mcols(Y_train[,first_k_train_ct[1:k_gene]])))/k_gene
  }
  return(knn_features_all)
}

# generate_global_avg <- function(chip, k_list, first_k_train_ct, y_train){
#   knn_features_all = granges(chip)
#   mcols(knn_features_all) <- data.frame(matrix(nrow=length(knn_features_all), ncol=length(k_list)))
#   names(mcols(knn_features_all)) <- k_list
#   for (ind in 1:length(k_list)){
#     k_gene = k_list[ind]
#     mcols(knn_features_all)[,ind] <- rowSums(as.data.frame(mcols(Y_train[,first_k_train_ct[1:k_gene]])))/k_gene
#   }
#   return(knn_features_all)
# }

generate_local_features <- function(gene, gene_col, feature, nearest_celltype_list, train_ct, k_list, y_train_df, new_col_name){
  rows_gene = which(mcols(feature)[,gene_col] == gene)
  nearest_train_ct = nearest_celltype_list[[gene]][,train_ct]
  local_features = data.frame(matrix(nrow=length(rows_gene), ncol=length(k_list)))
  for (ind in 1:length(k_list)){
    k_gene = k_list[ind]
    local_features[,ind] <- rowSums(as.data.frame(y_train_df[rows_gene,nearest_train_ct[1:k_gene]]))/k_gene
  }
  names(local_features) <- paste0(new_col_name, k_list)
  local_features$rows <- rows_gene
  local_features[,gene_col] <- gene
  return(local_features)
}

perform_rpart <- function(formula, data, minbucket=1, plot_file_name, save_file_name, save_plot=TRUE, save_model=TRUE){
  regression_tree <- rpart(formula, data=data, method = "anova", control=rpart.control(minbucket=minbucket, cp=0))
  best_reg <- regression_tree$cptable[which.min(regression_tree$cptable[,"xerror"]),"CP"]
  pruned_regression_tree <- prune(regression_tree, cp=best_reg)

  if(save_plot){
    png(filename=plot_file_name, width=4, height=3, unit="in", res=600)
    prp(pruned_regression_tree,
        faclen=0, #use full names for factor labels
        extra=1, #display number of obs. for each terminal node
        roundint=F, #don't round to integers in output
        digits=5) #display 5 decimal places in output
    dev.off()
  }

  if(save_model){
    saveRDS(pruned_regression_tree, compress=TRUE, save_file_name)
  }
}


perform_glmnet <- function(predictor, target, features, family = "gaussian", plot_file_name, save_file_name, save_plot=TRUE, save_model=TRUE){
  glmnet_pred <- cv.glmnet(x=predictor, y=target, family=family)
  
  ## Get coefficients
  glmnet_coefficients <- coef(glmnet_pred)
  glmnet_coefficients <- as.numeric(glmnet_coefficients)
  names(glmnet_coefficients) <- c('(Intercept)', features)
  glmnet_ss_coefficients <- glmnet_coefficients[which(glmnet_coefficients != 0)]

  lambda = glmnet_pred$lambda.1se

  if(save_plot){
    png(filename=plot_file_name, width=8, height=6, unit="in", res=600)
    n <- length(glmnet_ss_coefficients)
    barplot(glmnet_ss_coefficients, main = sprintf("Lasso Regression Coefficients for Full Model, Lambda = %s", lambda), xlab = "Coefficient", ylab = "Estimate", xaxt = "n")
    axis(1, at = 1:n, labels = names(glmnet_ss_coefficients), las = 2, cex.axis = 0.5)
    dev.off()
  }

  if(save_model){
    saveRDS(glmnet_ss_coefficients, compress=TRUE, save_file_name)
  }
}


make_testing_pred_df <- function(nrow_test, testing_celltypes, data, colname_in_data){
  testing_pred_df <- data.frame(matrix(nrow=nrow_test, ncol=length(testing_celltypes)))
  colnames(testing_pred_df) <- testing_celltypes

  for(test_ct in testing_celltypes){
    print(test_ct)
    pred_signal = data[data$Cell_type == test_ct,colname_in_data]
    testing_pred_df[,test_ct] <- pred_signal
  }
  return(testing_pred_df)
}

pred_from_coef <- function(coefficients, df){
  coef_df <- data.frame(names = names(coefficients), coef = unname(coefficients))
  if (length(coefficients) > 1){
    features <- df[,names(coefficients)[2:length(names(coefficients))]]
    coef_df$values <- c(1, unname(unlist(features)))
    return(sum(coef_df$coef * coef_df$values))
  } else {
    return(unname(coefficients))
  }
}