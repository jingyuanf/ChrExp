mc_readRDS <- function(file, threads=parallel::detectCores()) {
  con <- pipe(paste0("pigz -d -c -p",threads," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}



subsetter <- function(g, col) {
  elementMetadata(g)[, col]
}



make_knn_pred = function(k = 1, training, predicting, y, true_y) {
  pred <- FNN::knn.reg(train = training, 
                       test = predicting, 
                       y = y, k = k)$pred
  rmse(predicted = pred, actual = true_y)
}



reduce_to_tss <- function(x) { 
  return( resize(x, 1, fix="start", use.names=TRUE, ignore.strand=FALSE) ) 
}



extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*")) warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  x
}



cluster_loci <- function(x, k = 500, iter = 30, return_stats = TRUE, seed = 42){
  km <- ClusterR::KMeans_arma(as.data.frame(mcols(x)),
                              clusters = k,
                              n_iter = iter,
                              seed_mode = "random_subset",
                              CENTROIDS = NULL, verbose = TRUE, seed = seed)
  pr <- ClusterR::predict_KMeans(as.data.frame(mcols(x)), km) 
  gr <- GenomicRanges::GRanges(seqnames = seqnames(x), ranges = ranges(x), CLUSTER = as.vector(pr))
  
  if(return_stats == TRUE){
    tmp <- cbind.data.frame(mcols(gr), mcols(x))
    tmp %>%
      dplyr::mutate(MEAN = apply(dplyr::select(tmp, -CLUSTER), 1, mean),
                    MEDIAN = apply(dplyr::select(tmp, -CLUSTER), 1, median),
                    MAX = apply(dplyr::select(tmp, -CLUSTER), 1, max)) %>%
      tibble::as_tibble() %>%
      dplyr::select(CLUSTER, MEAN, MEDIAN, MAX) -> mcols(gr)
    
    as.data.frame(mcols(gr)) %>%
      dplyr::group_by(CLUSTER) %>%
      dplyr::summarise(
        N_LOCI = dplyr::n(),
        MEAN = mean(MEAN),
        MEDIAN = mean(MEDIAN),
        MAX = max(MAX)
      ) -> stats
  } else {
    stats <- NULL
  }
  return(list(CLUSTER = gr, STATS = stats))
}



plot_model_results <- function(x, bins, title = ""){
  mse <- mean((x$True - x$Predicted)^2)
  fit <- lm(True ~ Predicted, data = x)
  interval <- predict(fit, interval="confidence")
  x$inside <- ifelse(x$True < interval[,"upr"] & x$True > interval[,"lwr"], "", as.character(x$Cell_type))
  P <- ggplot(x, aes(Predicted, True)) +
    geom_point() +
    geom_smooth(method = "lm", colour = "forestgreen", se=TRUE) + 
    ggrepel::geom_text_repel(aes(label = inside), size=1.5, segment.size = 0.2) +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank()) + geom_abline(slope=1, intercept=0, colour="black", linetype="dotted", size=0.5) +
    labs(title = title,
         subtitle = paste0("Cluster ", unique(bins$CLUSTER), " with ", length(bins), " loci",
                           "\nPearson Corr: ", round(cor(x$Predicted, x$True, method = "pearson"), 5),
                           "\nMSE: ", round(mse, 5)))
  
  return(P)
}



read_signal <- function(x, format='bw', genome='mm9') {
  require(rtracklayer)
  if (format == 'bw') {
    f <- import.bw(x, format='bw')
  } else if (format == 'wig') {
    f <- import.wig(x, format='wig', genome=genome)
  } else {
    message('INVALID FORMAT')
  }
  toRet <- keepSeqlevels(f, paste0('chr', seq(1,19)), pruning.mode="coarse")
  return(toRet)
}



range_norm <- function(x){(x-min(x))/(max(x)-min(x))}



# RMSE <- function(predicted, true){
#   sqrt(mean((predicted - true)^2))
# }


rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}



# normalize_counts <- function(x, sf = 10000){
#   log2((sweep(x, 2, apply(x, 2, function(x) sum(x)), "/") * sf) + 1)
# }



make_scrna_aggregrate <- function(x, labels){
  scrna_aggr <- do.call(cbind.data.frame, lapply(sort(unique(labels)), function(l) {
    cells <- colnames(x[, labels %in% l])
    aggr <- rowSums(x[, colnames(x) %in% cells, drop = FALSE])
  }))
  names(scrna_aggr) <- paste("Cluster", sort(unique(labels)), sep = "_")
  return(scrna_aggr)
}



counts_to_tpm <- function(counts, len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}



write_signal <- function(gr, mcols_name, out, seqinfo=readRDS(file.path("data", "encode-h3k27ac", "mm10-seqinfo.rds"))){
  tmp <- granges(gr)
  mcols(tmp)[["score"]] <- mcols(gr)[[mcols_name]]
  seqinfo(tmp) <- seqinfo
  rtracklayer::export.bw(tmp, con = out)  
  message(sprintf("Wrote %s", out))
}



##########################################################################################################

# The functions below use parallelized versions of gzip, xz, and bzip2 to
# improve compression/decompression performance of RDS serialization in R.
# Each function searches for the appropriate program (based on the required
# compression format) and if found, offloads the compression handling to the
# external program and therefore leaves R free to do the data import/export.
# The two main functions (saveRDS and readRDS) mask R's native read and write
# functions. The functions have been only tested on macOS, but they must work
# on any Linux/Unix.
#
# Requires the following packages: pxz, pbzip2, and pigz.
#
# Run the following line at the command prompt before using the functions.
#
#     brew install pigz pbzip2 pigz
#

library(parallel)

saveRDS.xz <-
  function(object, file, threads = parallel::detectCores()) {
    pxzAvail <- any(grepl("(XZ Utils)", system("pxz -V", intern = TRUE)))
    if (pxzAvail) {
      con <- pipe(paste0("pxz -T", threads, " > ", file), "wb")
      base::saveRDS(object, file = con)
      close(con)
    } else {
      base::saveRDS(object, file = file, compress = "xz")
    }
  }

readRDS.xz <- function(file, threads = parallel::detectCores()) {
  pxzAvail <- any(grepl("(XZ Utils)", system("pxz -V", intern = TRUE)))
  if (pxzAvail) {
    con <- pipe(paste0("pxz -d -k -c -T", threads, " ", file))
    object <- base::readRDS(file = con)
    close(con)
  } else {
    object <- base::readRDS(file)
  }
  return(object)
}

saveRDS.gz <-
  function(object,
           file,
           threads = parallel::detectCores(),
           compression_level = 6) {
    pigzAvail <- any(grepl("pigz", system("pigz -V 2>&1", intern = TRUE)))
    if (pigzAvail) {
      con <-
        pipe(paste0("pigz -c", compression_level, " -p", threads, " > ", file),
             "wb")
      base::saveRDS(object, file = con)
      close(con)
    } else {
      base::saveRDS(object, file = file, compress = "gzip")
    }
  }

readRDS.gz <- function(file, threads = parallel::detectCores()) {
  pigzAvail <- any(grepl("pigz", system("pigz -V 2>&1", intern = TRUE)))
  if (pigzAvail) {
    con <- pipe(paste0("pigz -d -c -p", threads, " ", file))
    object <- base::readRDS(file = con)
    close(con)
  } else {
    object <- base::readRDS(file)
  }
  return(object)
}

saveRDS.bz2 <-
  function(object,
           file,
           threads = parallel::detectCores(),
           compression_level = 6) {
    pbz2Avail <-
      any(grepl("Parallel BZIP2", system("pbzip2 -V 2>&1", intern = TRUE)))
    if (pbz2Avail) {
      con <-
        pipe(paste0("pbzip2 -c", compression_level, " -p", threads, " > ", file),
             "wb")
      base::saveRDS(object, file = con)
      close(con)
    } else {
      base::saveRDS(object, file = file, compress = "bzip2")
    }
  }

readRDS.bz2 <- function(file, threads = parallel::detectCores()) {
  pbz2Avail <-
    any(grepl("Parallel BZIP2", system("pbzip2 -V 2>&1", intern = TRUE)))
  if (pbz2Avail) {
    con <- pipe(paste0("pbzip2 -d -c -p", threads, " ", file))
    object <- base::readRDS(file = con)
    close(con)
  } else {
    object <- base::readRDS(file)
  }
  return(object)
}

readRDS <- function(file, threads = parallel::detectCores()) {
  if (!file.exists(file)) {
    stop(paste0(file, " does not exist!"))
  }
  fileDetails <- system2("file", args = file, stdout = TRUE)
  selector <-
    sapply(c("gzip", "XZ", "BZ"), function (x) {
      grepl(x, fileDetails)
    })
  format <- names(selector)[selector]
  if (length(format) == 0) {
    format <- "not found"
  }
  if (format == "gzip") {
    object <- readRDS.gz(file, threads = threads)
  } else if (format == "XZ") {
    object <- readRDS.xz(file, threads = threads)
  } else if (format == "bzip2") {
    object <- readRDS.bz2(file, threads = threads)
  } else {
    object <- base::readRDS(file)
  }
  return(object)
}

saveRDS <- function(object,
                    file = "",
                    compress = TRUE) {
  if (compress %in% c(TRUE, "gz", "gzip")) {
    saveRDS.gz(object, file)
  } else if (compress %in% c("bzip", "bzip2", "bz", "bz2")) {
    saveRDS.bz2(object, file)
  } else if (compress %in% c("xz", "7zip", "7z")) {
    saveRDS.xz(object, file)
  } else if (compress == FALSE) {
    base::saveRDS(object, file)
  } else {
    stop(paste0(compress, " is not a recognized compression method!"))
  }
}

##########################################################################################################