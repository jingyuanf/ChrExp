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
  ## TODO: Check other clustering methods such as hierarchical clustering
  
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


findKNN <- function(query, subject, k=5L, ignore.overlaps = FALSE,
                    ignore.strand = FALSE)
{
    seqlevels(subject) <- seqlevels(query)
    
    starts <- with(subject, GRanges(seqnames, IRanges(start, width=1L), strand))
    ends <- with(subject, GRanges(seqnames, IRanges(end, width=1L), strand))

    if (ignore.strand) {
        starts <- unstrand(starts)
        ends <- unstrand(ends)
    }
    
    start_ord <- order(starts)
    end_ord <- order(ends)

    starts <- starts[start_ord]
    ends <- ends[end_ord]

    phits <- precede(query, starts, ignore.strand=ignore.strand)
    fhits <- follow(query, ends, ignore.strand=ignore.strand)

    if (!ignore.strand) {
        exchange <- decode(strand(query) == "-")
        tmp <- phits[exchange]
        phits[exchange] <- fhits[exchange]
        fhits[exchange] <- tmp
    }

    findPart <- function(x, w) {
        S4Vectors:::findIntervalAndStartFromWidth(x, w)[["interval"]]
    }
    
    if (!ignore.strand) {
        b <- width(disjoin(c(ranges(seqnames(starts)), ranges(strand(starts)))))
    } else {
        b <- runLength(seqnames(starts))
    }

which(expr_metadata$gene_name == "NEK2")
which(expr_metadata$gene_name == "LPGAT1")


    # b <- runLength(seqnames(starts))
    
    seqends <- end(seqnames(starts))[findPart(phits, b)]
    phits[is.na(phits)] <- 1L
    seqends[is.na(seqends)] <- 0L
    pwindows <- restrict(IRanges(phits, width = k), end=seqends, keep.all.ranges=TRUE)

    seqstarts <- start(seqnames(ends))[findPart(fhits, b)]
    seqstarts[is.na(seqstarts)] <- 1L
    fhits[is.na(fhits)] <- 0L
    fwindows <- restrict(IRanges(end=fhits, width = k), seqstarts, keep.all.ranges=TRUE)
    
    dist <- pc(extractList(start(starts), pwindows) - end(query), ### pwindows=5
               end(query) - extractList(end(ends), fwindows)) ### Try to write another function that disregards preceding/follows or somethingelse, just look at absolute distance regarding TSS. 

    ans <- pc(extractList(start_ord, pwindows), extractList(end_ord, fwindows))

    if (!ignore.overlaps) {
        hits <- findOverlaps(query, subject, ignore.strand=ignore.strand)
        hitsList <- as(hits, "List")
        dist <- pc(dist, relist(rep(0L, length(hits)), hitsList))
        ans <- pc(ans, hitsList)
    }

    ans[heads(order(dist), k)]
}


findKNN <- function(query, subject, k=5L, ignore.overlaps = FALSE,
                    ignore.strand = FALSE)
{
    seqlevels(subject) <- seqlevels(query)
    
    starts <- with(subject, GRanges(seqnames, IRanges(start, width=1L), strand))
    ends <- with(subject, GRanges(seqnames, IRanges(end, width=1L), strand))

    if (ignore.strand) {
        starts <- unstrand(starts)
        ends <- unstrand(ends)
    }
    
    start_ord <- order(starts)
    end_ord <- order(ends)

    starts <- starts[start_ord]
    ends <- ends[end_ord]

    phits <- precede(query, starts, ignore.strand=ignore.strand)
    fhits <- follow(query, ends, ignore.strand=ignore.strand)

    if (!ignore.strand) {
        exchange <- decode(strand(query) == "-")
        tmp <- phits[exchange]
        phits[exchange] <- fhits[exchange]
        fhits[exchange] <- tmp
    }

    findPart <- function(x, w) {
        S4Vectors:::findIntervalAndStartFromWidth(x, w)[["interval"]]
    }
    
    if (!ignore.strand) {
        b <- width(disjoin(c(ranges(seqnames(starts)), ranges(strand(starts)))))
    } else {
        b <- runLength(seqnames(starts))
    }

# which(expr_metadata$gene_name == "NEK2")
# which(expr_metadata$gene_name == "LPGAT1")


    # b <- runLength(seqnames(starts))
    
    seqends <- end(seqnames(starts))[findPart(phits, b)]
    phits[is.na(phits)] <- 1L
    seqends[is.na(seqends)] <- 0L
    pwindows <- restrict(IRanges(phits, width = k), end=seqends, keep.all.ranges=TRUE)

    seqstarts <- start(seqnames(ends))[findPart(fhits, b)]
    seqstarts[is.na(seqstarts)] <- 1L
    fhits[is.na(fhits)] <- 0L
    fwindows <- restrict(IRanges(end=fhits, width = k), seqstarts, keep.all.ranges=TRUE)
    
    dist <- pc(extractList(start(starts), pwindows) - end(query), ### pwindows=5
               end(query) - extractList(end(ends), fwindows)) ### Try to write another function that disregards preceding/follows or somethingelse, just look at absolute distance regarding TSS. 

    ans <- pc(extractList(start_ord, pwindows), extractList(end_ord, fwindows))

    if (!ignore.overlaps) {
        hits <- findOverlaps(query, subject, ignore.strand=ignore.strand)
        hitsList <- as(hits, "List")
        dist <- pc(dist, relist(rep(0L, length(hits)), hitsList))
        ans <- pc(ans, hitsList)
    }

    ans[heads(order(dist), k)]
}