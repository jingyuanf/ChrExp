## Functions
subset_to_chr <- function(gr, chr){
  if (chr == "all"){
    gr = gr
  } else {
    gr = gr[seqnames(gr) == chr]
  }
  return(gr)
}

make_min_zero <- function(gr){
  mat <- as.matrix(mcols(gr))
  mat[mat < 0] <- 0
  mcols(gr) <- mat
  return(gr)
}

process_across <- function(across_dir, proj_name_across, chr){
  across_gr <- readRDS(file.path(across_dir, proj_name_across, "final-pred", "testing_pred_reg_full.rds"))
  across_gr <- make_min_zero(across_gr)
  across_gr <- subset_to_chr(across_gr, chr)

  dir.create(file.path(across_dir, proj_name_across, "final-pred", chr), recursive=TRUE)
  saveRDS(across_gr, file.path(across_dir, proj_name_across, "final-pred", chr, "across_gr.rds"))
  message("Across gr saved!")
  return(across_gr)
}

process_within <- function(pred_dir, proj_name_wit, chr){
  within_file_path <- file.path("bss-predictions", "within-model", proj_name_wit, "final_pred", "2", "with_distance", num_train_ct, "bw_files", "raw")
  within_files <- list.files(file.path(within_file_path), pattern = "*.bw", full.names = FALSE)

  within <- lapply(within_files, function(x){
    gr <- rtracklayer::import.bw(file.path(within_file_path, x), format = "bw")
    names(mcols(gr)) <- sapply(strsplit(x, "-", fixed=TRUE), function(x) (x[3]))
    gr
  })

  within_gr <- granges(within[[1]])
  mcols(within_gr) <- do.call(cbind, lapply(within, mcols))

  within_gr <- subset_to_chr(within_gr, chr)
  within_gr <- make_min_zero(within_gr)

  saveRDS(within_gr, file.path("bss-predictions", "within-model", proj_name_wit,"final_pred", "2", "with_distance", num_train_ct, "bw_files", sprintf("within_gr_log2_%s.rds", chr)))
  message("Within gr saved!")
  return(within_gr)
}

process_chrom <- function(chr_dir, proj_name_chr, chr){
  chrom_file_path <- file.path(chr_dir, proj_name_chr, "tracks", "raw")
  chrom_files <- list.files(file.path(chrom_file_path), pattern = "*.bw", full.names = FALSE)
  #across_files_raw <- list.files(file.path("across-model", proj_name, "predictions", proj_name, num_train_ct, "raw"), pattern = "*.bw", full.names = FALSE)

  chrom <- lapply(chrom_files, function(x){
    gr <- rtracklayer::import.bw(file.path(chrom_file_path, x), format = "bw")
    names(mcols(gr)) <- sapply(strsplit(x, "_", fixed=TRUE), function(x) (x[1])) # change to 5 for multlin setting
    gr
  })

  chrom_gr <- granges(chrom[[1]])
  mcols(chrom_gr) <- do.call(cbind, lapply(chrom, mcols))

  chrom_gr <- subset_to_chr(chrom_gr, chr)
  chrom_gr <- make_min_zero(chrom_gr)

  saveRDS(chrom_gr, file.path(chr_dir, proj_name_chr, sprintf("chrom_gr_%s.rds", chr)))
  message("Chrom gr saved!")
  return(chrom_gr)
}

process_int <- function(pred_dir, proj_name_int, chr, thrd){
  int_gr <- readRDS(sprintf("./bss-predictions/integration/%s/integration-acr-within-thrd-%s-new.rds", proj_name_int, thrd))

  int_gr <- subset_to_chr(int_gr, chr)
  int_gr <- make_min_zero(int_gr)
  return(int_gr)
}

deconvolve <- function(gr, cells_to_convolve, convolved){
  convolved_gr_new <- rowSums(as.data.frame(mcols(gr)[,cells_to_convolve]))
  gr_decon <- granges(gr)

  for(cell in cells_to_convolve){
    test <- mcols(gr)[,cell]
    test <- test / convolved_gr_new * convolved$score
    mcols(gr_decon) <- cbind(mcols(gr_decon), test)
  }

  df_gr_decon <- as.data.frame(mcols(gr_decon))
  df_gr_decon[is.na(df_gr_decon)] <- 0
  mcols(gr_decon) <- df_gr_decon
  names(mcols(gr_decon)) <- cells_to_convolve
  return(gr_decon)
}
