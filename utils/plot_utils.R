### BASIC HELPER FUNCTIONS ###
set_colors <- function(df, col, colors="default"){
    if(colors == "default"){
        n <- length(unique(df[,col]))
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        set.seed(10)
        cols=sample(col_vector, n)
    } else {
        cols=colors
    }
    return(cols)
}


### PART 1: RELATED PLOTTING FUNCTIONS FOR TOTAL CROSS-LOCI CORRELATIONS

boxplot_total_cross_loci_corr <- function(corr_df, mark_col, method_col, corr_col, x_name, y_name, w, h, colors, plot_proj_name, filename, num_row, num_col, main){
    marks <- unique(corr_df[,mark_col])
    methods <- unique(corr_df[,method_col])

    colors = set_colors(corr_df, method_col, colors)

    p1 <- ggplot(corr_df, aes(x=as.factor(eval(as.name(paste(method_col)))),y=eval(as.name(corr_col)),colour=factor(eval(as.name(paste(method_col))))))+
    geom_boxplot() + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    facet_wrap(facets = vars(factor(eval(as.name(paste(mark_col))), levels=marks)), nrow=num_row, ncol=num_col) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),) +
    guides(col= guide_legend(title= "Method"))

    dir.create(sprintf("./plots/%s", plot_proj_name), recursive=TRUE)
    ggsave(plot=p1, filename = sprintf("./plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}

### RELATED PLOTTING FUNCTIONS FOR CT-SPECIFIC CROSS-LOCI CORRELATIONS
smoothplot_grid_ct_spec_cross_loci <- function(corr_df, main, x, y, mark_col, x_name, y_name, plot_proj_name, filename, colors, method_col="Method", num_row, num_col, w=11, h=5){
    marks <- unique(corr_df[,mark_col])
    methods <- unique(corr_df[,method_col])
    colors = set_colors(corr_df, method_col, colors)

    p1 <- ggplot(corr_df, aes(x=as.numeric(eval(as.name(x))),y=eval(as.name(y)),colour=factor(eval(as.name(paste(method_col))))))+
    geom_smooth(se=FALSE) + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    facet_wrap(facets = vars(factor(eval(as.name(paste(mark_col))), levels=marks)), nrow=num_row, ncol=num_col) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    guides(col= guide_legend(title= "Method"))

    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}

### PART 2: 

###############################
plot_cross_loci_ct_spec_smooth <- function(cross_loci_corr_df, main, x, y, x_name, y_name, plot_proj_name, filename, colors, col_name_method="method",w=11, h=5){
    p1 <- ggplot(cross_loci_corr_df, aes(x=as.numeric(eval(as.name(x))),y=eval(as.name(y)),colour=factor(eval(as.name(paste(col_name_method))))))+
    geom_smooth(se=FALSE) + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    guides(col= guide_legend(title= "Method"))

    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}

##############################
plot_boxplot_dist <- function(df_corr, colors, pred_list, method, title, x_name, y_name, filename, w=20, h=6){
    df_full <- df_corr
    df_filter <- df_full[complete.cases(df_full),]
    df_long <- as.data.frame(pivot_longer(df_full, cols=all_of(pred_list), names_to="Feature", values_to=method))

    df_long_ss <- df_long[df_long$Feature %in% pred_list,]
    p1 <- ggplot(df_long_ss, aes(x=as.factor(dist),y=eval(as.name(paste(method))),colour=factor(Feature)))+
        #geom_point()+
        geom_boxplot() + 
        labs(title=title) +
        labs(x = x_name, y = y_name) +
        theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
        guides(col= guide_legend(title= "Method")) +
        scale_colour_manual(values=colors)

    ggsave(plot=p1, filename = filename, width = w, height = h)
}

metaplot_dist_individual_ct <- function(celltypes, df_long, meta, exp, chip, l=10000, res=200, plot_proj_name, colors, method="spearman", file_ext, pred_list){
    for (celltype in celltypes[1:length(celltypes)]){
        knn_features_ct_ss <- df_long[df_long$Cell_type == celltype,]
        gr <- granges(chip)
        mcols(gr) <- knn_features_ct_ss
        knn_all_features_ct_ss <- gr

        test_exp <- exp[,celltype, drop = FALSE]
        test_exp$gene_name <- row.names(test_exp)

        tss <- reduce_to_tss(meta)
        tss_extended <- extend(tss, upstream = 500, downstream = 500)
        common_genes <- intersect(tss$gene_name, rownames(test_exp))
        tss <- tss[tss$gene_name %in% common_genes]
        test_exp <- test_exp[common_genes,]
        message(paste("Row names of test_exp and gene names are identical:", identical(rownames(test_exp), tss$gene_name)))

        exp_granges <- tss

        message(celltype)
        test_chip <- chip[,celltype]
        seqlevels(exp_granges) <- seqlevels(test_chip)
        seqlengths(test_chip) <- seqlengths(exp_granges)
        seqlengths(knn_all_features_ct_ss) <- seqlengths(exp_granges)

        n_w <- ((l/res)*2+1)
        exp_granges <- extend(exp_granges, upstream = l, downstream = l)
        overlap <- data.frame(findOverlaps(exp_granges, test_chip, type="any"))

        overlap_pred <- data.frame(findOverlaps(exp_granges, knn_all_features_ct_ss, type="any"))

        chip_df <- as.data.frame(test_chip)
        knn_all_features_ct_ss_df <- as.data.frame(knn_all_features_ct_ss)

        vector_of_scores <- chip_df[,celltype][overlap$subjectHits]
        scores_around_tss_df <- do.call(rbind.data.frame, split(vector_of_scores, ceiling(seq_along(vector_of_scores)/n_w)))
        names(scores_around_tss_df) <- c(paste0('w', seq(-1*((n_w-1)/2),0,1)), paste0('w', seq(1,(n_w-1)/2,1)))

        features <- pred_list
        cor_data_to_plot_feature_l <- lapply(features, function(feature){
            vector_of_scores_feature <- knn_all_features_ct_ss_df[,feature][overlap_pred$subjectHits]
            scores_around_tss_df_feature <- do.call(rbind.data.frame, split(vector_of_scores_feature, ceiling(seq_along(vector_of_scores_feature)/n_w)))
            names(scores_around_tss_df_feature) <- c(paste0('w', seq(-1*((n_w-1)/2),0,1)), paste0('w', seq(1,(n_w-1)/2,1)))

            corr_values_true_feature <- sapply(names(scores_around_tss_df), function(x){
                cor(scores_around_tss_df[,x], scores_around_tss_df_feature[,x], method = method)
            })

            cor_data_to_plot_feature <- data.frame(Window=names(scores_around_tss_df_feature), Feature=corr_values_true_feature)
            names(cor_data_to_plot_feature)[2] <- feature
            cor_data_to_plot_feature <- melt(cor_data_to_plot_feature, id.vars = 'Window')
            cor_data_to_plot_feature$Window <- factor(cor_data_to_plot_feature$Window, levels=names(scores_around_tss_df_feature))

            cor_data_to_plot_feature
        })

        cor_data_full <- do.call(rbind, cor_data_to_plot_feature_l)
        saveRDS(cor_data_full, file.path(sprintf("plots/%s/", plot_proj_name), sprintf("cor_data_full_%s_%s.rds", celltype, file_ext)))

        cor_data_full <- readRDS(file.path(sprintf("plots/%s/", plot_proj_name), sprintf("cor_data_full_%s_%s.rds", celltype, file_ext)))
        features_ss <- pred_list
        cor_data_ss <- cor_data_full[cor_data_full$variable %in% features_ss,]

        plot_dist_to_tss(cor_data_full, method=method, proj_name=plot_proj_name, ct=celltype, l=l, file_ext=file_ext)
        plot_dist_to_tss(cor_data_ss, method=method, proj_name=plot_proj_name, ct=celltype, l=l, file_ext=paste0(file_ext, "-ss"))
    }
}

metaplot_dist_per_loci <- function(df_long, meta, exp, chip, l=10000, res=200, plot_proj_name, colors, method="spearman", file_ext, pred_list){
    celltypes = intersect(names(mcols(chip)), names(exp))
    celltype=celltypes[1]

    # knn_features_ct_ss <- df_long[df_long$locus %in% loci,]
    knn_features_ct_ss <- df_long
    gr <- granges(chip)
    mcols(gr) <- knn_features_ct_ss
    knn_all_features_ct_ss <- gr

    test_exp <- exp[,celltype, drop=FALSE]
    test_exp$gene_name <- row.names(test_exp)

    tss <- reduce_to_tss(meta)
    tss_extended <- extend(tss, upstream = 500, downstream = 500)
    common_genes <- intersect(tss$gene_name, rownames(test_exp))
    tss <- tss[tss$gene_name %in% common_genes]
    test_exp <- test_exp[common_genes,]
    message(paste("Row names of test_exp and gene names are identical:", identical(rownames(test_exp), tss$gene_name)))

    exp_granges <- tss

    message(celltype)
    test_chip <- chip[,celltype]
    seqlevels(exp_granges) <- seqlevels(test_chip)
    seqlengths(test_chip) <- seqlengths(exp_granges)
    seqlengths(knn_all_features_ct_ss) <- seqlengths(exp_granges)

    n_w <- ((l/res)*2+1)
    exp_granges <- extend(exp_granges, upstream = l, downstream = l)
    overlap <- data.frame(findOverlaps(exp_granges, test_chip, type="any"))
    overlap_pred <- data.frame(findOverlaps(exp_granges, knn_all_features_ct_ss, type="any"))

    chip_df <- as.data.frame(test_chip)
    knn_all_features_ct_ss_df <- as.data.frame(knn_all_features_ct_ss)

    vector_of_loci <- overlap$subjectHits
    loci_around_tss_df <- do.call(rbind.data.frame, split(vector_of_loci, ceiling(seq_along(vector_of_loci)/n_w)))
    names(loci_around_tss_df) <- c(paste0('w', seq(-1*((n_w-1)/2),0,1)), paste0('w', seq(1,(n_w-1)/2,1)))

    loci_around_tss_df$group = seq(1, nrow(loci_around_tss_df))
    loci_around_tss_df_long = loci_around_tss_df %>% pivot_longer(!group, names_to="Window", values_to="locus") %>% as.data.frame()

    knn_features_ct_ss_merged = merge(knn_features_ct_ss, loci_around_tss_df_long, by="locus")
    saveRDS(knn_features_ct_ss_merged, file.path(sprintf("plots/%s/", plot_proj_name), sprintf("cor_ct_data_full_%s_%s.rds", "all", file_ext)))

    knn_features_ct_ss_merged <- readRDS(file.path(sprintf("plots/%s/", plot_proj_name), sprintf("cor_ct_data_full_%s_%s.rds", "all", file_ext)))
    knn_features_ct_ss_merged_grouped <- knn_features_ct_ss_merged %>% group_by(Window, method) %>% dplyr::summarize(mean_cor=mean(value)) %>% as.data.frame()

    w_order = c(paste0('w', seq(-1*((n_w-1)/2),0,1)), paste0('w', seq(1,(n_w-1)/2,1)))
    plot_dist_to_tss(knn_features_ct_ss_merged_grouped, level_order = w_order, method=method, variable = "method", score="mean_cor", proj_name=plot_proj_name, ct="all", l=l, file_ext=file_ext)
}

plot_dist_to_tss <- function(df, method, level_order, variable, proj_name, ct, l=10000, file_ext, score="value"){
    p1 <- ggplot(df, aes(factor(Window, levels=level_order), eval(as.name(paste(score))))) +
        geom_line(aes(group=eval(as.name(paste(variable))), colour=eval(as.name(paste(variable))))) +
        labs(y=paste0(method,' Correlation'), x='Bin', title=paste0('Correlation of predicted signal intensity with true signal intensity as a function of distance to TSS, ', ct), subtitle=paste0('Window = ', l)) +
        theme_bw(base_size = 16) +
        theme(panel.grid = element_blank(), 
            title = element_text(size=12),  
            axis.text.x = element_text(angle = 90, hjust = 1, size=4), 
            legend.position = c(0.8,0.8)) +
        scale_colour_manual(values=colors)
    ggsave(plot=p1, filename=paste0(getwd(), sprintf("/plots/%s/metaplot_corr_w_true_combined_w_int/", proj_name), paste0('Window=', l, ", ", ct, sprintf("_all_methods_%s_%s", method, file_ext)), '.png'), width=12, height=9)
}

plot_p_val <- function(p_val_df, main, x, y, x_name, y_name, plot_proj_name, filename, colors, w=11, h=5, rotate=90, scale_max=999){
    if(scale_max == 999){
        scale_max = max(p_val_df[,y])
    } else {
        p_val_df[,y][p_val_df[,y] > scale_max] = scale_max
    }
    
    p1 <- ggplot(p_val_df, aes(x=as.factor(eval(as.name(x))),y=as.numeric(eval(as.name(y)))))+
    geom_point() + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    guides(col= guide_legend(title= "Method")) +
    theme(axis.text.x = element_text(angle = rotate, vjust = 0.5, hjust=1))+
    scale_y_continuous(limits = c(0, scale_max))
    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}


plot_cross_loci_ct_spec_smooth <- function(cross_loci_corr_df, main, x, y, x_name, y_name, plot_proj_name, filename, colors, col_name_method="method",w=11, h=5){
    p1 <- ggplot(cross_loci_corr_df, aes(x=as.numeric(eval(as.name(x))),y=eval(as.name(y)),colour=factor(eval(as.name(paste(col_name_method))))))+
    geom_smooth(se=FALSE) + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    guides(col= guide_legend(title= "Method"))

    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}



plot_cross_loci_ct_spec_box <- function(cross_loci_corr_df, main, x, y, x_name, y_name, plot_proj_name, filename, colors="default", col_name_method="method",w=11, h=5){

    colors = set_colors(cross_loci_corr_df, x, colors)

    p1 <- ggplot(cross_loci_corr_df, aes(x=as.factor(eval(as.name(x))),y=eval(as.name(y)),colour=factor(eval(as.name(paste(col_name_method))))))+
    geom_boxplot() + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    guides(col= guide_legend(title= "Method"))

    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}

plot_cross_loci_ct_spec_box_grid <- function(cross_loci_corr_df, main, grid_col, grid_col_level, x, y, x_name, y_name, plot_proj_name, filename, colors="default", w=11, h=5, num_row=10, num_col=10){

    colors = set_colors(cross_loci_corr_df, x, colors)

    p1 <- ggplot(cross_loci_corr_df, aes(x=as.factor(eval(as.name(paste(x)))),y=eval(as.name(y)),colour=factor(eval(as.name(paste(x))))))+
    geom_boxplot() + 
    labs(title=main) +
    scale_colour_manual(
        values = colors,
        aesthetics = c("colour")) +
    facet_wrap(facets = vars(factor(eval(as.name(paste(grid_col))), levels=grid_col_level)), nrow=num_row, ncol=num_col) +
    labs(x = x_name, y = y_name) +
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    guides(col= guide_legend(title= "Method"))

    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}

# pie(rep(1,n), col=sample(col_vector, n))

plot_cross_loci_general_box <- function(cross_loci_corr_df, main, x, y, level_order, x_name, y_name, plot_proj_name, filename, colors="default", col_name_method="method",w=11, h=5){
    if(colors == "default"){
        n <- length(unique(cross_loci_corr_df[,x]))
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        set.seed(10)
        cols=sample(col_vector, n)
    } else {
        cols=colors
    }


    p1 <- ggplot(cross_loci_corr_df, aes(x=factor(eval(as.name(x)), levels=level_order),y=eval(as.name(y)),fill=factor(eval(as.name(x)), levels=level_order)))+
        geom_boxplot() + 
        labs(title=main) +
        stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
        labs(x = x_name, y = y_name) +
        theme(axis.text.x = element_text(angle = 30), axis.title.x = element_blank(), legend.position = "none", plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
        guides(col= guide_legend(title= "Method")) +
        scale_colour_manual(
            values = cols,
            aesthetics = c("fill"))

    dir.create(sprintf("plots/%s/", plot_proj_name), recursive=TRUE)
    ggsave(plot=p1, filename = sprintf("plots/%s/%s", plot_proj_name, filename), width = w, height = h)
}

plot_cross_ct_ct_spec_box <- function(cross_ct_corr_df, main, x, y, x_name, y_name, plot_proj_name, filename, colors, col_name_method="method",w=11, h=5){
    p2 <- ggplot(pcc_df_acr_loci, aes(x=factor(method, levels=level_order),y=p_corr,fill=factor(method, levels=level_order)))+
    geom_boxplot() + 
    labs(title="Cross-loci correlation Pred-True Before and After Deconvolution (Using 30 training celltypes)") +
    facet_grid(cols = vars(factor(decon, levels=c("Before Deconvolution", "After Deconvolution")))) +
    stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
    theme(axis.text.x = element_text(angle = 30), axis.title.x = element_blank(), legend.position = "none", plot.margin = unit(c(0.5,0.5,0,0.5), "cm")) +
    labs(x = "Methods", y = "Pearson Correlation") +
    scale_colour_manual(
        aesthetics = c("fill")) +
    labs(x = "Methods", y = "Pearson Correlation")
}



