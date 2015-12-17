#Author: Tomislav Ilicic
#Date: 2/07/15
#Organisation: WTSI
#Description: SCRIPT TO EXTRACT BIOLOGICAL AND TECHNICAL FEATURES FORM HT-SEQ COUNT DATA. 

#Code reviewed and edited by Davis McCarthy, December 2015

#' Extracts biological and technical features for given dataset
#' 
#' @param counts_nm Gene expression counts dataframe (genes x cells). Either 
#' normalised by library size or TPM values
#' @param read_metrics Dataframe with mapping statistics produced by python 
#' pipeline
#' @param prefix Prefix of outputfiles
#' @param output_dir Output directory of files
#' @param common_features Subset of features that are applicable within one 
#' species, but across cell types
#' @param GO_terms DataFrame with gene ontology term IDs, that will be used in 
#' feature extraction
#' @param extra_genes Additional genes used for feature extraction
#' @param Organism The target organism to generate the features for
#' 
#' @details This function takes a combination of gene counts and mapping 
#' statistics to extract biological and technical features, which than can be 
#' used for quality data analysis
#' 
#' @return a list with two elements, one providing all features, and one 
#' providing common features.
#' 
#' @export
#' 
extract_features <- function(counts_nm, read_metrics, prefix, output_dir, 
                             common_features, GO_terms, extra_genes, organism) {
    ## define output directory and create if needed
    info("Extracting features")
    o <- paste(output_dir, prefix, sep = "/")
    print(output_dir)
    dir.create(o, showWarnings = TRUE, recursive = TRUE)
    ## define genes
    genes <- rownames(counts_nm)
    if (is.null(genes) | length(genes) == 0) {
        info("Please annotate your expression matrix with genes identifiers as rownames")
        return(NULL)
    }
    #GENERATE ALL FEATURES       
    features_all <- feature_generation(counts_nm, read_metrics, GO_terms, 
                                       extra_genes, organism)
    info(paste0("FEATURES generated."))
    sds <- apply(features_all, 2, sd)
    #REMOVE 0-VARIANCE VALUE FEATURES
    features_all <- features_all[,sds != 0]
    types <- c("all", "common")
    ## define common features
    features_common <- features_all[, which(colnames(features_all) %in% 
                                                common_features)]
    ## write features to file
    ## NB may be better to use file.path() for defining output files
    f_all <- paste0(o, "/", prefix, ".", types[1], ".features")
    f_common <- paste0(o, "/", prefix, ".", types[2], ".features")
    write.table(features_common, f_common)
    write.table(features_all, f_all)
    info(paste0("Features saved: ", f_all))
    info(paste0("Features saved: ", f_common))
    ## return features in list
    return(list(features_all, features_common))
}

################################################################################

#' Helper Function to create all features
#' 
#' @param counts_nm Gene expression counts dataframe (genes x cells). Either 
#' normalised by library size or TPM values
#' @param read_metrics Dataframe with mapping statistics produced by python 
#' pipeline
#' @param common_features Subset of features that are applicable within one 
#' species, but across cell types
#' @param GO_terms DataFrame with gene ontology term IDs, that will be used in 
#' feature extraction
#' @param extra_genes Additional genes used for feature extraction
#' @param Organism The target organism to generate the features for
#' 
#' @importFrom AnnotationDbi Term
#' @importFrom topGO annFUN.org
#' @export
#' 
#' @return Returns the entire set of features in a data.frame
#'
#' 
feature_generation <- function(counts_nm, read_metrics, GO_terms, extra_genes, 
                               organism) {
    ## initialise features list  
    features <- list()
    
    #REMOVE ALL 0 GENES
    counts_nm <- data.frame(counts_nm)
    genes_mean <- rowMeans(counts_nm)
    genes_zero <- which(genes_mean == 0)
    genes_mean <- genes_mean[-genes_zero]
    counts_nm_mean <- counts_nm[-genes_zero,] / genes_mean
    
    ########################################
    ####TECHINCAL FEATURES##################
    ########################################
    
    ## only consider reads mapped to genes (excluding ERCCs)
    ercc_counts = read_metrics$ercc
    if ( is.null(ercc_counts) ) {
        ercc_counts <- 0
    }
    number_mapped_reads_prop <- ((read_metrics$mapped - ercc_counts) / 
                                     read_metrics$total)
    multi_mapped_prop <- read_metrics[, 4] / read_metrics$mapped
    
    #HOPE THAT THIS WORKS ALSO WHEN REGRESSION NORMALIZATION HAS BEEN APPLIED
    #AS SOME REGRESSION METHODS ALSO PUSH 0 TO ANOTHER VALUE
    detected_genes <- apply(counts_nm, 2, function(x) {
        return(length(which(x > 0)))
    })
    
    counts_nm_mean_log <- log(counts_nm_mean + 0.001)
    genes_var <- apply(counts_nm_mean_log,1,var)
    genes_means_log <- log(genes_mean)
    
    #SELECT ONLY HIGHLY VARIABLES GENES WITH STRONG EXPRESSION BASED ON LAST QUANTILE
    i <- which(genes_var > quantile(genes_var)[4] & 
                   genes_means_log > quantile(genes_means_log)[4])
    counts_nm_mean_log_high_var_mean <- counts_nm_mean_log[i,]
    
    #VARIANCE ACROSS HIGHLY EXPRESSED GENES
    transcriptome_variance <- apply(counts_nm_mean_log_high_var_mean, 2, var)
    
    #NUMBER OF HIGHLY EXPRESSED GENES
    m <- quantile(counts_nm_mean[i,2])[4]
    number_of_highly_expressed_variable_genes <- 
        apply(counts_nm_mean[i,], 2, function(x) {return(sum(x > m))})
    
    #NUMBER OF LOW TO HIGH EXPRESED AND VARIABLE GENES PER INTERVAL
    num_of_high_var_exp_genes_interval <-
        apply(counts_nm_mean_log_high_var_mean, 2, function(x) {
            hst <- hist(x, breaks = c(-1000, -4, -2, 0, 2, 1000), plot = FALSE)
            hst$counts
        })
    num_of_high_var_exp_genes_interval <- t(num_of_high_var_exp_genes_interval)
    colnames(num_of_high_var_exp_genes_interval) <- 
        paste0("num_of_high_var_exp_genes_interval_", 
               1:ncol(num_of_high_var_exp_genes_interval))
    
    cell_to_mean_corr <- cor(counts_nm, rowMeans(counts_nm),
                             method = "spearman")
    
    techincal_features <- cbind(read_metrics$total, number_mapped_reads_prop,
                                multi_mapped_prop, detected_genes, 
                                transcriptome_variance,
                                num_of_high_var_exp_genes_interval, 
                                cell_to_mean_corr,
                                number_of_highly_expressed_variable_genes)
    tech_names <- c("#Total reads", "Mapped %", "Multi-mapped %", 
                    "#Detected genes", "Transcriptome variance", 
                    paste0("High expr.var genes intv.", 
                           1:ncol(num_of_high_var_exp_genes_interval)), 
                    "Cell-to-mean corr", "#High exp.var genes")
    
    if ( ncol(read_metrics) > 23 ) {
        additional_mapping_stats <- read_metrics[, 24:(ncol(read_metrics))]
        additional_mapping_stats_prop <- (additional_mapping_stats / 
                                              read_metrics$mapped)
        
        #ADDITIONAL MAPPING STATISTICS PROVIDED BY HT-SEQ
        additional_mapping_stats_prop_names <- sapply(
            colnames(additional_mapping_stats_prop), function(x) {
                if ( x == "__not_aligned") return("Unmapped reads %")
                if ( x == "__no_feature") return("Non-exonic reads %")
                if ( x == "__ambiguous") return( "Ambigious-gene reads %")
                if ( x == "__alignment_not_unique") { 
                    return( "Multi-mapped reads (HTseq) %")
                }  
                if ( x == "__too_low_aQual" ) return( "Low aligQual %")
                if ( x == "ercc" ) return( "ERCC/mapped reads")
            })
        techincal_features <- cbind(techincal_features, 
                                    additional_mapping_stats_prop)
        tech_names <- c(tech_names, additional_mapping_stats_prop_names)
    }
    
    colnames(techincal_features) <- tech_names
    
    ########################################
    ####BIOLOGICAL FEATURES##################
    ########################################
    
    #ASSUME IT IS MOUSE
    GO_BP <- topGO::annFUN.org("BP", mapping = "org.Mm.eg.db", ID = "ensembl")
    GO_CC <- topGO::annFUN.org("CC", mapping = "org.Mm.eg.db", ID = "ensembl")
    
    if ( organism == "human" ) {
        GO_BP <- topGO::annFUN.org("BP", mapping = "org.Hs.eg.db", 
                                   ID = "ensembl")  
        GO_CC <- topGO::annFUN.org("CC", mapping = "org.Hs.eg.db", 
                                   ID = "ensembl")
    }
    
    GO <- c(GO_BP, GO_CC)
    
    #PROPORTION OF MAPPED READS MAPPED TO THE GO TERM
    go_prop <- sapply(unlist(GO_terms), function(go_id) {
        prop <- sum_prop(counts_nm, unlist(GO[go_id])) 
        return(prop)
    }, simplify = FALSE)
    go_prop <- do.call(cbind, go_prop)
    go_names <- AmnnotationDbi::Term(unlist(GO_terms))
    go_names <- sapply(go_names, simple_cap)
    colnames(go_prop) <- go_names 
    
    #PROPORTION OF MAPPED READS MAPPED TO SPECIFIC GENES
    extra_genes_prop <- sapply(extra_genes, function(extra_g) {
        #print(extra_g[1])
        prop <- sum_prop(counts_nm, extra_g[-1]) 
        return(prop)
    }, simplify = FALSE)
    extra_genes_prop <- do.call(cbind, extra_genes_prop)
    extra_genes_colnames <- sapply(extra_genes, function(extra_g) {
        return(extra_g[1])
    }, simplify = FALSE)
    colnames(extra_genes_prop) <- unlist(extra_genes_colnames)
    
    biological_features <- cbind(go_prop, extra_genes_prop)
    colnames(biological_features) <- paste0(colnames(biological_features), " %")
    
    features <- data.frame(techincal_features, biological_features)
    colnames(features) <- c(colnames(techincal_features), 
                            colnames(biological_features))
    rownames(features) <- colnames(counts_nm)
    return(features)
}

################################################################################
## sum_prop

#' Sums up normalised values of genes to groups. 
#' 
#' Supports TPM and proportion of mapped reads.
#' 
#' @param counts Normalised gene expression count matrix
#' @param gene_interest dataframe of genes of interest to merge
#' 
#' @export
#' 
sum_prop <- function(counts, genes_interest) {
    genes_interest_i <- which(rownames(counts) %in% unlist(genes_interest))
    genes_interest_counts_prop <- colSums(counts[genes_interest_i,])
    return(genes_interest_counts_prop)
}

################################################################################
## simple_cap

#' Converts all first letters to capital letters
#' 
#' @param x string
#' @export
#' 
simple_cap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep = "", collapse = " ")
}

################################################################################
## assess_cell_quality_SVM

#' Assess quality of a cell - SVM version
#' 
#' @param training_set_features A training set containing features (cells x 
#' features) for prediction
#' @param training_set_labels Annotation of each individual cell if high or low quality (1 or 0 respectively)
#' @param test_set_features Dataset to predict containing features (cells x 
#' features)
#' @param ensemble_param Dataframe of parameters for SVM
#' @return Returns a dataframe indicating which cell is low or high quality (0 
#' or 1 respectively)
#' 
#' @details This function takes a traning set + annotation to predict a test 
#' set. It requires that hyper-parameters have been optimised. 
#' 
#' @return data.frame with decision on quality of cells
#' 
#' @importFrom e1071 svm
#' @export
#' 
assess_cell_quality_SVM <- function(training_set_features, training_set_labels,
                                    test_set_features, ensemble_param) {
  
  data_set <- data.frame(l = training_set_labels, 
                       unlist(as.matrix(training_set_features)))
  data_set$l <- as.factor(data_set$l)
  test_data <- data.frame(as.matrix(test_set_features))
  form <- formula("l ~ .")
  
  final_results <- sapply(1:nrow(ensemble_param), function(x) {
      parameters <- ensemble_param[x,]
      weights <- table(data_set$l) 
      weights[1] <- parameters[3]
      weights[2] <- 1
      kernel <- "radial"
      
      model <- e1071::svm(form, data = data_set, gamma = parameters[1], 
                          cost = parameters[2], kernel = kernel, 
                          class.weights = weights)
      
      pred_test <- predict(model, test_data)
      svm_test <- as.numeric(levels(pred_test))[pred_test]
      return(svm_test)
  }, simplify = FALSE)
  final_results <- do.call(cbind, final_results)
  
  #Voting scheme to determine final label
  final <- vote(final_results)
  final_df <- data.frame(cell = rownames(test_set_features), quality = final)
  return(final_df)
}

################################################################################
## assess_cell_quality_PCA

#' ASSESS CELL QUALITY USING PCA AND OUTLIER DETECTION
#' 
#' @param test_set_features Input dataset containing features (cell x features)
#' @param output_dir Output directory
#' @param prefix Prefix of output 
#' 
#' @details This function applies PCA on features and uses outlier detection to
#'  determine which cells are low and which are high quality
#' @return Returns a dataframe indicating which cell is low or high quality (0 
#' or 1 respectively)
#' 
#' @importFrom mvoutlier pcout
#' @importFrom mvoutlier uni.plot
#' @export
#' 
assess_cell_quality_PCA <- function(test_set_features, output_dir, prefix) {
    ## perform PCA
    pca <- prcomp(test_set_features, scale = TRUE, center = TRUE)
    pca_var_explained <- summary(pca)
    pcout_c <- mvoutlier::pcout(pca$x[, 1:2])
    dimens <- min(10, ncol(test_set_features))
    low_qual_i <- which(pcout_c$wfinal01 == 0)
    #uni_2 <- mvoutlier::uni.plot(pca$x[, 1:dimens])
    #low_qual_i = which(uni_2$outliers == TRUE)
    
    mtdna <- t.test(test_set_features[,6][low_qual_i], 
                    test_set_features[,6][-low_qual_i], 
                    alternative = "greater")$p.value
    mapped_prop <- t.test(test_set_features[,2][low_qual_i], 
                          test_set_features[,2][-low_qual_i], 
                          alternative = "less")$p.value
    
    #Determine which cluster is low and which high quality
    if ((mapped_prop) > 0.5) {
        types <- rep(0, nrow(test_set_features))
        types[low_qual_i] <- 1
        
    } else{
        types <- rep(1, nrow(test_set_features))
        types[low_qual_i] = 0
    }
    
    ## define data frame with cell types
    data_frame <- data.frame(type = as.character(types), pca$x)  
    col <- c("0" = "red","1" = "darkgreen")
    
    ## plot PCs
    plot <- ggplot2::ggplot(data_frame, aes(x = PC1, y = PC2)) + 
        geom_point(aes(colour = type)) + 
        scale_colour_manual(values = col) +
        theme_bw()  + 
        theme(axis.line = element_blank(), axis.text.x = element_text(), 
              axis.text.y = element_text(), axis.ticks.length = unit(0, "mm"),
              axis.title.x = element_blank(), axis.title.y = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(fill = NA, color = "black", 
                                          linetype = "solid"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background = element_blank()) 
    
    #Plot top 3 features for PC1
    top5_PC1_i <- order(abs(pca$rotation[,1]), decreasing = TRUE)
    top5_f_pc1 <- names(pca$rotation[(top5_PC1_i[1:3]), 1])
    names <- colnames(test_set_features)
    test_set_features_top5_pc1_i <- which(names %in% top5_f_pc1)
    size <- 0.5
    text_size <- 12
    border_size <- 1
    plotsPC1 <- sapply(test_set_features_top5_pc1_i, function(f) {
        feature <- test_set_features[,f]
        df <- data.frame(counts = log(feature + 0.0001), type = types)
        plot <- ggplot2::ggplot(df, aes(x = type)) + 
            geom_boxplot(aes(colour = factor(type), y = counts), trim = TRUE, 
                         alpha = 0.3, adjust = 1, size = size, 
                         outlier.size = 0) +
            ggtitle(names[f]) + 
            theme_bw() + 
            theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                  axis.text.y = element_text(size = text_size),
                  axis.ticks.length = unit(0, "mm"), 
                  axis.title.x = element_blank(), axis.title.y = element_blank(),
                  legend.position = "none",
                  plot.title = element_text(size = text_size),
                  panel.background = element_blank(),
                  panel.border = element_rect(
                      fill = NA, color = "black", size = border_size, 
                      linetype = "solid"), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_blank()) + 
            scale_color_manual(values = col)
        return(plot)
    }, simplify = FALSE)
    
    #Plot top 3 features for PC2
    top5_PC2_i <- order(abs(pca$rotation[,2]), decreasing = TRUE)
    top5_f_pc2 <- names(pca$rotation[(top5_PC2_i[1:3]),2])
    names <- colnames(test_set_features)
    test_set_features_top5_pc2_i <- which(names %in% top5_f_pc2)
    
    plotPC2 <- sapply(test_set_features_top5_pc2_i, function(f) {
        feature <- test_set_features[,f]
        df <- data.frame(counts = log(feature + 0.0001), type = types)
        plot <- ggplot2::ggplot(df, aes(x = type)) + 
            geom_boxplot(aes(colour = factor(type), y = counts), trim = TRUE, 
                         alpha = 0.3, adjust = 1, size = size, 
                         outlier.size = 0) +
            ggtitle(names[f]) + 
            theme_bw() + 
            theme(axis.line = element_blank(), axis.text.x = element_blank(),
                  axis.text.y = element_text(size = text_size),
                  axis.ticks.length = unit(0, "mm"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none",
                  panel.background = element_blank(),
                  plot.title = element_text(size = text_size),
                  panel.border = element_rect(fill = NA, color = "black", 
                                              size = border_size, 
                                             linetype = "solid"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.background = element_blank()) + 
            scale_color_manual(values = col)
        return(plot)
    }, simplify = FALSE)
    
    #ARRANGE TOP FEATURES ONTO A GRID
    #PLOT PCA IN THE MIDDLE AND FEATURES LEFT AND BOTTOM
    l <- matrix(c(2, 2, 3, 3, 4, 4, rep(1,5), 5, rep(1,5), 5, rep(1,5), 6, 
                  rep(1,5), 6,rep(1,5), 7), nrow = 6)
    l <- cbind(l, c(rep(1, 5), 7))
    
    tp <- matrix(c(2, 2, 3, 3, 4, 4, 10, rep(1, 5), 5, 5, rep(1, 5), 5, 5,
                  rep(1, 5), 6, 6, rep(1, 5), 6, 6, rep(1, 5), 7, 7), nrow = 7)
    tp <- cbind(t, c(rep(1, 5), 7, 7))
    
    o <- paste0(output_dir, "/", prefix, "_pca_quality_check.pdf")
    pdf(o, width = 10, height = 7)
    multiplot(plot,  plotlist = c(plotsPC1, plotPC2), layout = l)
    dev.off()
    annot <- data.frame(cell = rownames(test_set_features), quality = types)
    return(annot)
}

################################################################################
## multiplot

#' Internal multiplot function to combine plots onto a grid
#' 
#' @param ... individual plots to combine into a single plot
#' @param plotlist a vector with names of plots to use in the plot
#' @param file string giving filename to which pdf of plots is to be saved
#' @param cols integer giving number of columns for the plot
#' @param layout matrix defining the layout for the plots
#' 
#' @returns a plot object
#' 
multiplot <- function(..., plotlist = NULL, file, cols = 6, layout = NULL) {
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(
            layout = grid::grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], 
                  vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
        }
    }
}

#' Internal voting function to get final labels
#' 
#' @param predicitions Predicted labels
#' 
vote <- function(predictions) {
    freq <- apply(predictions, 1, function(x) {
        f <- table(x)
        return(names(f)[which.max(f)])
    })
    
    freq <- as.numeric(freq)
    return(freq)
}

################################################################################
## info

#' Internal function to print info string
#' 
#' @param string a string to print as an info message
#' 
info <- function(string) { 
    print(paste0("[INFO]:", string))
}

################################################################################
## normalise_by_factor

#' Internal function to normalize by library size
#' 
#' @param counts matrix of counts
#' @param norm_factor vector of normalisation factors
#' 
normalise_by_factor <- function(counts, norm_factor) { 
    t(t(counts) / norm_factor)
}
