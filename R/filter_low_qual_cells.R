# Major functionality for the package

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
#' @param organism The target organism to generate the features for
#' 
#' @details This function takes a combination of gene counts and mapping 
#' statistics to extract biological and technical features, which than can be 
#' used for quality data analysis
#' 
#' @return a list with two elements, one providing all features, and one 
#' providing common features.
#' 
#' @importFrom utils data write.table
#' 
#' @export
#' @examples 
#' data(sample_counts)
#' data(sample_stats)
#' sample_counts_nm <- normalise_by_factor(sample_counts, colSums(sample_counts))
#' sample_features <- extract_features(sample_counts_nm, sample_stats)
extract_features <- function(counts_nm, read_metrics, prefix="", output_dir="", 
                             common_features=NULL, GO_terms=NULL, extra_genes=NULL, organism="mouse") {
    
    feature_info <- get("feature_info")
    if (is.null(common_features)) {
        common_features <- feature_info[[2]]
    }
    
    if (is.null(GO_terms)) {
        GO_terms <- feature_info[[1]]
    }
    
    if ( organism == "human" ||  organism == "org.Hs.eg.db") {
        organism <- "org.Hs.eg.db"
    } else if ( organism == "mouse" || organism == "org.Mm.eg.db") {
        organism <- "org.Mm.eg.db"
    } else{
        print.warnings("You have specified a different organism than mouse or human.\n
           This might work, but you have to make sure you have specified the appropiate database as organism (e.g. org.Hs.eg.db), and you have it also installed.\n 
           Also, pleae note that extra_genes need to match the organism of interest.")
    }
    
    if (is.null(extra_genes)) {
        if ( organism == "org.Hs.eg.db" ) {
            extra_genes <- get("extra_human_genes")
        } else if ( organism == "org.Mm.eg.db" ) {
            extra_genes <- get("extra_mouse_genes")
        }
    }
    
    .info("Extracting features")
    
    ## define genes
    genes <- rownames(counts_nm)
    if (is.null(genes) | length(genes) == 0) {
        .info("Please annotate your expression matrix with genes identifiers as rownames")
        return(NULL)
    }
    
    #GENERATE ALL FEATURES       
    features_all <- feature_generation(counts_nm, read_metrics, GO_terms, 
                                       extra_genes, organism)
    .info(paste0("Features extracted."))
    sds <- apply(features_all, 2, sd)
    #REMOVE 0-VARIANCE VALUE FEATURES
    features_all <- features_all[,sds != 0]
    types <- c("all", "common")
    ## define common features
    features_common <- features_all[, which(colnames(features_all) %in% 
                                                common_features)]
    ## write features to file if required
    if (prefix != "" && output_dir != "") {
        ## define output directory and create if needed
        o <- paste(output_dir, prefix, sep = "/")
        print(output_dir)
        dir.create(o, showWarnings = TRUE, recursive = TRUE)
        
        ## NB may be better to use file.path() for defining output files
        f_all <- file.path(o, paste0(prefix, ".", types[1], ".features"))
        f_common <- file.path(o, paste0(prefix, ".", types[2], ".features"))
        write.table(features_common, f_common)
        write.table(features_all, f_all)
        .info(paste0("Features saved: ", f_all))
        .info(paste0("Features saved: ", f_common))
    }
    ## return features in list
    return(list(features_all, features_common))
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
#' @examples 
#' data(param_mES_all)
#' data(training_mES_features)
#' data(training_mES_labels)
#' data(mES1_features)
#' data(mES1_labels)
#' mES1_features_all <- mES1_features[[1]]
#' training_mES_features_all <- training_mES_features[[1]]
#' mES1_quality_SVM <- assess_cell_quality_SVM( training_mES_features_all, 
#' training_mES_labels[,2], param_mES_all, mES1_features_all)
assess_cell_quality_SVM <- function(training_set_features, training_set_labels,
                                    ensemble_param, test_set_features) {
    
    train_f <- as.character(colnames(training_set_features))
    test_f <- as.character(colnames(test_set_features))
    
    feature_inter <- intersect(train_f, test_f)
    train_f_not_in <- which((train_f %in% feature_inter)  == TRUE)
    test_f_not_in <- which((test_f %in% feature_inter)  == FALSE)
    
    #Remove features that are not the same
    if (length(feature_inter) != length(train_f)) {
        print("Following features are not compatible and will be removed:")
      if (length(train_f_not_in) > 0) {
        print(train_f[train_f_not_in])
        training_set_features = training_set_features[-train_f_not_in]
      }
      
      if (length(test_f_not_in) > 0) {
        print(test_f[test_f_not_in])
        test_set_features = test_set_features[-test_f_not_in]
      }
    }
    
    training_set_features<-training_set_features[,order(colnames(training_set_features))]
    test_set_features<-test_set_features[,order(colnames(test_set_features))]
    
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
    final <- .vote(final_results)
    final_df <- data.frame(cell = rownames(test_set_features), quality = final)
    return(final_df)
}


################################################################################
## assess_cell_quality_PCA

#' ASSESS CELL QUALITY USING PCA AND OUTLIER DETECTION
#' 
#' @param features Input dataset containing features (cell x features)
#' @param file  Output_file where plot is saved
#' 
#' @details This function applies PCA on features and uses outlier detection to
#'  determine which cells are low and which are high quality
#' @return Returns a dataframe indicating which cell is low or high quality (0 
#' or 1 respectively)
#' 
#' @importFrom mvoutlier pcout
#' @export
#' @examples 
#' data(training_mES_features)
#' training_mES_features_all <- training_mES_features[[1]]
#' training_quality_PCA_allF <- assess_cell_quality_PCA(training_mES_features_all)
assess_cell_quality_PCA <- function(features, file="") {
    
    ## perform PCA
    pca <- prcomp(features, scale = TRUE, center = TRUE)
    pca_var_explained <- summary(pca)
    pcout_c <- mvoutlier::pcout(pca$x[, 1:2])
    dimens <- min(10, ncol(features))
    low_qual_i <- which(pcout_c$wfinal01 == 0)
    uni_2 <- (uni.plot(pca$x[, 1:2]))
    low_qual_i <- which(uni_2$outliers == TRUE)
    
    #DETERMINE WHICH OF TWO POPULATIONS ARE OUTLIERS.
    #RELY ON THAT IF MTDNA HIGH OR MAPPED PROP LOW, IT IS LOW QUALITY 
    #IF NOT AVAILABLE, ASSUME THAT CLUSTER WITH SMALLER NUMBER OF CELLS LOW QUALITY
    mtdna <- NA
    mtdna_i <- grep("mtDNA", colnames(features))
    if(length(mtdna_i) > 0) {
        mtdna <- t.test(features[,mtdna_i][low_qual_i], 
                        features[,mtdna_i][-low_qual_i], 
                        alternative = "greater")$p.value
    } 
    mapped_prop_i <- grep("Mapped", colnames(features))
    
    
    mapped_prop <- NA
    if(length(grep("Mapped", colnames(features))) > 0) {
        
        mapped_prop <- t.test(features[,mapped_prop_i][low_qual_i], 
                              features[,mapped_prop_i][-low_qual_i], 
                              alternative = "less")$p.value
    } 
    types <- rep(1, nrow(features))
    #Determine which cluster is low and which high quality
    if (!is.na(mapped_prop) && !is.na(mtdna)) {
        if (mapped_prop > 0.5 && mtdna > 0.5) {
            types[-low_qual_i] <- 0
        } else{
            types[low_qual_i] <- 0
        }
    } else {
        popul_1 <- length(low_qual_i) 
        popul_2 <- nrow(features) - length(low_qual_i)
        if (popul_1 <  popul_2) {
            types[low_qual_i] <- 0
        } else {
            types[-low_qual_i] <- 0
        }
    }
    annot <- data.frame(cell = rownames(features), quality = types)
    
    
    #PLOT PCA + MOST INFORMATIVE FEATURES
    if (file != "") {
        ## define data frame with cell types
        col <- c("0" = "red","1" = "darkgreen")
        plot_pca(features, as.character(types), pca, col, output_file = file)
    }
    return(annot)
}


################################################################################
## normalise_by_factor

#' Internal function to normalize by library size
#' 
#' @param counts matrix of counts
#' @param norm_factor vector of normalisation factors
#' @return a matrix with normalized gene counts
#' 
#' @export
#' @examples 
#' data(sample_counts)
#' data(sample_stats)
#' sample_counts_nm <- normalise_by_factor(sample_counts, colSums(sample_counts))
normalise_by_factor <- function(counts, norm_factor) { 
    return(t(t(counts) / norm_factor))
}

################################################################################

#' Helper Function to create all features
#' 
#' @param counts_nm Gene expression counts dataframe (genes x cells). Either 
#' normalised by library size or TPM values
#' @param read_metrics Dataframe with mapping statistics produced by python 
#' pipeline
#' @param GO_terms DataFrame with gene ontology term IDs, that will be used in 
#' feature extraction
#' @param extra_genes Additional genes used for feature extraction
#' @param organism The target organism to generate the features for
#' 
#' @importFrom AnnotationDbi Term
#' @importFrom topGO annFUN.org
#' @importFrom graphics hist par
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @importFrom stats cor formula mahalanobis prcomp predict qchisq quantile sd t.test var
#' 
#' @return Returns the entire set of features in a data.frame
#' 
feature_generation <- function(counts_nm, read_metrics, GO_terms, extra_genes, 
                               organism) {
    ## initialise features list  
    features <- list()
    
    read_metrics <- data.frame(read_metrics)
    
    #REMOVE ALL 0 GENES
    counts_nm <- data.frame(counts_nm)
    genes_mean <- rowMeans(counts_nm)
    genes_zero <- which(genes_mean == 0)
    
    if(length(genes_zero) > 0) {
      genes_mean <- genes_mean[-genes_zero]
      counts_nm_mean <- counts_nm[-genes_zero,] / genes_mean
    } else{
      counts_nm_mean <- counts_nm / genes_mean
    }
    ########################################
    ####TECHINCAL FEATURES##################
    ########################################
    
    ## only consider reads mapped to genes (excluding ERCCs)
    ercc_counts <- read_metrics$ercc
    if ( is.null(ercc_counts) ) {
        ercc_counts <- 0
    }
    number_mapped_reads_prop <- ((read_metrics$mapped - ercc_counts) / 
                                     read_metrics$total)
    
    #HOPE THAT THIS WORKS ALSO WHEN REGRESSION NORMALIZATION HAS BEEN APPLIED
    #AS SOME REGRESSION METHODS ALSO PUSH 0 TO ANOTHER VALUE
    detected_genes <- apply(counts_nm, 2, function(x) {
        return(length(which(x > 0)))
    })
    
    counts_nm_mean_log <- log(counts_nm_mean + 0.001)
    genes_var <- apply(counts_nm_mean_log,1,var)
    genes_means_log <- log(genes_mean)
    
    
    cell_to_mean_corr_spearman <- cor(counts_nm, rowMeans(counts_nm),
                                      method = "spearman")
    
    #SELECT ONLY HIGHLY VARIABLES GENES WITH STRONG EXPRESSION BASED ON LAST QUANTILE
    i <- which(genes_var > quantile(genes_var)[4] & 
                   genes_means_log > quantile(genes_means_log)[4])
    counts_nm_mean_log_high_var_mean <- counts_nm_mean_log[i,]
    
    transcriptome_variance <- matrix(0, ncol(counts_nm_mean_log_high_var_mean))
    number_of_highly_expressed_variable_genes <- matrix(0, ncol(counts_nm_mean_log_high_var_mean))
    num_of_high_var_exp_genes_interval <- matrix(0, ncol(counts_nm_mean_log_high_var_mean))
    if(length(i) > 100) {
        
        #VARIANCE ACROSS HIGHLY EXPRESSED GENES
        transcriptome_variance <- apply(counts_nm_mean_log_high_var_mean, 2, var)
        
        
        #NUMBER OF HIGHLY EXPRESSED GENES
        for (j in 1:ncol(counts_nm_mean)) {
          m <- mean(counts_nm_mean[i,j])
          number_of_highly_expressed_variable_genes <- 
              apply(counts_nm_mean[i,], 2, function(x) {return(sum(x > m))})
        }
        
        #NUMBER OF LOW TO HIGH EXPRESED AND VARIABLE GENES PER INTERVAL
        num_of_high_var_exp_genes_interval <-
            apply(counts_nm_mean_log_high_var_mean, 2, function(x) {
                hst <- hist(x, breaks = c(-100, -4, -2, 0, 2, 100), plot = FALSE)
                hst$counts
            })
        num_of_high_var_exp_genes_interval <- t(num_of_high_var_exp_genes_interval)
        colnames(num_of_high_var_exp_genes_interval) <- 
            paste0("num_of_high_var_exp_genes_interval_", 
                   1:ncol(num_of_high_var_exp_genes_interval))
    }
    
    
    mean_ex <- apply(counts_nm, 1, mean)
    i <- order(mean_ex, decreasing = FALSE)
    mean_ex <- mean_ex[i]
    lowl_expr <- mean_ex[1:(length(mean_ex)*0.01)]
    
    #ONLY LOWLY EXPRESSED
    l_i <- which(rownames(counts_nm) %in% names(lowl_expr))
    cell_to_mean_corr_spearman_low_ex = matrix(0, ncol(counts_nm))
    if (length(l_i) > 100) {
        cell_to_mean_corr_spearman_low_ex <- cor(counts_nm[l_i,], rowMeans(counts_nm[l_i,]),
                                                 method = "spearman")
    }
    
    techincal_features <- cbind(number_mapped_reads_prop,
                                read_metrics[, 6:11], detected_genes, cell_to_mean_corr_spearman, cell_to_mean_corr_spearman_low_ex,
                                transcriptome_variance,
                                num_of_high_var_exp_genes_interval, 
                                number_of_highly_expressed_variable_genes)
    tech_names <- c("Mapped %", "Multi-mapped %", "Intergenic %", "Intragenic %", "Exonic %", "Intronic %", "Ambigious %", 
                    "#Detected genes",  "Cell-to-mean", "Cell-to-mean lowE", "Transcriptome variance", 
                    paste0("High expr.var genes intv.", 1:ncol(num_of_high_var_exp_genes_interval)),  "#High exp + var genes")
    
    colnames(techincal_features) <- tech_names
    
    ########################################
    ####BIOLOGICAL FEATURES##################
    ########################################
    
    
    #ASSUME IT IS MOUSE
    GO_BP <- topGO::annFUN.org("BP", mapping = organism, ID = "ensembl")
    GO_CC <- topGO::annFUN.org("CC", mapping = organism, ID = "ensembl")
    
    
    GO <- c(GO_BP, GO_CC)
    
    #PROPORTION OF MAPPED READS MAPPED TO THE GO TERM
    go_prop <- sapply(unlist(GO_terms), function(go_id) {
        prop <- sum_prop(counts_nm, unlist(GO[go_id])) 
        return(prop)
    }, simplify = FALSE)
    go_prop <- do.call(cbind, go_prop)
    go_names <- Term(unlist(GO_terms))
    go_names <- sapply(go_names, simple_cap)
    colnames(go_prop) <- go_names 
    
    #CYTOPLASM AND MEMBRANE PRESENT IN GO TERMS
    m_i <- which(GO_terms[,1]  == "GO:0016020")
    c_i <- which(GO_terms[,1]  ==  "GO:0005737")
    
    volume_surface_ratio <- matrix(0, ncol(counts_nm))
    if (length(m_i) > 0 && length(c_i) > 0 && sum(go_prop[, c_i]) > 0) {
        volume_surface_ratio <- go_prop[, m_i]/go_prop[, c_i]
    }
    
    #PROPORTION OF MAPPED READS MAPPED TO SPECIFIC GENES
    extra_genes_prop <- sapply(extra_genes, function(extra_g) {
        prop <- sum_prop(counts_nm, extra_g) 
        return(prop)
    }, simplify = FALSE)
    extra_genes_prop <- do.call(cbind, extra_genes_prop)
    colnames(extra_genes_prop) <- unlist(names(extra_genes))
    
    biological_features <- cbind(go_prop, extra_genes_prop)
    colnames(biological_features) <- paste0(colnames(biological_features), " %")
    
    features <- data.frame(techincal_features, biological_features, volume_surface_ratio)
    colnames(features) <- c(colnames(techincal_features), 
                            colnames(biological_features), "Volume-surface ratio")
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
#' @param genes_interest dataframe of genes of interest to merge
#' @return a vector of sums per group
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
#' @return a character vector in title case
#' 
simple_cap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep = "", collapse = " ")
}


################################################################################
## uni.plot

#' Internal function to detect outliers from the mvoultier pacakge
#' Modified slightly so that plots are not printed
#' @param x A matrix containing counts
#' @param symb Symbols
#' @param quan quan
#' @param alpha alpha
#' @importFrom robustbase covMcd
#' @importFrom mvoutlier arw
#' @importFrom grDevices dev.off pdf rainbow
#' @import grid
#' 
#' @return a list of outlier indicators
#' 
uni.plot <- function(x, symb = FALSE, quan = 1/2, alpha = 0.025)  {
    if (!is.matrix(x) && !is.data.frame(x)) 
        stop("x must be matrix or data.frame")
    if (ncol(x) < 2) 
        stop("x must be at least two-dimensional")
    if (ncol(x) > 10) 
        stop("x should not be more than 10-dimensional")
    rob <- covMcd(x, alpha = quan)
    xarw <- mvoutlier::arw(x, rob$center, rob$cov, alpha = alpha)
    if (xarw$cn != Inf) {
        alpha <- sqrt(c(xarw$cn, qchisq(c(0.75, 0.5, 0.25), ncol(x))))
    }
    else {
        alpha <- sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(x)))
    }
    dist <- mahalanobis(x, center = rob$center, cov = rob$cov)
    sx <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for (i in 1:ncol(x)) sx[, i] <- (x[, i] - xarw$m[i])/sqrt(xarw$c[i, 
                                                                     i])
    r <- range(sx)
    if (symb == FALSE) {
        for (i in 1:ncol(x)) {
            
            o <- (sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, 
                                                              dim(x)[2]))))
            l <- list(outliers = o, md = sqrt(dist))
        }
    }
    if (symb == TRUE) {
        rd <- sqrt(dist)
        lpch <- c(3, 3, 16, 1, 1)
        lcex <- c(1.5, 1, 0.5, 1, 1.5)
        lalpha <- length(alpha)
        xs <- scale(x) - min(scale(x))
        eucl <- sqrt(apply(xs ^ 2, 1, sum))
        rbcol <- rev(rainbow(nrow(x), 
                             start = 0, end = 0.7))[
                                 as.integer(cut(eucl, nrow(x), labels = 1:nrow(x)))]
        
        o <- (sqrt(dist) > min(sqrt(xarw$cn), sqrt(qchisq(0.975, 
                                                          dim(x)[2]))))
        l <- list(outliers = o, md = sqrt(dist), euclidean = eucl)
    }
    par(yaxt = "s")
    l
}


################################################################################
## plot_pca

#' Plots PCA of all features. Colors high and low quality cells based on outlier detection. 
#' 
#' @param features Input dataset containing features (cell x features)
#' @param annot Matrix annotation of each cell
#' @param pca PCA of features
#' @param col color code indicating what color high and what low quality cells
#' @param output_file where plot is stored
#' 
#' @details This function plots PCA of all features + most informative features
#' @return Plots of PCA
#' 
#' @import ggplot2 
#' 
plot_pca <- function(features, annot, pca, col, output_file){
    
    feature_names <- colnames(features)
    ## define data frame 
    data_frame <- data.frame(type = as.character(annot), pca$x)  
    ## plot PC1 vs PC2
    plot <- ggplot2::ggplot(data_frame, ggplot2::aes_string(x = "PC1", y = "PC2")) + 
        ggplot2::geom_point(ggplot2::aes_string(colour = "type")) + 
        ggplot2::scale_colour_manual(values = col) +
        ggplot2::theme_bw()  + 
        ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(), 
                       axis.text.y = ggplot2::element_text(), axis.ticks.length = grid::unit(0, "mm"),
                       axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border = ggplot2::element_rect(fill = NA, color = "black", 
                                                            linetype = "solid"),
                       panel.grid.major = ggplot2::element_blank(), 
                       panel.grid.minor = ggplot2::element_blank(),
                       plot.background = ggplot2::element_blank()) 
    
    
    #Plot top 3 features for PC1
    top_PC1_i <- order(abs(pca$rotation[,1]), decreasing = TRUE)
    top_f_pc1 <- names(pca$rotation[(top_PC1_i[1:3]), 1])
    top_features_pc1_i <- which(feature_names %in% top_f_pc1)
    
    size <- 0.5
    text_size <- 12
    border_size <- 1
    
    plotsPC1 <- sapply(top_features_pc1_i, function(f) {
        feature <- features[,f]
        df <- data.frame(counts = log(feature + 0.0001), type = annot)
        plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "type")) + 
            ggplot2::geom_boxplot(ggplot2::aes_string(colour = "factor(type)", 
                                                      y = "counts"),
                                  alpha = 0.3, size = size, 
                                  outlier.size = 0) +
            ggplot2::ggtitle(feature_names[f]) + 
            ggplot2::theme_bw() + 
            ggplot2::theme(axis.line = ggplot2::element_blank(), 
                           axis.text.x = ggplot2::element_blank(), 
                           axis.text.y = ggplot2::element_text(size = text_size),
                           axis.ticks.length = grid::unit(0, "mm"), 
                           axis.title.x = ggplot2::element_blank(), 
                           axis.title.y = ggplot2::element_blank(),
                           legend.position = "none",
                           plot.title = ggplot2::element_text(size = text_size),
                           panel.background = ggplot2::element_blank(),
                           panel.border = ggplot2::element_rect(
                               fill = NA, color = "black", size = border_size, 
                               linetype = "solid"), panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background = ggplot2::element_blank()) + 
            ggplot2::scale_color_manual(values = col)
        return(plot)
    }, simplify = FALSE)
    
    #Plot top 3 features for PC2
    top_PC2_i <- order(abs(pca$rotation[,2]), decreasing = TRUE)
    top_f_pc2 <- names(pca$rotation[(top_PC2_i[1:3]),2])
    top_features_pc2_i <- which(feature_names %in% top_f_pc2)
    
    plotPC2 <- sapply(top_features_pc2_i, function(f) {
        feature <- features[,f]
        df <- data.frame(counts = log(feature + 0.0001), type = annot)
        plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "type")) + 
            ggplot2::geom_boxplot(ggplot2::aes_string(colour = "factor(type)", 
                                                      y = "counts"), 
                                  alpha = 0.3,  size = size, 
                                  outlier.size = 0) +
            ggplot2::ggtitle(feature_names[f]) + 
            ggplot2::theme_bw() + 
            ggplot2::theme(axis.line = ggplot2::element_blank(), 
                           axis.text.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_text(size = text_size),
                           axis.ticks.length = grid::unit(0, "mm"),
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           legend.position = "none",
                           panel.background = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(size = text_size),
                           panel.border = ggplot2::element_rect(fill = NA, 
                                                                color = "black", 
                                                                size = border_size, 
                                                                linetype = "solid"),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           plot.background = ggplot2::element_blank()) + 
            ggplot2::scale_color_manual(values = col)
        return(plot)
    }, simplify = FALSE)
    
    #ARRANGE TOP FEATURES ONTO A GRID
    #PLOT PCA IN THE MIDDLE AND FEATURES LEFT AND BOTTOM
    l <- matrix(c(2, 2, 3, 3, 4, 4, rep(1,5), 5, rep(1,5), 5, rep(1,5), 6, 
                  rep(1,5), 6,rep(1,5), 7), nrow = 6)
    l <- cbind(l, c(rep(1, 5), 7))
    
    tp <- matrix(c(2, 2, 3, 3, 4, 4, 10, rep(1, 5), 5, 5, rep(1, 5), 5, 5,
                   rep(1, 5), 6, 6, rep(1, 5), 6, 6, rep(1, 5), 7, 7), nrow = 7)
    tp <- cbind(tp, c(rep(1, 5), 7, 7))
    
    pdf(output_file, width = 10, height = 7)
    multiplot(plot,  plotlist = c(plotsPC1, plotPC2), layout = l)
    dev.off()
    
    multiplot(plot,  plotlist = c(plotsPC1, plotPC2), layout = l)
    
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
#' @return a plot object
#' 
multiplot <- function(..., plotlist = NULL, file, cols = 6, layout = NULL) {
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots <- length(plots)
    
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
#' @param predictions Predicted labels
#' @return numeric vector
#' 
.vote <- function(predictions) {
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
#' @return a string object of information (invisible)
#' 
.info <- function(string) { 
    print(paste0("[INFO]:", string))
}
