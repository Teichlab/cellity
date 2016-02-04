#' @name feature_info
#' @title Information which genes and GO categories should be included as features. Also defines which features are cell-type independent (common features)
#' @description This list contains metadata information 
#' that is used to extract features from in the function extract_features
#' @return NULL, but makes available a list with  metadata
#' @docType data
#' @usage feature_info
#' @format a list with 2 elements (GO_terms,common_features).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name extra_human_genes
#' @title Additional human genes that are used in feature extraction
#' @description This list contains human genes that are used for feature extraction of biological features
#' @return NULL, but makes available a list with  metadata
#' @docType data
#' @usage extra_human_genes
#' @format a list containing vectors of genes. Name indicates which GO category.
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name extra_mouse_genes
#' @title Additional mouse genes that are used in feature extraction
#' @description This list contains mouse genes that are used for feature extraction of biological features
#' @return NULL, but makes available a list with  metadata
#' @docType data
#' @usage extra_mouse_genes
#' @format a list containing vectors of genes. Name indicates which GO category.
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL


#' @name mES1_features
#' @title Real test dataset containing all and common features from the paper (mES1)
#' @description This list contains 2 dataframes
#' where each contains features per cell (cell X features) that can be used for classification.
#' @return NULL, but makes available a list with 2 dataframes
#' @docType data
#' @usage mES1_features
#' @format a list with 2 elements (all_features, common_features).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name mES1_labels
#' @title Real test dataset containing annotation of cells
#' @description This data frame has 2 columns:
#' First showing cell names, the second indicating if cell is of low (0) or high (1) quality
#' @return NULL, but makes available a dataframe with cell annotations
#' @docType data
#' @usage mES1_labels
#' @format a dataframe with 2 columns (cell_names, label).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL


#' @name param_mES_all
#' @title Parameters used for SVM classification
#' @description This data frame has 3 columns:
#' gamma, cost, class.weights  and is optimised for all features and our training data
#' @return NULL, but makes available a dataframe with parameters
#' @docType data
#' @usage param_mES_all
#' @format a dataframe with 3 columns (gamma, cost, class.weights).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name param_mES_common
#' @title Parameters used for SVM classification
#' @description This data frame has 3 columns:
#' gamma, cost, class.weights and is optimised for common features and our training data
#' @return NULL, but makes available a dataframe with parameters
#' @docType data
#' @usage param_mES_common
#' @format a dataframe with 3 columns (gamma, cost, class.weights).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL


#' @name sample_counts
#' @title Sample gene expression data containing 40 cells 
#' @description This data frame contains genes (rows)
#' and cells (columns) showing raw read counts
#' @return NULL, but makes available a dataframe with raw read counts
#' @docType data
#' @usage sample_counts
#' @format a dataframe with genes x cells
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name sample_stats
#' @title Sample read statistics data containing 40 cells 
#' @description This data frame contains read metrics (columns)
#' and cells (rows)
#' @return NULL, but makes available a dataframe with read statistics
#' @docType data
#' @usage sample_stats
#' @format a dataframe with cells x metrics
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name training_mES_features
#' @title Original training dataset containing all and common features from the paper (training mES)
#' @description This list contains 2 dataframes
#' where each contains features per cell (cell X features) that can be used for classification.
#' @return NULL, but makes available a list with 2 dataframes
#' @docType data
#' @usage training_mES_features
#' @format a list with 2 elements (all_features, common_features).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL

#' @name training_mES_labels
#' @title Original training dataset containing annotation of cells
#' @description This data frame has 2 columns:
#' First showing cell names, the second indicating if cell is of low (0) or high (1) quality
#' @return NULL, but makes available a dataframe with cell annotations
#' @docType data
#' @usage training_mES_labels
#' @format a dataframe with 2 columns (cell_names, label).
#' @source Wellcome Trust Sanger Institute
#' @author Tomislav Ilicic & Davis McCarthy, 2015-03-05
NULL
