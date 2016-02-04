## ----knitr-options, echo=FALSE, warning=FALSE----------------------------
## To render an HTML version that works nicely with github and web pages, do:
## rmarkdown::render("vignettes/vignette.Rmd", "all")
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')
#knitr::opts_chunk$set(echo=FALSE, fig.path='cellity/plot-', cache=TRUE)
library(ggplot2)
theme_set(theme_bw(12))

## ---- eval=TRUE, message=FALSE, warning=FALSE----------------------------
library(cellity)
data(sample_counts)
data(sample_stats)

## ---- eval=TRUE, message=FALSE, warning=FALSE, results='hide', error=FALSE----
sample_counts_nm <- normalise_by_factor(sample_counts, colSums(sample_counts))

## ---- eval=TRUE, message=FALSE, warning=FALSE, results='hide', error=FALSE----
sample_features <- extract_features(sample_counts_nm, sample_stats)

## ---- eval=TRUE, message=FALSE, warning=FALSE, error=FALSE---------------
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  #MAKE SURE YOU HAVE THE APPROPRIATE ORGANISM INSTALLED 
  #You would instsall the library below by: 
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("org.Hs.eg.db")
  library(org.Hs.eg.db)
  data("extra_human_genes")
  data("feature_info")
  GO_terms <- feature_info[[1]]
  common_features <- feature_info[[2]]
  features_human <- extract_features(
    sample_counts_nm, sample_stats, common_features = common_features, 
    GO_terms = GO_terms, extra_genes = extra_human_genes, organism = "human")
}

## ---- eval=TRUE----------------------------------------------------------
sample_features_all <- sample_features[[1]]
sample_qual_pca <- assess_cell_quality_PCA(sample_features_all)

## ---- eval=TRUE----------------------------------------------------------
data(training_mES_features)
training_mES_features_all <- training_mES_features[[1]]
training_quality_PCA_allF <- assess_cell_quality_PCA(
  training_mES_features_all, file = "~/training_quality_PCA_allF.pdf")

## ---- eval=TRUE, warning=FALSE, message=FALSE----------------------------
if (requireNamespace("caret", quietly = TRUE)) {
  library(caret)
  data(training_mES_labels)
  lvs <- c("0", "1")
  truth <- factor(training_mES_labels[,2],levels = rev(lvs))
  pred <- factor(training_quality_PCA_allF[,2], levels = rev(lvs))
  confusionMatrix(pred, truth)
}

## ---- eval=TRUE----------------------------------------------------------
training_mES_features_common <- training_mES_features[[2]]
training_quality_PCA_commonF <- assess_cell_quality_PCA(
  training_mES_features_common, file = "~/training_quality_PCA_commonF.pdf")

## ---- eval=TRUE, , warning=FALSE, message=FALSE--------------------------
if (requireNamespace("caret", quietly = TRUE)) {
  pred <- factor(training_quality_PCA_commonF[,2], levels = rev(lvs))
  confusionMatrix(pred, truth)
}

## ---- eval=TRUE----------------------------------------------------------
data(mES1_features)
data(mES1_labels)

## ---- eval=TRUE----------------------------------------------------------
data(param_mES_all)
mES1_features_all <- mES1_features[[1]]
mES1_quality_SVM <- assess_cell_quality_SVM(
  training_mES_features_all, training_mES_labels[,2], param_mES_all, 
  mES1_features_all)

## ---- eval=TRUE----------------------------------------------------------
if (requireNamespace("caret", quietly = TRUE)) {
  truth <- factor(mES1_labels[,2],levels = rev(lvs))
  pred <- factor(mES1_quality_SVM[,2], levels = rev(lvs))
  confusionMatrix(pred, truth)
}

## ---- eval=TRUE----------------------------------------------------------
data(param_mES_common)
training_mES_features_common <- training_mES_features[[2]]
mES1_features_common <- mES1_features[[2]]
mES1_quality_SVM_common <- assess_cell_quality_SVM(
  training_mES_features_common, training_mES_labels[,2], param_mES_common,
  mES1_features_common)

## ---- eval=TRUE----------------------------------------------------------
if (requireNamespace("caret", quietly = TRUE)) {
  truth <- factor(mES1_labels[,2],levels = rev(lvs))
  pred <- factor(mES1_quality_SVM_common[,2], levels = rev(lvs))
  confusionMatrix(pred, truth)
}

## ---- eval=TRUE----------------------------------------------------------
#PCA QUALITY
mES1_quality_PCA<-assess_cell_quality_PCA(mES1_features_all)

if (requireNamespace("caret", quietly = TRUE)) {
  truth <- factor(mES1_labels[,2],levels = rev(lvs))
  pred <- factor(mES1_quality_PCA[,2], levels = rev(lvs))
  confusionMatrix(pred, truth)
}

