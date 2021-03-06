% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_low_qual_cells.R
\name{assess_cell_quality_SVM}
\alias{assess_cell_quality_SVM}
\title{Assess quality of a cell - SVM version}
\usage{
assess_cell_quality_SVM(training_set_features, training_set_labels,
  ensemble_param, test_set_features)
}
\arguments{
\item{training_set_features}{A training set containing features (cells x 
features) for prediction}

\item{training_set_labels}{Annotation of each individual cell if high or low quality (1 or 0 respectively)}

\item{ensemble_param}{Dataframe of parameters for SVM}

\item{test_set_features}{Dataset to predict containing features (cells x 
features)}
}
\value{
Returns a dataframe indicating which cell is low or high quality (0 
or 1 respectively)

data.frame with decision on quality of cells
}
\description{
Assess quality of a cell - SVM version
}
\details{
This function takes a traning set + annotation to predict a test 
set. It requires that hyper-parameters have been optimised.
}
\examples{
data(param_mES_all)
data(training_mES_features)
data(training_mES_labels)
data(mES1_features)
data(mES1_labels)
mES1_features_all <- mES1_features[[1]]
training_mES_features_all <- training_mES_features[[1]]
mES1_quality_SVM <- assess_cell_quality_SVM( training_mES_features_all, 
training_mES_labels[,2], param_mES_all, mES1_features_all)
}

