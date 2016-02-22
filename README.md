# cellity: Classification of low quality cells in scRNA-seq data using R

The `cellity` package contains functions to help to identify low quality cells in scRNA-seq data. 
It extracts biological and technical features from gene expression data that help to detect low quality cells.

Input Requirements: 
`cellity` requires 1xgene expression matrix (Genes x Cells) and 1x Read statistics matrix (Cells x Metrics) that can be obtained by processing your data first with [`celloline`]. This will map your data and generate a counts-table + statistics about reads. For further details about [`celloline`] please check (https://github.com/ti243/celloline).

The package features:

* ability to extract meaningful biological and technical features from gene expression data
* PCA-based visualisation of low quality cells and illustration of most informative features
* SVM-based classification of low quality cells

Future versions of `cellity` may also incorporate:

* Wider range of biological and technical features that could help to distinguish high from low quality cells
* Alternatives to SVM and PCA for detection and visualisation 
* More organisms and cell types that are supported by default

See below for information about installation, getting started and highlights of the package.

## Installation
The `cellity` package has  recently been submitted to
[Bioconductor](http://bioconductor.org/) and is currently under review. Once the 
package is accepted into Bioconductor, the most reliable way to install the 
package will be to use the usual Bioconductor method:

```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("scater")
```

The `cellity` package should become available in the next Bioconductor release in
April 2016. In the meantime, `cellity` can be installed directly from GitHub as
described below. We recommend using Hadley Wickham's `devtools` package to 
install `cellity` directly from GitHub. If you don't have `devtools` installed, 
then install that from CRAN (as shown below) and then run the call below to 
install `cellity`:

**If you are using the development version of R, 3.3:**
```{r}
install.packages("devtools")
devtools::install_github("ti243/cellity", build_vignettes = TRUE)
```

As the `cellity` package has recently been submitted to Bioconductor, 
development of the package is proceeding with the development version of R 
(version 3.3). As such, using `cellity` with the current release version of R is
not supported. Assuming eventual acceptance at Bioconductor, `cellity` will be 
made available through Bioconductor in the next release in April, at which point
both a "release" version (that will operate with the release version of R) and a
"devel" version (which will depend on the appropriate development version of R) 
will be available and supported.

If you would like to use `cellity` with the current version of R (3.2.3 at the 
time of writing), please leave a note in the Issues section of this repository 
and the authors will consider making this possible.

There are several other packages from CRAN and Bioconductor that `cellity` uses,
so you will need to have these packages installed as well. The CRAN packages
should install automatically when `cellity` is installed.

Not all of the following are strictly necessary, but they enhance the
functionality of `cellity` and are good packages in their own right. The commands
below should help with package installations.

CRAN packages:

```{r}
install.packages(c("ggplot2", "knitr","testthat", "mvoutlier"))
```


You might also like to install `dplyr` for convenient data manipulation:

```{r}
install.packages("dplyr")
```


## Getting started

<!---
The best place to start is the [vignette](http://htmlpreview.github.io/?http://github.com/davismcc/scater/blob/master/vignettes/vignette.html).
-->

The best place to start is the vignette. From inside an R session, load `cellity`
and then browse the vignettes:

```{r}
library(cellity)
browseVignettes("cellity")
```

There is a detailed HTML document available that introduces the main features
and functionality of `cellity`.

## `cellity` workflow

The diagram below provides an overview of the functionallity of `cellity` and its partner tool `celloline`(https://github.com/ti243/celloline).

![Diagram outlining the cellity workflow](inst/cellity_overview.png)


## Highlights

The `cellity` package allows you:

1. ability to extract meaningful biological and technical features from gene expression data
2. PCA-based visualisation of low quality cells and illustration of most informative features
3. SVM-based classification of low quality cells

For details of how to use these functions, please consult the **vignette** and **package documentation**.  The plots shown use the example data included with the package and as shown require only one or two lines of code to generate.

### Feature extraction

Use the `extract_features` to extract both, technical and biological features from your data. To use the function you need pre-processed data. You need to provide a matrix containing gene expression levels (Genes x Cells) that either raw reads normalised by total counts or TPM values (e.g. converted from FPKM cufflinks). Also, you will need a matrix with read statistics of your data (Cells x Metrics). 

To get both pieces of data, use `celloline`(https://github.com/ti243/celloline) which is a python package that can pre-process your data and outputs both, a gene expression and read statistics matrix.

### Visualisation of low quality cells `assess_cell_quality_PCA`

The `assess_cell_quality_PCA` function provides the ability to visualise low quality cells by performing a PCA on features. It than uses an outlier detection algorithm to determine what are high and low quality cells. Moreover, it outputs the most informative features that can indicate what is going on in these cells in terms of quality.

### Classification of low quality cells `assess_cell_quality_SVM`

Using `assess_cell_quality_SVM` one predict low quality cells in their data by using either our original mES data (960 cells, available within the package) that contain all, and common (applicable to most cell types) features. Or, alternatively one can use part of their own data (if prior quality annotation is available) to train the SVM and predict on the remainder of data.


The package is currently in an Beta state. The major functionality of the 
package is settled, but it is still under development so may change from time 
to time. Please do try it and contact me with bug reports, feedback, feature 
requests, questions and suggestions to improve the package.

Tomislav Ilicic & Davis McCarthy, February 2016
