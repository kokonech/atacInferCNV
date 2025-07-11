% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareAnalysis.R
\name{prepareAtacInferCnvInput}
\alias{prepareAtacInferCnvInput}
\title{Prepare analysis for the CNV calling from ATAC data}
\usage{
prepareAtacInferCnvInput(
  dataPath = "",
  annPath = "",
  resDir = "",
  inObj = NULL,
  sId = "sample",
  targColumn = "CellType",
  ctrlGrp = "Normal",
  ctrlObj = NULL,
  binSize = NULL,
  chromLength = NULL,
  metaCells = FALSE
)
}
\arguments{
\item{dataPath}{Path to the input data in 10X format}

\item{annPath}{Path to annotation of the cells in the target column}

\item{resDir}{Path to the result directory}

\item{inObj}{Precomputed Seurat/Signac object with required input data}

\item{sId}{Result name. Default: "Sample"}

\item{targColumn}{Name of the target column in annotation. Default: "CellType"}

\item{ctrlGrp}{Name for the reference control cell type. Could be several names, separated by comma. Default: "Normal"}

\item{ctrlObj}{Seurat/Signac object to use as non-tumor control. Default: NULL}

\item{binSize}{Apply custom bin size to combine signals in windows for CNV calling
e.g. 500000 for 500 KBp. Default: NULL (not use this option)}

\item{chromLength}{Numeric vector of chromosome sizes, specific for genome. Default: NULL}

\item{metaCells}{Set TRUE to use meta cells (n=5 cells will be used to merge) or assign a number of cells. Default: FALSE}
}
\description{
Prepare analysis for the CNV calling from ATAC data
}
