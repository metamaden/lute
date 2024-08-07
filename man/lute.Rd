% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lute_framework.R
\name{lute}
\alias{lute}
\title{lute framework}
\usage{
lute(
  singleCellExperiment = NULL,
  referenceExpression = NULL,
  bulkExpression = NULL,
  bulkSummarizedExperiment = NULL,
  cellScaleFactors = NULL,
  returnInfo = FALSE,
  markersPerType = 20,
  assayName = "counts",
  cellTypeVariable = "celltype",
  typemarkerAlgorithm = NULL,
  deconvolutionAlgorithm = "nnls",
  verbose = TRUE
)
}
\arguments{
\item{singleCellExperiment}{Object of type SingleCellExperiment. Optional (see argument z).}

\item{referenceExpression}{Signature matrix of cell type-specific signals. Optional (see 
argument singleCellExperiment).}

\item{bulkExpression}{Bulk mixed signals matrix of samples, which can be matched to 
single-cell samples. Optional (see argument y.se).}

\item{bulkSummarizedExperiment}{SummarizedExperiment or similar data type containing the bulk signals
matrix in its assays (e.g. accessible with assays(y.se)[[assayName]] using
the provided assayName argument). Optional (see argument y).}

\item{cellScaleFactors}{Cell size factor transformations of length equal to the K cell types 
to deconvolve. Optional, if not provided, uses equal weights for types.}

\item{returnInfo}{Whether to return metadata and original method outputs 
with predicted proportions.}

\item{markersPerType}{Number of top markers to get per cell type.}

\item{assayName}{Name of expression matrix in singleCellExperiment, and optionally y.se, 
assays. Optional (e.g. "counts"; see arguments singleCellExperiment, y.se).}

\item{cellTypeVariable}{Name of cell type variable in singleCellExperiment coldata.}

\item{typemarkerAlgorithm}{Which type-specific marker selection algorithm to 
use. If NULL, skips type marker analyses.}

\item{deconvolutionAlgorithm}{Where deconvolution algorithm to use. If NULL, 
skips deconvolution.}

\item{verbose}{Whether to show verbose status messages.}
}
\value{
A list containing results returned from type marker selection and
deconvolution runs, with additional information returned if
 \code{returnInfo == TRUE}.
}
\description{
Obtain cell type markers and proportion predictions from various algorithms. 
Allows flexible data types and standard application of cell size scale 
factors.
}
\details{
Main function to use the lute deconvolution framework. Manages data
conversions and mappings to deconvolution experiment steps, including setup,
gene marker identification, and main deconvolution runs.

Support is provided for \linkS4class{SummarizedExperiment}-type or 
matrix-type inputs for the Z signature matrix (see referenceExpression 
argument) and Y bulk signals matrix (see bulkExpression arguments). Note,
both Z and Y need to be provided or derivable in order to run deconvolution.

Key framework steps are controlled with the arguments: 1. typemarkerAlgorithm 
and 2. deconvolutionAlgorithm, where the first argument is parsed prior to 
the second. Either step may be skipped by setting it to "NULL".
}
\examples{
# get example bulk data
bulkExpression <- getDeconvolutionExampleData()$reference

# get example singleCellExperiment
singleCellExperiment <- randomSingleCellExperiment()[seq(10),]

# get framework results
experiment.results <- lute(
singleCellExperiment=singleCellExperiment, 
bulkExpression=bulkExpression, typemarkerAlgorithm=NULL
)

}
