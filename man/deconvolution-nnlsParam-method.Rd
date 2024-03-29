% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nnlsParam-class.R
\name{deconvolution,nnlsParam-method}
\alias{deconvolution,nnlsParam-method}
\title{Deconvolution method for nnlsParam}
\usage{
\S4method{deconvolution}{nnlsParam}(object)
}
\arguments{
\item{object}{An object of class \linkS4class{nnlsParam} (see 
\code{?nnlsParam}).}
}
\value{
Either a vector of predicted proportions, or a list containing 
predictions, metadata, and original outputs.
}
\description{
Defines the deconvolution method for \linkS4class{nnlsParam}.
}
\details{
Takes an object of class \linkS4class{nnlsParam} as input, returning 
either a list containing proportions, return info, and metadata, or a vector 
of predicted cell type proportions. 

The key term mappings for this method include:
* \code{A} : \code{bulkExpression}, bulk signals matrix (Y).
* \code{b} : \code{referenceExpression}, signature matrix (Z).
}
\examples{
exampleList <- getDeconvolutionExampleData()
param <- nnlsParam(
cellScaleFactors=exampleList[["cellScaleFactors"]], 
bulkExpression=exampleList[["bulkExpression"]],
referenceExpression=exampleList[["referenceExpression"]])

## return only predicted proportions
deconvolution(param)

# return full results
param@returnInfo <- TRUE
names(deconvolution(param))

}
\references{
Katharine M. Mullen and Ivo H. M. van Stokkum (2012). "nnls: The 
Lawson-Hanson algorithm for non-negative least squares (NNLS)." CRAN, R 
package version 1.4. URL: 
https://cran.r-project.org/web/packages/nnls/index.html
}
