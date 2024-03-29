% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bisqueParam-class.R
\docType{class}
\name{bisqueParam-class}
\alias{bisqueParam-class}
\alias{BisqueParam-class}
\title{bisqueParam-class}
\value{
New object of class \linkS4class{bisqueParam}.
}
\description{
Applies the BisqueRNA::ReferenceBasedDecomposition() implementation of the 
Bisque deconvolution algorithm.
}
\details{
Main constructor for class \linkS4class{bisqueParam}.
}
\examples{
## get data
exampleList <- getDeconvolutionExampleDataBisque()
bulkExpressionSet <- exampleList[["bulkExpressionSet"]][,seq(10)]
bulkExpression <- exprs(exampleList[["bulkExpressionSet"]])
bulkExpression <- bulkExpression[,c(11:ncol(bulkExpression))]

## get param object
newBisqueParameter <- bisqueParam(bulkExpressionSet=bulkExpressionSet, 
                     bulkExpressionIndependent=bulkExpression,
                     scData=exampleList[["singleCellExpressionSet"]], 
                     batchVariable="SubjectName", 
                     cellTypeVariable="cellType", 
                     useOverlap=FALSE)

## get predicted proportions
res <- deconvolution(newBisqueParameter)

}
\references{
Brandon Jew and Marcus Alvarez (2021). BisqueRNA: Decomposition of Bulk 
Expression with Single-Cell Sequencing. CRAN, R package version 1.0.5.
URL: https://CRAN.R-project.org/package=BisqueRNA

Brandon Jew et al. Accurate estimation of cell composition in bulk 
expression through robust integration of single-cell information. 
Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
}
\seealso{
\linkS4class{deconvolutionParam}, 
\linkS4class{referencebasedParam}, 
\linkS4class{independentbulkParam}
}
