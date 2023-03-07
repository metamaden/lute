% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bisqueParam-class.R
\docType{class}
\name{bisqueParam-class}
\alias{bisqueParam-class}
\alias{BisqueParam-class}
\title{bisqueParam-class}
\description{
Applies the BisqueRNA::ReferenceBasedDecomposition() implementation of the 
Bisque deconvolution algorithm.
}
\details{
Main constructor for class \linkS4class{bisqueParam}.
}
\examples{
# example
lexample <- .get_decon_example_data_bisque()
param <- bisqueParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])

# return only predicted proportions
deconvolution(param)

# return full results
param@return.info <- T
names(deconvolution(param))
# [1] "predictions" "result.info" "metadata"

}
\seealso{
\linkS4class{deconParam}
}