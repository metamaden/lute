% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolutionParam-class.R
\docType{class}
\name{deconvolutionParam-class}
\alias{deconvolutionParam-class}
\alias{DeconvolutionParam-class,}
\alias{DeconParam-class,}
\alias{deconParam-class}
\title{deconvolutionParam-class}
\value{
New deconvolutionParam object.
}
\description{
Defines the principal parent class for all deconvolution method parameters.
}
\details{
Defines the parent class for deconvolution method parameters. Since all
deconvolution runs require a \code{y}  signals matrix, whether from 
experiment data or simulations such as pseudobulking, this parent class 
manages the bulk signals matrix. For this class, the 
\code{deconvolution} generic performs basic summaries of the bulk 
signals matrix.
}
\examples{
param <- new("deconvolutionParam")
deconvolution(param)

}
\seealso{
\code{deconvolution}
}
