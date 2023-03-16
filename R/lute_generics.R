#' deconvolution
#'
#' Get predicted cell type proportions using a deconvolution method.
#'
#' @param object A \linkS4class{deconvolutionParam}-type object.
#' @return 
#' By default, return named numeric vector of predicted proportions for each
#' cell type.
#' 
#' If \code{return.info==T}, instead returns a list including proportions, 
#' results object returned from specified method, and additional metadata.
#' 
#' @details 
#' This generic maps standard deconvolution inputs to the parameters of the
#' specified deconvolution method for which a subclass of type 
#' \code{\link{deconParam}} exists. This generic uses a similar approach to
#' the \code{bluster} R/Bioconductor package.
#' 
#' @seealso
#' \linkS4class{deconvolutionParam}, \linkS4class{referencebasedParam}, 
#' \linkS4class{independentbulkParam}, \linkS4class{nnlsParam}, 
#' \linkS4class{musicParam}, \linkS4class{epicParam}, 
#' \linkS4class{deconrnaseqParam}, \linkS4class{scdcParam},
#' \linkS4class{bisqueParam}, \linkS4class{music2Param}
#' 
#' @author Sean Maden
#' 
#' @references 
#' 
#' Aaron Lun. bluster: Clustering Algorithms for Bioconductor. (2022) 
#' Bioconductor, R package version 1.6.0.
#' 
#' @aliases 
#' deconvolute, Deconvolution, Deconvolute
#'
#' @export
setGeneric("deconvolution", function(object) standardGeneric("deconvolution"))

