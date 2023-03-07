#' deconvolution
#'
#' Get predicted cell type proportions using a deconvolution method.
#'
#' @param object A \code{\linkS4class{deconParam}}-type object.
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
#' \code{\link{deconParam}} exists.
#' 
#' @seealso
#' \linkS4class{deconParam}, \linkS4class{nnlsParam}, \linkS4class{musicParam}, 
#' \linkS4class{epicParam}
#' 
#' @author Sean Maden
#' 
#' @aliases 
#' deconvolute
#'
#' @export
setGeneric("deconvolution", function(object) standardGeneric("deconvolution"))

