#!/usr/bin/env R

# Author: Sean Maden

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
#' \linkS4class{deconvolutionParam} exists. This generic uses a similar approach to
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
#' @importFrom methods callNextMethod is new slot
#' @importFrom stats rnbinom rpois
#'
#' @export
setGeneric("deconvolution", function(object) standardGeneric("deconvolution"))

#' typemarkers
#'
#' Get cell type gene markers using standard accessors to supported functions.
#'
#' @param object A \linkS4class{tyepmarkersParam}-type object.
#' @return 
#' By default, return a vector of marker genes.
#' 
#' If \code{return.info==T}, provides detailed results, including original 
#' outputs.
#' 
#' @details 
#' This generic manages tasks for marker gene identification. In particular, it
#' takes a specified amount of marker genes to return per type.
#' 
#' @seealso
#' \linkS4class{tyepmarkersParam}
#' 
#' @author Sean Maden
#' 
#' @aliases 
#' Typemarkers, TypeMarkers
#'
#' @export
setGeneric("typemarkers", function(object) standardGeneric("typemarkers"))
