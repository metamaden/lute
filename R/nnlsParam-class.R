#!/usr/bin/env R

### Author: Sean Maden

#' nnlsParam-class
#'
#' Uses nnls::nnls().
#' 
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#' 
#' @details Main constructor for class \linkS4class{nnlsParam}.
#' @rdname nnlsParam-class
#' @seealso \linkS4class{deconParam}
#' 
#' @examples 
#' lexample <- get_decon_example_data()
#' param <- nnlsParam(s=lexample[["s"]], y=lexample[["y"]], 
#'                     z=lexample[["z"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@return.info <- TRUE
#' names(deconvolution(param))
#'
#' @returns New object.
#' 
#' @aliases 
#' NNLSParam-class
#' 
setClass("nnlsParam", contains="referencebasedParam")

#' Make new object of class nnlsParam
#'
#' Main constructor for class \linkS4class{nnlsParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, 
#' can be computed from a provided \linkS4class{ExpressionSet} containing 
#' single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell 
#' types to deconvolve.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples 
#' lexample <- get_decon_example_data()
#' param <- nnlsParam(s=lexample[["s"]], y=lexample[["y"]], 
#'                     z=lexample[["z"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@return.info <- TRUE
#' names(deconvolution(param))
#' 
#' @returns Object of class \linkS4class{nnlsParam}
#' 
#' @seealso \linkS4class{referencebasedParam}, \linkS4class{deconvolutionParam}
#'
#' @details Main parameter class for mapping inputs to the non-negative least 
#' squares (NNLS) deconvolution algorithm, implemented as \code{nnls::nnls()}.
#' 
#' @export
nnlsParam <- function(y, z, s, return.info=FALSE) {
  new("nnlsParam", s=s, y=y, z=z, return.info=return.info)
}

#' Deconvolution method for nnlsParam
#'
#' Defines the deconvolution method for \linkS4class{nnlsParam}.
#'
#' @param object An object of class \linkS4class{nnlsParam} (see 
#' \code{?nnlsParam}).
#'
#' @details Takes an object of class \linkS4class{nnlsParam} as input, returning 
#' either a list containing proportions, return info, and metadata, or a vector 
#' of predicted cell type proportions. 
#' 
#' The key term mappings for this method include:
#' * \code{A} : \code{y}, bulk signals matrix.
#' * \code{b} : \code{z}, signature matrix.
#' 
#' @examples 
#' lexample <- get_decon_example_data()
#' param <- nnlsParam(s=lexample[["s"]], y=lexample[["y"]], 
#'                     z=lexample[["z"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@return.info <- TRUE
#' names(deconvolution(param))
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#' 
#' @references 
#' 
#' Katharine M. Mullen and Ivo H. M. van Stokkum (2012). "nnls: The Lawson-Hanson 
#' algorithm for non-negative least squares (NNLS)." CRAN, R package version 1.4. 
#' URL: https://cran.r-project.org/web/packages/nnls/index.html
#'
#' @export
setMethod("deconvolution", signature(object="nnlsParam"), function(object){
  lparam <- callNextMethod()
  input_y <- lparam[["y"]]; input_z <- lparam[["z"]]; input_s <- lparam[["s"]]
  bulk.samples.index.vector <- seq(ncol(input_y))
  result <- lapply(
    bulk.samples.index.vector, function(index){
      
      nnls::nnls(A=input_z, b=input_y[,index])
  
      }
    )
  
  names(result) <- colnames(input_y)
  predictions <- lapply(result, function(iter){iter$x})
  return_list <- .parse_deconvolution_predictions_results(predictions, 
                                                 colnames(input_z), 
                                                 colnames(input_y))
  if(object[["return.info"]]){
    
    return_list <- list(predictions=predictions, 
               result.info=result, 
               metadata=lparam[["metadata"]])
    
    }
  return(lr)
})

#' Show generic behavior for object of class \linkS4class{nnlsParam}
#' @param object An object of class \linkS4class{nnlsParam} (see 
#' \code{?nnlsParam}).
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{nnlsParam}
#' 
#' @examples
#' lexample <- get_decon_example_data()
#' param <- nnlsParam(s=lexample[["s"]], y=lexample[["y"]], z=lexample[["z"]])
#' param
#' 
#' @returns Shows object summaries.
#' 
#' @export
setMethod("show", signature(object="nnlsParam"), function(object){
  show(object)
})
