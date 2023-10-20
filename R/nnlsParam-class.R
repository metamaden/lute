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
#' exampleList <- getDeconvolutionExampleData()
#' param <- nnlsParam(cellScaleFactors=exampleList[["cellScaleFactors"]], 
#' bulkExpression=exampleList[["bulkExpression"]],
#' referenceExpression=exampleList[["referenceExpression"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@returnInfo <- TRUE
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
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples.
#' @param referenceExpression Signature matrix of cell type-specific signals. If not provided, 
#' can be computed from a provided \linkS4class{ExpressionSet} containing 
#' single-cell data.
#' @param cellScaleFactors Cell size factor transformations of length equal to the K cell 
#' types to deconvolve.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' param <- nnlsParam(cellScaleFactors=exampleList[["cellScaleFactors"]], 
#' bulkExpression=exampleList[["bulkExpression"]],
#' referenceExpression=exampleList[["referenceExpression"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@returnInfo <- TRUE
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
nnlsParam <- function(
    bulkExpression, referenceExpression, cellScaleFactors, returnInfo=FALSE) {
  new("nnlsParam", cellScaleFactors=cellScaleFactors, 
      bulkExpression=bulkExpression, referenceExpression=referenceExpression, 
      returnInfo=returnInfo)
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
#' * \code{A} : \code{bulkExpression}, bulk signals matrix (Y).
#' * \code{b} : \code{referenceExpression}, signature matrix (Z).
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' param <- nnlsParam(
#' cellScaleFactors=exampleList[["cellScaleFactors"]], 
#' bulkExpression=exampleList[["bulkExpression"]],
#' referenceExpression=exampleList[["referenceExpression"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@returnInfo <- TRUE
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
  parametersList <- callNextMethod()
  bulkExpression <- parametersList[["bulkExpression"]]
  referenceExpression <- parametersList[["referenceExpression"]]
  cellScaleFactors <- parametersList[["cellScaleFactors"]]
  bulkSamplesIndexVector <- seq(ncol(bulkExpression))
  result <- lapply(
    bulkSamplesIndexVector, function(index){
      
      nnls::nnls(A=referenceExpression, b=bulkExpression[,index])
  
      }
    )
  
  names(result) <- colnames(bulkExpression)
  predictions <- lapply(result, function(iter){iter$x})
  returnList <- parseDeconvolutionPredictionsResults(
    predictions, colnames(referenceExpression), colnames(bulkExpression))
  if(object[["returnInfo"]]){
    
    returnList <- list(
      predictions=predictions, 
      result.info=result, 
      metadata=parametersList[["metadata"]]
    )
    
  }
  return(returnList)
})

#' Show generic behavior for object of class \linkS4class{nnlsParam}
#' 
#' @param object An object of class \linkS4class{nnlsParam} (see 
#' \code{?nnlsParam}).
#' 
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{nnlsParam}
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' param <- nnlsParam(
#' cellScaleFactors=exampleList[["cellScaleFactors"]], 
#' bulkExpression=exampleList[["bulkExpression"]],
#' referenceExpression=exampleList[["referenceExpression"]])
#' 
#' ## return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@returnInfo <- TRUE
#' names(deconvolution(param))
#' 
#' @returns Shows object summaries.
#' 
#' @export
setMethod("show", signature(object="nnlsParam"), function(object){
  show(object)
})
