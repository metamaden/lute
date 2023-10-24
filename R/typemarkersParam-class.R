#!/usr/bin/env R

### Author: Sean Maden

#' typemarkersParam-class
#'
#' Main constructor for class to manage mappings to the \code{typemarkers()} 
#' generic.
#' 
#' @include lute_generics.R
#'  
#' @details Main constructor for class \linkS4class{typemarkersParam}.
#' 
#' @rdname typemarkersParam-class
#' 
#' @seealso meanratiosParam
#' 
#' @param markersPerType Number of top markers to get per cell type.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#'
#' @returns New object.
#' 
#' @aliases 
#' TypemarkersParam-class, TypeMarkersParam-class
#'
setClass("typemarkersParam", slots=c(markersPerType="numeric", 
                                     returnInfo="logical"))

#' Make new object of class typemarkersParam
#'
#' Main constructor for class \linkS4class{typemarkersParam}.
#'
#' @param markersPerType Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param returnInfo Whether to return metadata and original marker selection 
#' method outputs with predicted proportions.
#'
#' @returns New object of class \linkS4class{typemarkersParam}.
#'
#' @details This is the main parent class for cell type gene marker 
#' identification methods. Currently supported methods and their child classes include:
#' 
#' 1. Mean Ratios: The method DeconvoBuddies::get_mean_ratios2(), supported by the
#' class meanratiosParam.
#' 
#' @examples
#' example.data <- getDeconvolutionExampleData()
#' 
#' @export
typemarkersParam <- function(markersPerType=20, returnInfo=FALSE) {
  new("typemarkersParam", markersPerType=markersPerType, 
      returnInfo=returnInfo)
}

#' Method for class \linkS4class{typemarkersParam}
#' 
#' @param object An object of class \linkS4class{typemarkersParam}.
#' 
#' @returns Info related to gene markers for cell types.
#'
#' @examples
#' example.data <- getDeconvolutionExampleData()
#'
#' @export
setMethod("typemarkers", signature(object="typemarkersParam"), function(object){
  parametersList <- callNextMethod()
  ## instantiate and format objects
  markersPerType <- parametersList[["markersPerType"]]
  returnInfo <- parametersList[["returnInfo"]]
})

#' Inspect slot in \linkS4class{typemarkersParam} object
#' @param x Object to access.
#' @param i Slot to access.
#' @returns Contents of specified slot.
#' @details Inspect slot in \linkS4class{typemarkersParam} object
#' 
#' @examples
#' example.data <- getDeconvolutionExampleData()
#' 
#' @export
setMethod("[[", "typemarkersParam", function(x, i) {slot(x, i)})

#' Show generic behavior for object of class \linkS4class{typemarkersParam}
#' @param object An object of class \linkS4class{typemarkersParam} (see 
#' \code{?typemarkersParam}).
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{typemarkersParam}
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' 
#' @returns Shows object summaries.
#' 
#' @export
setMethod("show", signature(object="typemarkersParam"), function(object){
  show(object)
})