#!/usr/bin/env R

### Author: Sean Maden

#' independentbulkParam-class
#'
#' Class and methods for managing methods requiring independent bulk samples.
#'
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#'
#' @param bulkExpressionIndependent Bulk mixed signals matrix of independent 
#' samples, which should not overlap samples in y.
#'
#' @details The main purpose of this class is to compare bulk sample data 
#' between the passed objects y and yi. Since we assume yi contains the 
#' independent bulk samples, it should not have overlapping sample IDs 
#' (colnames), and it should have overlapping marker IDs (rownames) compared to 
#' the reference bulk samples y.
#'
#' @seealso \linkS4class{deconParam}, \linkS4class{referencebasedParam}
#' 
#' @examples 
#' new("independentbulkParam")
#' 
#' @returns New object.
setClass("independentbulkParam", contains="referencebasedParam", 
         slots=c(bulkExpressionIndependent="matrix"))

#' Make a new \linkS4class{independentbulkParam} object
#' 
#' Function to make a new object of class \linkS4class{independentbulkParam}
#'
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param bulkExpressionIndependent Bulk mixed signals matrix of independent 
#' samples, which should not overlap samples in y.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a provided ExpressionSet containing 
#' single-cell data.
#' @param cellScaleFactors Cell size scale factor transformations of 
#' length equal to the K cell types to deconvolve.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples 
#' new("independentbulkParam")
#'
#' @returns New object.
#'
#' @export
independentbulkParam <- function(
    bulkExpression=NULL, bulkExpressionIndependent=NULL, 
    referenceExpression=NULL, cellScaleFactors=NULL, returnInfo=FALSE) {
    if(is(bulkExpression, "NULL")){
      bulkExpression <- matrix(0)}
    if(is(referenceExpression, "NULL")){
      referenceExpression <- matrix(0)}
    if(is(bulkExpressionIndependent, "NULL")){
      bulkExpressionIndependent <- matrix(0)}
    if(is(cellScaleFactors, "NULL")){
      cellScaleFactors <- rep(1, ncol(referenceExpression))}
    param <- new("independentbulkParam", 
                 bulkExpression=bulkExpression, 
                 bulkExpressionIndependent=bulkExpressionIndependent, 
                 referenceExpression=referenceExpression, 
                 cellScaleFactors=cellScaleFactors, 
                 returnInfo=returnInfo)
    return(param)
}

#' Deconvolution method for class \linkS4class{independentbulkParam}
#'
#' Function to perform standard operations prior to deconvolution (a.k.a. 
#' "deconvolution prep") for an object of class 
#' \linkS4class{independentbulkParam}.
#'
#' @param object An object of class \linkS4class{independentbulkParam}.
#'
#' @details Takes an object of \linkS4class{independentbulkParam} class as 
#' input, and returns a list with the filtered/checked/parsed experiment 
#' objects.
#' 
#' @examples 
#' new("independentbulkParam")
#'
#' @returns Method results.
#'
#' @export
setMethod("deconvolution", "independentbulkParam", function(object) {
    parameterList <- callNextMethod()
    uniqueMarkerLabels <- uniqueSampleLabels <- NULL
    overlappingMarkerLabels <- overlappingSampleLabels <- NULL
    bulkExpression <- parameterList[["bulkExpression"]] # get bulk data
    bulkExpressionIndependent <- parameterList[["bulkExpressionIndependent"]]
    markersBulkExpression <- rownames(bulkExpression) # parse bulk marker IDs
    markersBulkExpressionIndependent <- rownames(bulkExpressionIndependent) 
    
    ## compare marker labels and subset yi on overlapping markers
    if(is(markersBulkExpression, "NULL")|
       is(markersBulkExpressionIndependent, "NULL")){
        message(
          paste0("Warning, no marker labels found in either Y ",
                 "(bulkExpression) or Yi (bulkExpressionIndependent)."))
    } else{
        uniqueMarkerLabels <- 
          unique(markersBulkExpression, markersBulkExpressionIndependent)
        overlappingMarkerLabels <- 
          intersect(markersBulkExpression, markersBulkExpressionIndependent)
        if(length(overlappingMarkerLabels) > 0){
          bulkExpressionIndependent <- 
            bulkExpressionIndependent[overlappingMarkerLabels,]
        }
    }
  
    ## compare sample labels and remove overlapping samples
    samplesBulkExpression <- colnames(bulkExpression) # parse bulk sample IDs
    samplesBulkExpressionIndependent <- colnames(bulkExpressionIndependent) 
    ## compare sample IDs
    if(is(samplesBulkExpression, "NULL")|
       is(samplesBulkExpressionIndependent, "NULL")){
        message(
          paste0("Warning, no sample labels found in either ",
                 "samplesBulkExpression or samplesBulkExpressionIndependent."))
    } else{
        uniqueSampleLabels <- unique(
          samplesBulkExpression, samplesBulkExpressionIndependent)
        overlappingSampleLabels <- intersect(
          samplesBulkExpression, samplesBulkExpressionIndependent)
        if(length(overlapping.samples) > 0){
            bulkFilter <- 
              !colnames(bulkExpressionIndependent) %in% overlappingSampleLabels
            bulkExpressionIndependent <- 
              bulkExpressionIndependent[, bulkFilter, drop=FALSE]
        }
    }
  
    ## parse return list
    ## get metadata to return
    metadataList <- list(
      uniqueMarkerLabels=uniqueMarkerLabels,
                uniqueSampleLabels=uniqueSampleLabels,
                overlappingMarkerLabels=overlappingMarkerLabels,
                overlappingSampleLabels=overlappingSampleLabels
    )
    returnList <- list(
      bulkExpression=bulkExpression, 
      bulkExpressionIndependent=bulkExpressionIndependent,
      object=object, metadata=metadataList
    )
    return(returnList)
})

#' Method for \linkS4class{independentbulkParam}
#'
#' @param object An object of class \linkS4class{independentbulkParam} (see 
#' \code{?independentbulkParam}).
#' @details Display data summaries for an object of class 
#' \linkS4class{independentbulkParam}.
#' 
#' @examples 
#' new("independentbulkParam")
#'
#' @returns Shows object summaries.
#'
#' @export
setMethod("show", "independentbulkParam", function(object) {
  bulkExpression <- object[["bulkExpression"]]
  bulkExpressionIndependent <- object[["bulkExpressionIndependent"]]  
  samplesBulkExpression <- colnames(bulkExpression)
  samplesBulkExpressionIndependent <- colnames(bulkExpressionIndependent) 
  markersBulkExpression <- rownames(bulkExpression)
  markersBulkExpressionIndependent <- rownames(bulkExpressionIndependent)
  ## print info summaries
  message("Summary of independentbulkParam data:")
  message("\tNumber of unique sample IDs : ", 
          length(unique(markersBulkExpression, 
                        markersBulkExpressionIndependent)), "\n")
  message("\tNumber of unique marker IDs : ", 
          length(unique(samplesBulkExpression, 
                        samplesBulkExpressionIndependent)), "\n")
  message("\tNumber of independent samples : ", 
          length(
            samplesBulkExpressionIndependent[
              !samplesBulkExpressionIndependent %in% 
                samplesBulkExpression]), "\n")
})