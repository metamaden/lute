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
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a provided \linkS4class{ExpressionSet} 
#' containing single-cell data.
#' @param cellScaleFactors Cell size factor transformations of length equal to 
#' the K cell types to deconvolve.
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
#' Katharine M. Mullen and Ivo H. M. van Stokkum (2012). "nnls: The 
#' Lawson-Hanson algorithm for non-negative least squares (NNLS)." CRAN, R 
#' package version 1.4. URL: 
#' https://cran.r-project.org/web/packages/nnls/index.html
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

#' Show generic behavior for object of class nnlsParam
#' @param object Object of class \linkS4class{nnlsParam} (see 
#' \code{?nnlsParam}).
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' referencebasedParam(
#' bulkExpression=exampleList$bulkExpression, 
#' referenceExpression=exampleList$referenceExpression, 
#' cellScaleFactors=exampleList$cellScaleFactors)
#' 
#' @returns Prints data summary messages to console.
#' @export
setMethod("show", "nnlsParam", function(object) {
  # nnlsParam inherits from deconvolutionParam -> referencebasedParam
  # needs to show standard properties for each parent class
  
  ## deconvolutionParam -- show properties
  bulkExpression <- object[["bulkExpression"]]
  message("Object of class deconvolutionParam")
  message("\nData summaries:")
  message("\tNumber of bulk markers: ", nrow(bulkExpression))
  message("\tNumber of bulk samples: ", ncol(bulkExpression))
  markers <- rownames(bulkExpression)
  if(length(markers) > 10){markers <- markers[seq(10)]}
  message("\tFirst bulk marker labels:\n", 
          paste0(rownames(bulkExpression), collapse="; "))
  samples <- colnames(bulkExpression)
  if(length(samples) > 10){samples <- samples[seq(10)]}
  message("\tFirst sample labels:\n", paste0(samples, collapse="; "), "\n\n")
  
  ## referencebasedParam -- show properties
  ## get metadata
  cellScaleFactors <- object[["cellScaleFactors"]]
  bulkExpression <- object[["bulkExpression"]]
  referenceExpression <- object[["referenceExpression"]]
  uniqueTypes <- try(colnames(object[["referenceExpression"]]))
  markersBulkExpression <- rownames(bulkExpression)
  markersReferenceExpression <- rownames(referenceExpression)
  uniqueMarkers <- unique(c(markersBulkExpression, markersReferenceExpression))
  overlappingMarkers <- 
    intersect(markersBulkExpression, markersReferenceExpression)
  markerGenes <- nrow(referenceExpression)
  bulkSamples <- ncol(bulkExpression)
  numberCellTypesK <- ncol(referenceExpression)
  metadataList <- list(
    markerGenes = markerGenes, bulkSamples = bulkSamples, 
    numberCellTypesK = numberCellTypesK, cellScaleFactors = cellScaleFactors, 
    uniqueTypes = uniqueTypes, markersBulkExpression = markersBulkExpression, 
    markersReferenceExpression = markersReferenceExpression)
  ## post console messages
  cat(paste0("class: ", class(object)[1], "\n\n"))
  cat("key deconvolution run info:\n")
  cat("\tmarker info:\n")
  cat("\tsignature markers (Gz): ", markerGenes, "\n")
  cat("\tunique marker labels (Gy | Gz): ", length(uniqueMarkers), "\n")
  cat("\toverlapping marker labels (Gy & Gz): ", 
      length(overlappingMarkers), "\n\n")
  cat("\tsamples info:\n")
  cat("\tnumber of bulk samples (J): ", ncol(object[["bulkExpression"]]), "\n")
  cat("\tsample labels: ", 
      paste0(colnames(bulkExpression), collapse = "; "), "\n")
  cat("\n")
  cat("\tcell size factor properties:\n")
  if(!is(cellScaleFactors, "NULL")){
    for(type in names(cellScaleFactors)){
      cat("\tscale factor for type ", 
          type, ": ", cellScaleFactors["type"], "\n")}
    if(length(cellScaleFactors) == ncol(referenceExpression)){
      referenceExpression <- 
        .zstransform(referenceExpression, cellScaleFactors)}
  }; cat("\n")
  cat("\ttypes info:\n")
  cat("\tnumber of types (K): ", ncol(object[["referenceExpression"]]), "\n")
  if(!(is(uniqueTypes, "NULL")|is(uniqueTypes, "try-error"))){
    uniqueTypes <- uniqueTypes[order(uniqueTypes)]
    cat("\tunique type labels: ", paste0(uniqueTypes, collapse = ";"), "\n")
  } else{
    cat(
      "\nWarning, object 'referenceExpression' has no type labels (colnames)\n")
  }; cat("\n")
  ## parse additional warnings
  if(is(markersBulkExpression, "NULL")){
    cat("Warning, object 'bulkExpression' has no marker labels (rownames)\n\n")}
  if(is(markersReferenceExpression, "NULL")){
    cat(paste0("Warning, object 'referenceExpression' has no marker labels",
               " (rownames)\n\n"))}
  
})
