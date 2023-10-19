#!/usr/bin/env R

### Author: Sean Maden

#' referencebasedParam-class
#'
#' Class and methods for managing reference-based deconvolution methods.
#' 
#' @include deconvolutionParam-class.R
#' 
#' @details This is a parent class to manage reference-based deconvolution 
#' algorithms. 
#' 
#' Child/sub-classes of this are distinguished by their use of
#' either an explicit or implied \code{z} signature matrix (i.e. Z[G,K] of
#' dimensions G markers by K cell types). These also have an implied cell size 
#' term for biases from systematic cell size differences. If no cell size 
#' transformation is intended, this is the equivalent of passing equal size 
#' scales, (e.g. a K-length vector of equal values). See 
#' `vignette(package="lute")` for details about experiment terms.
#' 
#' @examples 
#' exampleList <- get_decon_example_data()
#' referencebasedParam(
#' bulkExpression=lexample$bulkExpression, 
#' referenceExpression=lexample$referenceExpression, 
#' cellScaleFactors=lexample$cellScaleFactors)
#'
#' @returns New object.
#' 
setClass("referencebasedParam", contains="deconvolutionParam", 
         slots=c(z = "matrix", s = "numeric"))

#' Make new object of class referencebasedParam
#'
#' Main constructor for class \linkS4class{referencebasedParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, can 
#' be computed from a provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types 
#' to deconvolve.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#'
#' @examples 
#' exampleList <- get_decon_example_data()
#' referencebasedParam(
#' bulkExpression=lexample$bulkExpression, 
#' referenceExpression=lexample$referenceExpression, 
#' cellScaleFactors=lexample$cellScaleFactors)
#'
#' @returns New object of class \linkS4class{referencebasedParam}.
#'
#' @details Takes standard inputs for reference-based deconvolution algorithms.
#'
#' @returns New object.
#' 
#' @export
referencebasedParam <- function(y, z, s, return.info = FALSE) {
  new("referencebasedParam", y = y, z = z, s = s, return.info = return.info)
}

#' Deconvolution generic behavior for object of class 
#' \linkS4class{referencebasedParam}
#' 
#' @param object An object of class \linkS4class{referencebasedParam} (see 
#' \code{?referencebasedParam}).
#' @details Method for behavior of deconvolution generic when called for object 
#' of class \linkS4class{referencebasedParam}.
#' 
#' @examples 
#' exampleList <- get_decon_example_data()
#' referencebasedParam(
#' bulkExpression=lexample$bulkExpression, 
#' referenceExpression=lexample$referenceExpression, 
#' cellScaleFactors=lexample$cellScaleFactors)
#' 
#' @returns Method results.
#' @export
setMethod("deconvolution", "referencebasedParam", function(object) {
  ## get metadata
  cellScaleFactors <- object[[cellScaleFactors]]
  bulkExpression <- object[[bulkExpression]]
  referenceExpression <- object[[referenceExpression]]
  ## cell types in z, s
  if(is(cellScaleFactors, "NULL")){
    cellScaleFactors <- rep(1, ncol(referenceExpression))}
  uniqueTypes <- try(colnames(referenceExpression))
  condition.z.types <- is(uniqueTypes, "NULL")|is(uniqueTypes, "try-error")
  if(!condition.z.types){
    uniqueTypes <- uniqueTypes[order(uniqueTypes)]
    referenceExpression <- 
      referenceExpression[,order(colnames(referenceExpression), uniqueTypes)]
    condition.s.types <- is(names(cellScaleFactors), "NULL")
    if(!condition.s.types){
      filter.s.types <- names(cellScaleFactors) %in% uniqueTypes
      cellScaleFactors <- cellScaleFactors[filter.s.types]
      cellScaleFactors <- 
        cellScaleFactors[order(names(cellScaleFactors), uniqueTypes)]
    }
  }
  referenceExpression <- .zstransform(referenceExpression, cellScaleFactors)
  ## matching markers in y and z
  markersBulkExpression <- rownames(bulkExpression)
  markersReferenceExpression <- rownames(referenceExpression)
  if(!is(markersBulkExpression, "NULL") & 
     !is(markersReferenceExpression, "NULL")){
    uniqueMarkers <- unique(
      c(markersBulkExpression, markersReferenceExpression))
    overlappingMarkers <- intersect(
      markersBulkExpression, markersReferenceExpression)
    y.filter <- rownames(bulkExpression) %in% overlappingMarkers
    z.filter <- rownames(referenceExpression) %in% overlappingMarkers
    bulkExpression <- bulkExpression[y.filter,,drop=FALSE]
    referenceExpression <- referenceExpression[z.filter,,drop=FALSE]
    bulkExpression <- bulkExpression[
      order(match(rownames(bulkExpression), overlappingMarkers)),]
    referenceExpression <- referenceExpression[
      order(match(rownames(referenceExpression), overlappingMarkers)),]
  } else{
    message("Warning, rownames not provided in both y and z. ",
            "Can't match marker labels.")
  }
  ## parse additional warnings
  if(is(markersBulkExpression, "NULL")){
    message("Warning, object 'y' has no marker labels (rownames)\n")}
  if(is(markersReferenceExpression, "NULL")){
    message("Warning, object 'z' has no marker labels (rownames)\n")}
  ## get final metadata
  markerGenes <- nrow(referenceExpression)
  bulkSamples <- ncol(bulkExpression)
  numberCellTypesK <- ncol(referenceExpression)
  metadataList <- list(markerGenes = markerGenes, bulkSamples = bulkSamples, 
                        numberCellTypesK = numberCellTypesK, 
                        cellScaleFactors = cellScaleFactors, 
                        uniqueTypes = uniqueTypes,
                        markersBulkExpression = markersBulkExpression, 
                        markersReferenceExpression = markersReferenceExpression)
  ## return list
  return(
    list(bulkExpression = as.matrix(bulkExpression), 
         referenceExpression = as.matrix(referenceExpression), 
         cellScaleFactors = as.numeric(cellScaleFactors), 
         metadata = metadataList)
    )
})

#' Show generic behavior for object of class referencebasedParam
#' @param object Object of class \linkS4class{referencebasedParam} (see 
#' \code{?referencebasedParam}).
#' 
#' @examples 
#' exampleList <- get_decon_example_data()
#' referencebasedParam(
#' bulkExpression=lexample$bulkExpression, 
#' referenceExpression=lexample$referenceExpression, 
#' cellScaleFactors=lexample$cellScaleFactors)
#' 
#' @returns Prints data summary messages to console.
#' @export
setMethod("show", "referencebasedParam", function(object) {
  ## get metadata
  cellScaleFactors <- object[[cellScaleFactors]]
  bulkExpression <- object[[bulkExpression]]
  referenceExpression <- object[[referenceExpression]]
  uniqueTypes <- try(colnames(object[[referenceExpression]]))
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
  cat("\tnumber of bulk samples (J): ", ncol(object[[bulkExpression]]), "\n")
  cat("\tsample labels: ", 
      paste0(colnames(bulkExpression), collapse = "; "), "\n")
  cat("\n")
  cat("\tcell size factor properties:\n")
  if(!is(cellScaleFactors, "NULL")){
    for(type in names(cellScaleFactors)){
      cat("\tscale factor for type ", type, ": ", cellScaleFactors[type], "\n")}
    if(length(cellScaleFactors) == ncol(referenceExpression)){
      referenceExpression <- .zstransform(referenceExpression, cellScaleFactors)}
  }; cat("\n")
  cat("\ttypes info:\n")
  cat("\tnumber of types (K): ", ncol(object[[referenceExpression]]), "\n")
  if(!(is(uniqueTypes, "NULL")|is(uniqueTypes, "try-error"))){
    uniqueTypes <- uniqueTypes[order(uniqueTypes)]
    cat("\tunique type labels: ", paste0(uniqueTypes, collapse = ";"), "\n")
  } else{
    cat("\nWarning, object 'z' has no type labels (colnames)\n")
  }; cat("\n")
  ## parse additional warnings
  if(is(markersBulkExpression, "NULL")){
    cat("Warning, object 'y' has no marker labels (rownames)\n\n")}
  if(is(markersReferenceExpression, "NULL")){
    cat("Warning, object 'z' has no marker labels (rownames)\n\n")}
})
