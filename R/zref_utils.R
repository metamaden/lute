#!/usr/bin/env R

# Defines utilities for making and managing the Z types reference matrix.
#

#' compute_type_se
#'
#' Calculate a new types object, SummarizedExperimentTypes, from a 
#' SingleCellExperiment. Returns an object of either class 
#' SummarizedExperimentTypes or MultiAssayExperiment.
#'
#' @param sce Either a SummarizedExperiment or SingleCellExperiment object.
#' @param num.top.markers Number of top markers to retain by type (e.g. total 
#' markers is k*num.top.markers).
#' @param marker.varname Marker variable name from rowData(sce).
#' @param type.varname Type variable name from colData(sce).
#' @param bind.marker.results Whether to bind valid flat outputs from marker
#' searches to the colData in the returned results object.
#' @param marker.method Valid method to identify markers from sce object.
#' @param type.method Valid method to obtain type-level expression/signals.
#' @param return.type Class of results object to return.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments to pass for marker discovery method.
#' @returns Either a SummarizedExperimentTypes or or MultiAssayExperiment object.
#' @export
compute_type_se <- function(sce, marker.varname = NULL, type.varname = "cellType", 
                            num.top.markers = 20, bind.marker.results = TRUE, 
                            marker.method = "meanratio2", type.method = "mean", 
                            return.type = "SummarizedExperimentTypes", 
                            verbose = T, ...){
  
  # manage type-specific markers
  if(is(marker.varname, "NULL")){user.markers <- 0}
  while(user.markers == 0)
    marker.cond <- !marker.varname %in% colnames(rowData(sce))
    if(marker.cond){
      if(verbose){message("Didn't find markers variable.")}
      marker.varname <- "typeMarkers"
    } else{user.markers <- 1}
    if(user.markers == 0){
      valid.markers <- supported_marker_methods(verbose = verbose)
      if(marker.method %in% valid.markers){
        if(verbose){
          message("Calculating markers using method ", marker.method, "...")}
        sce <- sce_append_markers()
      } else{
        stop("Provided method ", marker.method, " is not a valid marker method.")
      }
    }
    repeat
  # get the type-specific marker signal summaries
  set.assays <- kexpr_sce()
  # form new set object
  set <- SummarizedExperimentTypes(set.assays = set.assays)
  return(set)
}

#' sce_append_markers
#'
#' Wrapper function for marker discovery procedures. Accepts and returns a 
#' SingleCellExperiment, with marker labels appended.
#'
#' @param marker.method Supported marker discovery methods.
#' @param verbose Whether to show verbose status messages.
#' 
sce_append_markers <- function(sce, marker.varname = "typeMarker", 
                               marker.method = c("meanratio2")){
  
  
  return(sce.new)
}
