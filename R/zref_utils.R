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
  
  
  if(is(marker.varname, "NULL")){user.markers <- 0}
  while(user.markers == 0)
    marker.cond <- !marker.varname %in% colnames(rowData(sce))
    if(marker.cond){
      if(verbose){message("Didn't find markers variable.")}
      marker.varname <- "typeMarkers"
    } else{user.markers <- 1}
    if(user.markers == 0){
      valid.markers <- get_valid_marker_method()
      if(marker.method %in% valid.markers){
        if(verbose){
          message("Calculating markers using method ", marker.method, "...")}
        sce <- sce_append_markers()
      } else{
        stop("Provided method ", marker.method, " is not a valid marker method.")
      }
    }
    repeat
  return(set)
}


sce_append_markers <- function(sce, marker.varname = "typeMarker", 
                               marker.method = c("meanratio2")){
  return(sce.new)
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

#' get_valid_marker_method
#'
#' Show supported marker identification methods. These are load from the 
#' provided table in the lute package.
#'
#' @seealso supported_strict_methods
#' @export
get_valid_marker_method <- function(verbose = FALSE){
  if(verbose){message("Loading data from ",fname,"...")}
  ltl <- get(load(data(fname)))
  if(!varname %in% colnames(ltl)){
    stop("Error, ",varname," is not a variable in dataset.")}
  return(ltl[,varname])
}

#' supported_strict_methods
#' 
#' Show supported strict deconvolution function. That is, functions which take
#' some reference/signature matrix Z and some bulk/convoluted signals matrix Y
#' and return the vector of predictions or proportions, p.
#' 
#' @param varname Variable/column in lute_transfer_learning.rda dataset with 
#' strict deconvolution method names.
#' @param fname Name of the file containing supported deconvolution methods 
#' info.
#' @examples
#' supported_strict_methods()
#' @returns List of supported strict deconvolution function names and 
#' descriptions.
#' @export 
supported_strict_methods <- function(varname = "strict_deconvolution_method_used",
                                     fname = "lute_transfer_learning",
                                     varbose = FALSE){
  if(verbose){message("Loading data from ",fname,"...")}
  ltl <- get(load(data(fname)))
  if(!varname %in% colnames(ltl)){
    stop("Error, ",varname," is not a variable in dataset.")}
  return(ltl[,varname])
}

