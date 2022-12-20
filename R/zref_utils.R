#!/usr/bin/env R

# Defines utilities for making and managing the Z types reference matrix.
#

#' get_types
#'
#' Get type data from a SingleCellExperiment, returning either an object of 
#' either class SummarizedExperimentTypes or MultiAssayExperiment.
#'
#' @param sce Either a SummarizedExperiment or SingleCellExperiment object.
#' @param num.top.markers Number of top markers to retain by type (e.g. total 
#' markers is k*num.top.markers).
#' @param bind.marker.results Whether to bind valid flat outputs from marker
#' searches to the colData in the returned results object.
#' @param marker.method Valid method to identify markers from sce object.
#' @param type.method Valid method to obtain type-level expression/signals.
#' @param return.type Class of results object to return.
#' @param verbose Whether to show verbose status messages.
#' @param ... Additional arguments to pass for marker discovery method.
#' @returns Either a SummarizedExperimentTypes or or MultiAssayExperiment object.
#' @export
get_types <- function(sce, top.markers = 20, bind.marker.results = TRUE,
                      marker.method = "meanratio2", type.method = "mean", 
                      return.type = "SummarizedExperimentTypes", 
                      verbose = T, ...){
  
}

#' get_markerse
#'
#' Wrapper function for marker discovery procedures.
#'
#' @param marker.method Supported marker discovery methods.
#' 
get_markers <- function(marker.method = c("meanratio2")){
  
}




