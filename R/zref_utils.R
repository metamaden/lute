#!/usr/bin/env R

# Author: Sean Maden
#
# Defines utilities for making and managing the Z types reference matrix.
#

#' get_types
#'
#' Get type data from a SingleCellExperiment, returning either an object of 
#' either class SummarizedExperimentTypes or MultiAssayExperiment.
#'
#' @param sce
#' @param top.markers
#' @param bind.marker.results
#' @param marker.method
#' @param type.method
#' @param return.type
#' @param verbose
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




