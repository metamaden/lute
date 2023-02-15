#!/usr/bin/env r

# Author: Sean Maden
#
# Main functions managing mappings between standard deconvolution interface and
# specific functions.
#

#' run_deconvolution
#'
#'
run_deconvolution <- function(method = "nnls", ...){
  mappings <- map_deconvolution_arguments(...)
  results <- get_deconvolution_predictions(method = method, mappings = mappings)
  return(results)
}

#' map_deconvolution_arguments
#'
#'
#'
map_deconvolution_arguments <- function(){
  
}

#' get_deconvolution_predictions
#'
#'
#'
get_deconvolution_predictions <- function(){
  
}

#-------------------
# main use functions
#-------------------

#' use_nnls
#'
#'
#'
use_nnls <- function(){
  
}

#' use_music
#'
#'
#'
use_music <- function(){
  
}

#' use_bisque
#'
#'
use_bisque <- function(){
  
}

