#!/usr/bin/env R

# Functions to view the supported methods.
#
#
#


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
#' @param verbose Whether to show verbose status messages.
#' @examples
#' supported_strict_methods()
#' @returns List of supported strict deconvolution function names and 
#' descriptions.
#' @export 
supported_strict_methods <- function(varname = "strict_deconvolution_method_used",
                                     fname = "lute_transfer_learning"){
  if(verbose){message("Loading data from ",fname,"...")}
  ltl <- get(load(data(fname)))
  if(!varname %in% colnames(ltl)){
    stop("Error, ",varname," is not a variable in dataset.")}
  return(ltl[,varname])
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