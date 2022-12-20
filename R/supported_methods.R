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
#' @returns Vector of supported strict deconvolution function names and 
#' descriptions.
#' @export 
supported_strict_methods <- function(varname = "strict_deconvolution_method_used",
                                     fname = "lute_deconvolution_transfer_learning",
                                     verbose = FALSE){
  if(verbose){message("Loading data from ",fname,"...")}
  ltl <- get(load(data("lute_deconvolution_transfer_learning")))
  if(!varname %in% colnames(ltl)){
    stop("Error, ",varname," is not a variable in dataset.")}
  return(ltl[,varname])
}

#' supported_marker_methods
#'
#' Show supported marker identification methods. These are load from the 
#' provided table in the lute package.
#'
#' @param varname Variable/column in lute_transfer_learning.rda dataset with 
#' strict deconvolution method names.
#' @param fname Name of the file containing supported deconvolution methods 
#' info.
#' @param verbose Whether to show verbose status messages.
#' @examples
#' supported_marker_methods()
#' @returns Vector of supported marker selection methods.
#' @seealso supported_strict_methods
#' @export
supported_marker_methods <- function(varname = "marker_method_name", 
                                     fname = "lute_marker-method_transfer_learning", 
                                     verbose = FALSE){
  if(verbose){message("Loading data from ",fname,"...")}
  # data("lute_marker-method_transfer_learning")
  eval(parse(text = paste0("data(", fname, ")")))
  if(!varname %in% colnames(dfmarker)){
    stop("Error, ",varname," is not a variable in dataset.")}
  return(dfmarker[,varname])
}
