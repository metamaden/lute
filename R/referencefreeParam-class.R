#' referencefreeParam-class
#'
#' Class and methods for managing reference-free deconvolution methods.
#' 
#' @inheritParams deconParam
#' @include deconParam-class.R
#' 
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
setClass("referencefreeParam", contains="deconParam", slots = c())