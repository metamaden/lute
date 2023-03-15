#' referencefreeParam-class
#'
#' Class and methods for managing reference-free deconvolution methods.
#' 
#' @details Reference-free deconvolution methods requently lack a signature 
#' matrix or use alternative means of prediction such as machine learning 
#' models. This parent class is meant to manage methods of this type.
#' 
#' @include deconvolutionParam-class.R
#' 
#' @examples 
#' param <- new("referencefreeParam)
#' deconvolution(param)
#' 
setClass("referencefreeParam", contains="deconvolutionParam", 
         slots = c(model.metadata = "list"))