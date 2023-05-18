#' referencefreeParam-class
#'
#' Class and methods for managing reference-free deconvolution methods.
#' 
#' @include deconvolutionParam-class.R
#' 
#' @details Reference-free deconvolution methods requently lack a signature 
#' matrix or use alternative means of prediction such as machine learning 
#' models. This parent class is meant to manage methods of this type.
#' 
#' @examples 
#' new("referencefreeParam")
#' 
setClass("referencefreeParam", contains="deconvolutionParam", 
         slots = c(model.metadata = "list"))