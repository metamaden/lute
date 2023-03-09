#' deconParam-class
#' 
#' Defines the parent class for deconvolution method parameters.
#' 
#' @include lute_generics.R
#' 
#' @details
#' Defines the parent class for deconvolution method parameters. For this class,
#' the \link{deconvolution} generic performs several useful housekeeping
#' operations for a standard deconvolution run, including inspecting and 
#' summarizing properties of the objects z, y, and s.
#'
#' @rdname deconParam-class
#' @seealso 
#' \link{\code{deconvolution}}
#'
#' @examples 
#' deconparam <- new("deconParam")
#' deconvolution(deconparam)
#'
#' @aliases 
#' DeconParam-class

#' @export
setClass("deconParam",  slots=c(y="matrix", return.info = "logical"))

#' @export
setMethod("[[", "deconParam", function(x, i) {slot(x, i)})

#' @export
setMethod("deconvolution", "deconParam", function(object) {})

#' @export
setMethod("show", "deconParam", function(object) {})