#!/usr/bin/env R

### Author: Sean Maden

#' deconvolutionParam-class
#' 
#' Defines the principal parent class for all deconvolution method parameters.
#' 
#' @include lute_generics.R
#' 
#' @details
#' Defines the parent class for deconvolution method parameters. Since all
#' deconvolution runs require a \code{y}  signals matrix, whether from 
#' experiment data or simulations such as pseudobulking, this parent class 
#' manages the bulk signals matrix. For this class, the 
#' \code{deconvolution} generic performs basic summaries of the bulk 
#' signals matrix.
#'
#' @rdname deconvolutionParam-class
#' @seealso 
#' \code{deconvolution}
#'
#' @examples 
#' param <- new("deconvolutionParam")
#' deconvolution(param)
#'
#' @aliases 
#' DeconvolutionParam-class, DeconParam-class, deconParam-class
#'
#' @returns New deconvolutionParam object.
#'
#' @export
setClass("deconvolutionParam",  slots=c(y="matrix", return.info="logical"))

#' Inspect slot in \linkS4class{deconvolutionParam} object
#' @param x Object to access.
#' @param i Slot to access.
#' @returns Contents of specified slot.
#' @details Inspect slot in \linkS4class{deconvolutionParam} object
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @returns Object slot contents.
#' 
#' @export
setMethod("[[", "deconvolutionParam", function(x, i) {slot(x, i)})

#' Deconvolution generic behavior for object of class \linkS4class{deconvolutionParam}
#' @param object An object of class \linkS4class{deconvolutionParam}.
#' @details Method for behavior of deconvolution generic when called for object of class 
#' \linkS4class{deconvolutionParam}
#' @returns Null method.
#' @export
setMethod("deconvolution", "deconvolutionParam", function(object) {})

#' Show generic behavior for object of class \linkS4class{deconvolutionParam}
#' @param object An object of class \linkS4class{deconvolutionParam}.
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{deconvolutionParam}
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @returns Shows object summaries.
#' 
#' @export
setMethod("show", "deconvolutionParam", function(object) {
  y <- object[["y"]]
  message("Object of class deconvolutionParam")
  message("\nData summaries:")
  message("\tNumber of bulk markers: ", nrow(y))
  message("\tNumber of bulk samples: ", ncol(y))
  markers <- rownames(y)
  if(length(markers) > 10){markers <- markers[seq(10)]}
  message("\tFirst bulk marker labels:\n", paste0(rownames(y), collapse="; "))
  samples <- colnames(y)
  if(length(samples) > 10){samples <- samples[seq(10)]}
  message("\tFirst sample labels:\n", paste0(samples, collapse="; "), "\n\n")
})
