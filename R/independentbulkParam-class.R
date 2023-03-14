#' independentbulkParam-class
#'
#' Class and methods for managing methods requiring independent bulk samples.
#' 
#' @include referencebasedParam-class.R
#'
#' @param yi Mixed signals matrix from bulk samples, independent from primary mixed signals matrix y.
#' 
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
setClass("independentbulkParam", contains="referencebasedParam", slots = c(yi = "matrix"))