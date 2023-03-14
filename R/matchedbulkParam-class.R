#' matchedbulkParam-class
#'
#' Class and methods for managing methods which allow for bulk sample data to be matched with signatures data.
#' 
#' @include referencebasedParam-class.R
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
setClass("matchedbulkParam", contains="referencebasedParam")