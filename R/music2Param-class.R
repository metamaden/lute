#' music2Param-class
#' 
#' Main constructor for class to manage mappings to the deconvolution
#' method function \code{MuSiC::music2()}.
#' 
#' @include lute_generics.R
#' @include deconParam-class.R
#' @include referencebasedParam-class.R
#' @include independentbulkParam-class.R
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
#' @aliases 
#' MuSiC2Param-class, Music2Param-class
#'
setClass("music2Param", contains="independentbulkParam", slots=c(return.info = "logical"))
