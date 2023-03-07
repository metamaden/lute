#' scdcParam-class
#' 
#' Main constructor for class to manage mappings to the deconvolution
#' method function \code{SCDC::SCDC_prop()}.
#' 
#' @include lute_generics.R
#' @include deconParam-class.R
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
#' @aliases 
#' SCDCParam-class, ScdcParam-class
#'
setClass("scdcParam", contains="deconParam", slots=c(return.info = "logical"))
