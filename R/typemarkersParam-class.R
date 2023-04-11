#' typemarkersParam-class
#'
#' Main constructor for class to manage mappings to the \code{typemarkers()} 
#' generic.
#' 
#' @include lute_generics.R
#'  
#' @details Main constructor for class \linkS4class{typemarkersParam}.
#' 
#' @rdname typemarkersParam-class
#' 
#' @seealso \linkS4class{meanratiosParam}
#' 
#' @param markers.per.type Number of top markers to get per cell type.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples
#' lexample <- lute:::.get_decon_example_data()
#' 
#' @aliases 
#' TypemarkersParam-class, TypeMarkersParam-class
#'
setClass("typemarkersParam", slots=c(markers.per.type = "numeric", 
                                     return.info = "logical"))

#' Make new object of class typemarkersParam
#'
#' Main constructor for class \linkS4class{typemarkersParam}.
#'
#' @param markers.per.type Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param return.info Whether to return metadata and original marker selection 
#' method outputs with predicted proportions.
#'
#' @returns New object of class \linkS4class{typemarkersParam}.
#'
#' @details This is the main parent class for cell type gene marker 
#' identification methods. Currently supported methods and their child classes include:
#' 
#' 1. Mean Ratios: The method DeconvoBuddies::get_mean_ratios2(), supported by the
#' class \linkS4class{meanratiosParam}.
#' 
#' @export
typemarkersParam <- function(markers.per.type = 20, return.info = FALSE) {
  new("typemarkersParam", markers.per.type = markers.per.type, 
      return.info = return.info)
}

#' Deconvolution method for class \linkS4class{deconrnaseqParam}
#' 
#' Main deconvolution method for the \linkS4class{deconrnaseqParam} to run the 
#' \code{DeconRNASeq::DeconRNASeq()} implementation of the DeconRNASeq algorithm.
#' 
#' @param object An object of class \linkS4class{deconrnaseqParam}.
#' 
#' @returns Either a vector of gene markers, or a list of detailed outputs that
#' includes such a marker vector.
#'
#' @export
setMethod("typemarkers", signature(object = "typemarkersParam"), function(object){
  lparam <- callNextMethod()
  # instantiate and format objects
  markers.per.type <- lparam[["markers.per.type"]]
  return.info <- lparam[["return.info"]]
})

#' Inspect slot in \linkS4class{typemarkersParam} object
#' @param x Object to access.
#' @param i Slot to access.
#' @returns Contents of specified slot.
#' @details Inspect slot in \linkS4class{typemarkersParam} object
#' @export
setMethod("[[", "typemarkersParam", function(x, i) {slot(x, i)})
