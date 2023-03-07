#' bisqueParam-class
#'
#' Applies the MuSiC::music.basic() implementation of the MuSiC deconvolution 
#' algorithm.
#' 
#' @include lute_generics.R
#' 
#' @details Main constructor for class \linkS4class{musicParam}.
#' @rdname musicParam-class
#' @seealso \linkS4class{deconParam}
#' 
#' @examples
#' # example
#' lexample <- .get_decon_example_data()
#' param <- bisqueParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@return.info <- T
#' names(deconvolution(param))
#' # [1] "predictions" "result.info" "metadata"
#' 
#' @aliases 
#' BisqueParam-class
#'
setClass("bisqueParam", contains="deconParam", slots=c(return.info = "logical"))

#' @export
bisqueParam <- function(y, z, s = NULL, return.info = FALSE) {
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
  new("bisqueParam", y = y, z = z, s = s, return.info = return.info)
}

#' @export
setMethod("deconvolution", signature(object = "bisqueParam"), function(object){
  require(BisqueRNA)
  require(Biobase)
  lparam <- callNextMethod()
  # instantiate objects
  y <- lparam[["y"]]
  z <- lparam[["z"]]
  s <- lparam[["s"]]
  # format objects
  y <- as.matrix(y)
  z <- as.matrix(z)
  s <- as.numeric(s)
  
  # get predictions
  predictions <- BisqueRNA::ReferenceBasedDecomposition()
  lr <- predictions
  
  
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = lparam[["metadata"]])}
  return(lr)
})


