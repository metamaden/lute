#' nnlsParam-class
#'
#' Uses nnls::nnls().
#' 
#' @include lute_generics.R
#' @include deconParam-class.R
#' 
#' @details Main constructor for class \linkS4class{nnlsParam}.
#' @rdname nnlsParam-class
#' @seealso \linkS4class{deconParam}
#' 
#' @examples 
#' lexample <- .get_decon_example_data()
#' param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], 
#' z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' deconvolution(param)
#' # type1      type2 
#' # 0.48908543 0.05896868
#' 
#' # return full results
#' param@return.info <- T
#' names(deconvolution(param))
#' # [1] "predictions" "result.info" "metadata"
#' 
#' @aliases 
#' NNLSParam-class
#' 
setClass("nnlsParam", contains="deconParam", slots=c(return.info = "logical"))

#' @export
nnlsParam <- function(y, z, s = NULL, return.info = FALSE) {
  if(is(s, "NULL")){s <- matrix(rep(1, ncol(z)), nrow = 1)}
  new("nnlsParam", s = s, y = y, z = z, return.info = return.info)
}

#' @export
setMethod("deconvolution", signature(object = "nnlsParam"), function(object){
  require(nnls)
  lparam <- callNextMethod()
  y <- lparam[["y"]]; z <- lparam[["z"]]; s <- lparam[["s"]]
  y <- as.matrix(y)
  z <- as.matrix(z)
  s <- as.numeric(s)
  result <- nnls::nnls(A = z, b = y)
  predictions <- result$x; names(predictions) <- colnames(z)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = lparam[["metadata"]])}
  return(lr)
})
