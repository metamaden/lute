#' Non-negative least squares
#'
#' Uses nnls::nnls().
#' 
#' @inheritParams deconParam
#' 
setClass("nnlsParam", contains="deconParam", 
         slots=c(s="matrix", y="matrix", z="matrix",
                 return.info = "logical"))

#' Function to get nnlsParam
#' @export
#' @rdname AffinityParam-class
nnlsParam <- function(y, z, s = NULL, return.info = FALSE) {
  if(is(s, "NULL")){s <- matrix(rep(1, ncol(z)), nrow = 1)}
  new("nnlsParam", s = s, y = as.matrix(y), z = as.matrix(z), 
      return.info = return.info)
}

#' Method for deconParam
#' @export
setMethod("deconvolution", signature(object = "nnlsParam"), function(object){
  require(nnls)
  lr <- callNextMethod()
  y <- lr[["y"]]; z <- lr[["z"]]; s <- lr[["s"]]
  result <- nnls::nnls(A = z, b = y)
  predictions <- result$x; names(predictions) <- colnames(z)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = lr[["metadata"]])}
  return(lr)
})
