#' Run MuSiC
#'
#' Applies the MuSiC::music.basic() implementation of the MuSiC deconvolution 
#' algorithm.
#' 
#' @include lute_generics.R
#' 
#' @examples
#' # example
#' lexample <- .get_decon_example_data()
#' param <- musicParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' deconvolution(param)
#' # type1     type2 
#' # 0.6770833 0.3229167
#' 
#' # return full results
#' param@return.info <- T
#' names(deconvolution(param))
#' # [1] "predictions" "result.info" "metadata"
#'
setClass("musicParam", contains="deconParam", 
         slots=c(sigma = "matrix", nu = "numeric", eps = "numeric", 
                 iter.max = "numeric", return.info = "logical"))

#' @export
musicParam <- function(y, z, s = NULL, sigma = NULL, nu = NULL, 
                       iter.max = NULL, eps = NULL, return.info = FALSE) {
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
  if(is(sigma, "NULL")){sigma <- matrix(0, ncol = 1, nrow = nrow(z))}
  if(is(nu, "NULL")){nu <- 1e-10}
  if(is(iter.max, "NULL")){iter.max <- 1000}
  if(is(eps, "NULL")){eps <- 0}
  new("musicParam", y = y, z = z, s = s, sigma = sigma, nu = nu, 
      iter.max = iter.max, eps = eps, return.info = return.info)
}

#' @export
setMethod("deconvolution", signature(object = "musicParam"), function(object){
  require(MuSiC)
  lparam <- callNextMethod()
  # instantiate objects
  nu <- object[["nu"]]
  iter.max <- object[["iter.max"]]
  eps <- object[["eps"]]
  sigma <- object[["sigma"]]
  y <- lparam[["y"]]
  z <- lparam[["z"]]
  s <- lparam[["s"]]
  # format objects
  y <- as.matrix(y)
  z <- as.matrix(z)
  s <- as.numeric(s)
  sigma <- as.matrix(sigma)
  result <- MuSiC::music.basic(X = z, 
                               Y = y, 
                               S = s, 
                               Sigma = sigma, 
                               nu = nu, 
                               iter.max = iter.max, 
                               eps = eps)
  predictions <- result$p.weight; names(predictions) <- colnames(z)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = lparam[["metadata"]])}
  return(lr)
})
