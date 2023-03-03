#' Non-negative least squares
#'
#' Uses nnls::nnls().
#' 
#' @inheritParams deconParam
#' 
#' @examples
#' set.seed(0)
#' y <- matrix(rnbinom(n=10, size=10, mu=10), ncol = 1)
#' z <- matrix(rnbinom(n=20, size=10, mu=10), ncol = 2)
#' colnames(z) <- paste0("type", seq(ncol(z)))
#' s <- c(1, 10)
#' param <- musicParam(s = NULL, y = y, z = z)
#' res <- deconvolution(param)
#'
setClass("musicParam", contains="deconParam", 
         slots=c(s="numeric", y="matrix", z="matrix", sigma = "matrix",
                 nu = "numeric", eps = "numeric", iter.max = "numeric",
                 return.info = "logical"))

#' Function to get nnlsParam
#' @export
#' @rdname AffinityParam-class
musicParam <- function(y, z, s = NULL, sigma = NULL, nu = NULL, 
                       iter.max = NULL, eps = NULL, return.info = FALSE) {
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
  if(is(sigma, "NULL")){sigma <- matrix(0, ncol = 1, nrow = nrow(z))}
  if(is(nu, "NULL")){nu <- 1e-10}
  if(is(iter.max, "NULL")){iter.max <- 1000}
  if(is(eps, "NULL")){eps <- 0}
  new("musicParam", y = as.matrix(y), z = as.matrix(z), s = as.numeric(s), 
      sigma = as.matrix(sigma), nu = as.numeric(nu), iter.max = iter.max, 
      eps = eps)
}

#' Method for deconParam
#' @export
setMethod("deconvolution", signature(object = "musicParam"), function(object){
  require(MuSiC)
  ldecon <- callNextMethod()
  y <- ldecon[["y"]]
  # z <- .zstransform(object[["z"]], object[["s"]])
  result <- MuSiC::music.basic(X = object[["z"]], 
                               Y = y, 
                               S = object[["s"]], 
                               Sigma = object[["sigma"]], 
                               nu = object[["nu"]], 
                               iter.max = object[["iter.max"]], 
                               eps = object[["eps"]])
  predictions <- result$p.weight; names(predictions) <- colnames(z)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = ldecon[["metadata"]])}
  return(predictions)
})

param <- musicParam(y = y, z = z)
res <- deconvolution(param)

