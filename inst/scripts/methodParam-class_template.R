#' methodParam-class
#' 
#' Main constructor for class to manage mappings to the deconvolution
#' method function \code{\link{LIBRARY::METHOD}}.
#' 
#' @include lute_generics.R
#' @include deconParam-class.R
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
#' @aliases 
#' MethodParam-class
#'
# setClass("methodParam", contains="deconParam", slots=c(return.info = "logical"))

# #' @export
# methodParam <- function(y, z, s = NULL, return.info = FALSE) {
#  if(is(s, "NULL")){s <- matrix(rep(1, ncol(z)), nrow = 1)}
#  new("methodParam", s = s, y = y, z = z, return.info = return.info)
#}

# #' @export
# setMethod("deconvolution", signature(object = "nnlsParam"), function(object){
#  require(nnls); lparam <- callNextMethod()
#  y <- lparam[["y"]]; z <- lparam[["z"]]; s <- lparam[["s"]]
#  y <- as.matrix(y); z <- as.matrix(z); s <- as.numeric(s)
#  result <- # deconvolution function goes here, e.g. nnls::nnls(A = z, b = y)
#  predictions <- # call to get pred proportions goes here, e.g. result$x
#  names(predictions) <- colnames(z); lr <- predictions
#  if(object[["return.info"]]){
#    lr <- list(predictions = predictions, result.info = result, 
#               metadata = lparam[["metadata"]])}
#  return(lr)
#})