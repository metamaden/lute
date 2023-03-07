#' Run DeconRNASeq
#'
#' Runs the DeconRNASeq::DeconRNASeq() deconvolution algorithm.
#' 
#' @include lute_generics.R
#' 
#' @examples
#' # example
#' lexample <- .get_decon_example_data()
#' param <- musicParam(s = lexample[["s"]], y = lexample[["y"]], 
#' z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' deconvolution(param)
#' # type1     type2 
#' # 0.9819837 0.0180163
#' 
#' # return full results
#' param@return.info <- T
#' names(deconvolution(param))
#' # [1] "predictions" "result.info" "metadata"
#'
setClass("deconrnaseqParam", contains="deconParam", 
         slots=c(use.scale = "logical", return.info = "logical"))

#' Function to get nnlsParam
#' @export
#' @rdname AffinityParam-class
deconrnaseqParam <- function(y, z, s = NULL, use.scale = FALSE, return.info = FALSE) {
  if(is(use.scale, "NULL")){use.scale <- FALSE}
  new("deconrnaseqParam", y = y, z = z, s = s, 
      use.scale = use.scale, return.info = return.info)
}

#' Method for deconParam
#' @export
setMethod("deconvolution", signature(object = "deconrnaseqParam"), function(object){
  require(DeconRNASeq)
  lparam <- callNextMethod()
  # instantiate and format objects
  y <- lparam[["y"]]
  z <- lparam[["z"]]
  s <- lparam[["s"]]
  proportions <- matrix(s, nrow = 1)
  z = as.data.frame(z)
  y = as.data.frame(cbind(y,y))
  use.scale = object[["use.scale"]]
  result <- DeconRNASeq::DeconRNASeq(datasets = y,
                                     signatures = z,
                                     use.scale = use.scale,
                                     proportions = NULL)
  predictions <- result$out.all[1,]; names(predictions) <- colnames(z)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = lparam[["metadata"]])}
  return(lr)
})
