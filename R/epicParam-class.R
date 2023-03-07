#' epicParam-class
#'
#' Runs the EPIC::EPIC() deconvolution algorithm.
#' 
#' @include lute_generics.R
#' 
#' @details Main constructor for class \linkS4class{epicParam}.
#' @rdname epicParam-class
#' @seealso \linkS4class{epicParam}
#' 
#' @examples
#' # example
#' lexample <- .get_decon_example_data()
#' param <- epicParam(s = lexample[["s"]], y = lexample[["y"]], 
#' z = lexample[["z"]])
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
#' EPICParam-class
#'
setClass("epicParam", contains="deconParam", 
         slots=c(z.var = "matrix", return.info = "logical"))

#' @export
epicParam <- function(y, z, s = NULL, z.var = NULL, return.info = FALSE) {
  if(is(z.var, "NULL")){
    message("Making default `z.var` variance matrix...")
    z.var <- .get_zvar(z)
  }
  new("epicParam", y = y, z = z, s = s, z.var = z.var, return.info = return.info)
}

#' @export
setMethod("deconvolution", signature(object = "epicParam"), function(object){
  require(EPIC)
  lparam <- callNextMethod()
  # instantiate and format objects
  y <- lparam[["y"]]
  z <- lparam[["z"]]
  z.var <- object[["z.var"]]
  s <- lparam[["s"]]
  s <- as.numeric(s)
  if(!"otherCells" %in% names(s)){
    message("Setting size/mRNA for missing label 'otherCells' to 0...")
    s["otherCells"] <- 0
  }
  z <- as.matrix(z)
  y <- as.data.frame(y)
  z.var <- as.matrix(z.var)
  # get reference object
  reference <- list(refProfiles = z,
                    sigGenes = rownames(z),
                    refProfiles.var = z.var)
  
  # get predictions
  result <- EPIC::EPIC(bulk = y, 
                       reference = reference,
                       mRNA_cell = s)
  predictions <- result$mRNAProportions
  names(predictions) <- colnames(z)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = list(lmd = lparam[["metadata"]],
                               z.final = lparam[["z"]]))}
  return(lr)
})
