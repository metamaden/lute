#' epicParam-class
#'
#' Runs the EPIC deconvolution algorithm.
#' 
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#' 
#' @details Main constructor for class \linkS4class{epicParam}.
#' @rdname epicParam-class
#' @seealso \linkS4class{epicParam}
#' 
#' @examples
#' # example
#' lexample <- lute:::.get_decon_example_data()
#' param <- epicParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' deconvolution(param)
#' 
#' # return full results
#' param@return.info <- TRUE
#' names(deconvolution(param))
#' # [1] "predictions" "result.info" "metadata"
#' 
#' @references 
#' 
#' Racle, J., Gfeller, D. (2020). EPIC: A Tool to Estimate the Proportions of 
#' Different Cell Types from Bulk Gene Expression Data. In: Boegel, S. (eds) 
#' Bioinformatics for Cancer Immunotherapy. Methods in Molecular Biology, vol 
#' 2120. Humana, New York, NY. https://doi-org.proxy1.library.jhu.edu/10.1007/978-1-0716-0327-7_17
#' 
#' Julien Racle and David Gfeller. EPIC: Estimate the Proportion of Immune and 
#' Cancer cells. (2022), GitHub, R package version 1.1.5. URL: https://github.com/GfellerLab/EPIC
#' 
#' @aliases 
#' EPICParam-class
#'
setClass("epicParam", contains="referencebasedParam", slots=c(z.var = "matrix"))

#' Make new object of class epicParam
#'
#' Main constructor for class \linkS4class{epicParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param z.var Signature variances matrix of same dimensions as z.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @returns New object of class \linkS4class{epicParam}.
#'
#' @details Takes standard inputs for the EPIC algorithm
#' 
#' @export
epicParam <- function(y, z, s = NULL, z.var = NULL, return.info = FALSE) {
  if(is(z.var, "NULL")){
    message("Making default `z.var` variance matrix...")
    z.var <- .get_zvar(z)
  }
  new("epicParam", y = y, z = z, s = s, z.var = z.var, return.info = return.info)
}

#' Deconvolution method for class \linkS4class{epicParam}
#' 
#' Main deconvolution method for the \linkS4class{epicParam} to run the 
#' \code{EPIC::EPIC()} implementation of the EPIC algorithm.
#' 
#' @param object An object of class \linkS4class{epicParam}.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#'
#' @references 
#' 
#' Racle, J., Gfeller, D. (2020). EPIC: A Tool to Estimate the Proportions of 
#' Different Cell Types from Bulk Gene Expression Data. In: Boegel, S. (eds) 
#' Bioinformatics for Cancer Immunotherapy. Methods in Molecular Biology, vol 
#' 2120. Humana, New York, NY. https://doi-org.proxy1.library.jhu.edu/10.1007/978-1-0716-0327-7_17
#' 
#' Julien Racle and David Gfeller. EPIC: Estimate the Proportion of Immune and 
#' Cancer cells. (2022), GitHub, R package version 1.1.5. URL: https://github.com/GfellerLab/EPIC
#'
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
  reference <- list(refProfiles = z,
                    sigGenes = rownames(z),
                    refProfiles.var = z.var)
  result <- EPIC::EPIC(bulk = y, 
                       reference = reference,
                       mRNA_cell = s)
  predictions <- result$mRNAProportions
  predictions <- apply(predictions, 1, function(ri){ri/sum(ri)})
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = list(lmd = lparam[["metadata"]],
                               z.final = lparam[["z"]]))}
  return(lr)
})
