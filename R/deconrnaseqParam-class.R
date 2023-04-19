#' deconrnaseqParam-class
#'
#' Main constructor for class to manage mappings to the deconvolution
#' method function \code{DeconRNASeq::DeconRNASeq()}.
#' 
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#' 
#' @details Main constructor for class \linkS4class{deconrnaseqParam}.
#' @rdname deconrnaseqParam-class
#' @seealso \linkS4class{deconParam}
#' 
#' @examples
#' lexample <- lute:::.get_decon_example_data() # get example data 
#' param <- deconrnaseqParam(s = lexample[["s"]], y = lexample[["y"]], 
#' z = lexample[["z"]])
#' 
#' # return only predicted proportions
#' deconvolution(param)
#' 
#' @aliases 
#' DeconRNASeqParam-class
#'
setClass("deconrnaseqParam", contains="referencebasedParam", 
         slots=c(use.scale = "logical"))

#' Make new object of class deconrnaseqParam
#'
#' Main constructor for class \linkS4class{deconrnaseqParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param use.scale Whether to rescale signature matrix data (loigical).
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @returns New object of class \linkS4class{deconrnaseqParam}.
#'
#' @details Takes standard inputs for the DeconRNASeq algorithm
#' 
#' @export
deconrnaseqParam <- function(y, z, s = NULL, use.scale = FALSE, return.info = FALSE) {
  if(is(use.scale, "NULL")){use.scale <- FALSE}
  new("deconrnaseqParam", y = y, z = z, s = s, 
      use.scale = use.scale, return.info = return.info)
}

#' Deconvolution method for class \linkS4class{deconrnaseqParam}
#' 
#' Main deconvolution method for the \linkS4class{deconrnaseqParam} to run the 
#' \code{DeconRNASeq::DeconRNASeq()} implementation of the DeconRNASeq algorithm.
#' 
#' @param object An object of class \linkS4class{deconrnaseqParam}.
#' 
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#'
#' @references 
#' 
#' Ting Gong and Joseph D. Szustakowski. DeconRNASeq: Deconvolution of 
#' Heterogeneous Tissue Samples for mRNA-Seq data. (2022), Bioconductor, 
#' R package version 1.38.0. DOI: 10.18129/B9.bioc.DeconRNASeq  
#' 
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
  result <- DeconRNASeq::DeconRNASeq(datasets = y, signatures = z,
                                     use.scale = use.scale, proportions = NULL)
  predictions <- matrix(result$out.all[1,], ncol = ncol(z))
  colnames(predictions) <- colnames(z)
  rownames(predictions) <- colnames(y)
  lr <- predictions
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = lparam[["metadata"]])}
  return(lr)
})
