#' bisqueParam-class
#'
#' Applies the BisqueRNA::ReferenceBasedDecomposition() implementation of the 
#' Bisque deconvolution algorithm.
#' 
#' @include lute_generics.R
#' @include deconParam-class.R
#' @include referencebasedParam-class.R
#' @include independentbulkParam-class.R
#' 
#' @details Main constructor for class \linkS4class{bisqueParam}.
#' @rdname bisqueParam-class
#' @seealso \linkS4class{deconParam}, \linkS4class{referencebasedParam}, \linkS4{independentbulkParam}
#' 
#' @examples
#' # example
#' lexample <- .get_decon_example_data_bisque()
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
setClass("bisqueParam", contains="independentbulkParam", 
         slots=c(y.eset = "ExpressionSet", sc.eset = "ExpressionSet", batch.variable = "character", 
                celltype.variable = "character", return.info = "logical"))

#' Make new object of class bisqueParam
#'
#' Main constructor for class \linkS4class{bisqueParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param yi Bulk mixed signals matrix of independent samples, which should not overlap samples in y.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param y.eset ExpressionSet of bulk mixed signals.
#' @param sc.eset ExpressionSet of single-cell transcriptomics data.
#' @param batch.variable Name of variable identifying the batches in sc.eset pData/coldata.
#' @param celltype.variable Name of cell type labels variable in sc.eset pData/coldata.
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @details Takes standard inputs for the Bisque method. If user provides matrices, will convert these
#' into ExpressionSet objects compatible with the main bisque method.
#' 
#' @export
bisqueParam <- function(y = NULL, yi = NULL, z = NULL, s = NULL, 
                        y.eset = NULL, sc.eset = NULL, batch.variable = "batch.id", 
                        celltype.variable = "celltype", return.info = FALSE) {
  
  
  # check y.eset/y
  if(is(y, "NULL")){
    if(is(y.eset, "NULL")){
      stop("Error, need to provide either y or bulk.eset.")
    } else{
      message("Getting y from provided bulk.eset...")
      y <- exprs(y.eset)
    }
  }
  if(is(y.eset, "NULL")){
    message("Making ExpressionSet from provided y...")
    y.assay <- y; y.colname <- colnames(y.assay)
    if(ncol(y) == 1){
      y.assay <- cbind(y, y)
      colnames(y.assay) <- c(y.colname, paste0(y.colname, "_rep1"))
    }
    df.y.pheno <- data.frame(SubjectName = colnames(y.assay))
    rownames(df.y.pheno) <- colnames(y)
    y.eset <- ExpressionSet(assayData = y.assay,
                            phenoData = AnnotatedDataFrame(df.y.pheno))
  }
  
  # check sc.eset
  if(is(sc.eset, "NULL")){
    stop("Error, no single-cell ExpressionSet provided.")  
    # add condition to call splatter simulations by default?
  }
  if(!batch.variable %in% colnames(pData(sc.eset))){
    stop("Error, didn't find batch id variable ",batch.variable,
         " in sc.eset pData/coldata.")
  }
  if(!celltype.variable %in% colnames(pData(sc.eset))){
    stop("Error, didn't find celltype id variable ", celltype.variable, 
         " in sc.eset pData/coldata.")
  }
  if(is(z, "NULL")){
    message("Getting z from sc.eset...")
    z <- .get_z_from_sce(SingleCellExperiment(sc.eset))
  }
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
  message("Checking batch ids in bulk and sc eset...")
  cond <- !batch.variable %in% colnames(pData(sc.eset))|
    !batch.variable %in% colnames(pData(y.eset))
  if(cond){stop("Error, didn't find batch variable in sc.eset or y.eset pData: ", batch.variable)}
  id.sc <- unique(sc.eset[[batch.variable]])
  id.bulk <- unique(y.eset[[batch.variable]])
  id.overlap <- intersect(id.sc, id.bulk)
  id.unique <- unique(c(id.sc, id.bulk))
  id.onlybulk <- id.bulk[!id.bulk %in% id.overlap]
  id.onlysc <- id.sc[!id.sc %in% id.overlap]
  message("Found ", length(id.unique), " unique batch ids...")
  message("Found ", length(id.overlap), " overlapping batch ids...")
  message("Found ", length(id.onlybulk), " bulk-only batch ids...")
  message("Found ", length(id.onlybulk), " sc-only batch ids...")
  if(length(id.overlap) == 0){
    stop("Error, no overlapping markers in y.eset and sc.eset.")
  } else if(length(id.unique) >= length(id.bulk)){
    stop("Error, no unique bulk samples provided which aren't also in sc.eset.")
  } else{
    message("Finished validating batch ids.")
  }
  new("bisqueParam", y = y, z = z, s = s, y.eset = y.eset, sc.eset = sc.eset,
      batch.variable = batch.variable, celltype.variable = celltype.variable,
      return.info = return.info)
}

#' Deconvolution method for bisqueParam
#'
#' Main method to access the Bisque deconvolution method from the main lute deconvolution genetic.
#'
#' @details Takes an object of class bisqueParam as input, returning a list.
#' @returns Either a vector of predicted proportions, or a list containing predictions, metadata, 
#' and original outputs.
#'
#' @export
setMethod("deconvolution", signature(object = "bisqueParam"), function(object){
  require(BisqueRNA); require(Biobase)
  lparam <- callNextMethod()
  # instantiate objects
  y.eset <- object[["y.eset"]]; sc.eset <- object[["sc.eset"]]
  # get predictions
  result <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = y.eset, sc.eset = z.eset)
  lr <- predictions <- results$bulk.props
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = list(lmd = lparam[["metadata"]], y.eset = y.eset, sc.eset = sc.eset))}
  return(lr)
})
