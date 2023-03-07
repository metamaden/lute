#' bisqueParam-class
#'
#' Applies the BisqueRNA::ReferenceBasedDecomposition() implementation of the 
#' Bisque deconvolution algorithm.
#' 
#' @include lute_generics.R
#' 
#' @details Main constructor for class \linkS4class{bisqueParam}.
#' @rdname bisqueParam-class
#' @seealso \linkS4class{deconParam}
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
setClass("bisqueParam", contains="deconParam", slots=c(return.info = "logical"))

#' @export
bisqueParam <- function(y = NULL, z = NULL, s = NULL, 
                        y.eset = NULL, sc.eset = NULL, 
                        batch.variable = "batch.id", 
                        celltype.variable = "celltype",
                        return.info = FALSE) {
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
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
  if(!batch.variable %in% colnames(colData(sc.eset))){
    stop("Error, didn't find batch id variable ",batch.variable,
         " in sc.eset coldata.")
  }
  if(!celltype.variable %in% colnames(colData(sc.eset)){
    stop("Error, didn't find celltype id variable ", celltype.variable, 
         " in sc.eset coldata.")
  })
  if(is(z, "NULL")){
    message("Getting z from sc.eset...")
    z <- .get_z_from_sce(SingleCellExperiment(sc.eset))
  }
  message("Checking batch ids in bulk and sc eset...")
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
  new("bisqueParam", y = y, z = z, s = s, y.eset = y.eset, sc.eset = sc.eset, 
      return.info = return.info)
}

#' @export
setMethod("deconvolution", signature(object = "bisqueParam"), function(object){
  require(BisqueRNA)
  require(Biobase)
  lparam <- callNextMethod()
  # instantiate objects
  y.eset <- object[["y.eset"]]
  sc.eset <- object[["sc.eset"]]
  # get predictions
  result <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = y.eset, 
                                                        sc.eset = z.eset)
  lr <- predictions <- results$bulk.props
  # lr <- proportions$bulk.props[,1]
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, 
               result.info = result, 
               metadata = list(lmd = lparam[["metadata"]],
                               y.eset = y.eset,
                               sc.eset = sc.eset))}
  return(lr)
})

