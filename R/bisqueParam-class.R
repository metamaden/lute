#' bisqueParam-class
#'
#' Applies the BisqueRNA::ReferenceBasedDecomposition() implementation of the 
#' Bisque deconvolution algorithm.
#' 
#' @include lute_generics.R
#' @include independentbulkParam-class.R
#' 
#' @details Main constructor for class \linkS4class{bisqueParam}.
#' @rdname bisqueParam-class
#' @seealso \linkS4class{deconvolutionParam}, \linkS4class{referencebasedParam}, 
#' \linkS4class{independentbulkParam}
#' 
#' @examples
#' # get data
#' lexample <- lute:::.get_decon_example_data_bisque()
#' 
#' # get param object
#' param <- bisqueParam(y.eset = lexample[["y.eset"]], 
#' sc.eset = lexample[["sc.eset"]], batch.variable = "SubjectName",
#' celltype.variable = "cellType")
#' 
#' # get predicted proportions
#' res <- deconvolution(param)
#' 
#' @references 
#' 
#' Brandon Jew and Marcus Alvarez (2021). BisqueRNA: Decomposition of Bulk 
#' Expression with Single-Cell Sequencing. CRAN, R package version 1.0.5.
#' URL: https://CRAN.R-project.org/package=BisqueRNA
#' 
#' Brandon Jew et al. Accurate estimation of cell composition in bulk 
#' expression through robust integration of single-cell information. 
#' Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
#' 
#' @aliases 
#' BisqueParam-class
#'
setClass("bisqueParam", contains="independentbulkParam", 
         slots=c(y.eset = "ExpressionSet", sc.eset = "ExpressionSet", assay.name = "character",
                  batch.variable = "character", celltype.variable = "character"))

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
#' @param assay.name Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batch.variable Name of variable identifying the batches in sc.eset pData/coldata.
#' @param celltype.variable Name of cell type labels variable in sc.eset pData/coldata.
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @returns New object of class \linkS4class{bisqueParam}.
#'
#' @details Takes standard inputs for the Bisque method. If user provides matrices, will convert these
#' into ExpressionSet objects compatible with the main bisque method.
#' 
#' @export
bisqueParam <- function(y = NULL, yi = NULL, z = NULL, s = NULL, 
                        y.eset = NULL, sc.eset = NULL, assay.name = "counts", 
                        batch.variable = "batch.id", 
                        celltype.variable = "celltype", return.info = FALSE) {
  require(Biobase)
  # check y.eset/y
  if(is(y, "NULL")){
    if(is(y.eset, "NULL")){
      stop("Error, need to provide either y or bulk.eset.")
    } else{
      message("Getting y from provided bulk.eset...")
      y <- as.matrix(exprs(y.eset))
    }
  } else{
      if(is(y.eset, "NULL")){
      message("Making ExpressionSet from provided y...")
      y.eset <- .make_eset_from_matrix(mat = y, batch.id = "SubjectName")
      # need at least 2 columns/samples to pass to bisque
      if(ncol(y.eset) == 1){
        sample.name <- colnames(y.eset)
        y.eset <- cbind(y.eset, y.eset)
        colnames(y.eset) <- c(sample.name, paste0(sample.name, "_rep1"))
      }
    }
  }

  # check sc.eset
  if(is(sc.eset, "NULL")){
    stop("Error, no single-cell ExpressionSet provided.")  
    # add condition to call splatter simulations by default?
  } else{
    if(!batch.variable %in% colnames(pData(sc.eset))){
    stop("Error, didn't find batch id variable ",batch.variable,
         " in sc.eset pData/coldata.")
    } else{
      id.sc <- unique(sc.eset[[batch.variable]])
    }
    if(!celltype.variable %in% colnames(pData(sc.eset))){
      stop("Error, didn't find celltype id variable ", celltype.variable, 
           " in sc.eset pData/coldata.")
    }
    if(is(z, "NULL")){
      message("Getting z from sc.eset...")
      sce <- .get_sce_from_eset(sc.eset)
      z <- .get_z_from_sce(sce = sce, celltype.variable = celltype.variable)
    }
  }

  # parse s
  unique.types <- colnames(z)
  unique.types <- unique.types[order(unique.types)]
  if(is(s, "NULL")){
    message("Setting equal cell size factors...")
    s <- rep(1, ncol(z))
    names(s) <- unique.types
  }
  
  # parse batch ids in bulk and sc
  message("Checking batch ids in bulk and sc eset...")
  if(cond <- !batch.variable %in% colnames(pData(y.eset))){
    stop("Error, didn't find batch variable in y.eset pData: ", batch.variable)
  } else{
    id.bulk <- unique(y.eset[[batch.variable]])
  }
  id.overlap <- intersect(id.sc, id.bulk)
  id.unique <- unique(c(id.sc, id.bulk))
  id.onlybulk <- id.bulk[!id.bulk %in% id.overlap]
  id.onlysc <- id.sc[!id.sc %in% id.overlap]
  message("Found ", length(id.unique), " unique batch ids...")
  message("Found ", length(id.overlap), " overlapping batch ids...")
  message("Found ", length(id.onlybulk), " bulk-only batch ids...")
  message("Found ", length(id.onlysc), " sc-only batch ids...")
  if(length(id.overlap) == 0){stop("Error, no overlapping markers in y.eset and sc.eset.")}
  
  # parse independent bulk samples
  if(length(id.onlybulk)==0){
    if(is(yi, "NULL")){
      stop("Error, no independent bulk samples found. ",
        "Provide either yi, or additional y samples.")
    } else{
      message("Using provided yi for independent bulk samples...")
    }
  } else{
    if(is(yi, "NULL")){
      message("Making yi from provided y bulk...")
      yi <- exprs(y.eset)[,colnames(y.eset) %in% id.onlybulk]
    } else{
      message("Using provided yi for independent bulk samples...")
    }
  }

  new("bisqueParam", y = y, yi = yi, z = z, s = s, y.eset = y.eset, 
      sc.eset = sc.eset, assay.name = assay.name, batch.variable = batch.variable, 
      celltype.variable = celltype.variable, return.info = return.info)
}

#' Deconvolution method for bisqueParam
#'
#' Main method to access the Bisque deconvolution method from the main lute 
#' \link{\code{deconvolution}} generic.
#'
#' @details Takes an object of class \linkS4class{bisqueParam} as input, 
#' returning a list.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#' 
#' @references 
#' 
#' Brandon Jew and Marcus Alvarez (2021). BisqueRNA: Decomposition of Bulk 
#' Expression with Single-Cell Sequencing. CRAN, R package version 1.0.5.
#' URL: https://CRAN.R-project.org/package=BisqueRNA
#' 
#' Brandon Jew et al. Accurate estimation of cell composition in bulk 
#' expression through robust integration of single-cell information. 
#' Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
#'
#' @export
setMethod("deconvolution", signature(object = "bisqueParam"), function(object){
  require(BisqueRNA); require(Biobase)
  lparam <- callNextMethod()
  # instantiate objects
  y.eset <- object[["y.eset"]]
  sc.eset <- lparam[["object"]]@sc.eset
  # get predictions
  result <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = y.eset, sc.eset = sc.eset)
  lr <- predictions <- result$bulk.props
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = list(lmd = lparam[["metadata"]], 
                y.eset = y.eset, sc.eset = sc.eset))}
  return(lr)
})
