#!/usr/bin/env R

# Author: Sean Maden

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
#' lexample <- get_decon_example_data_bisque()
#' y.eset <- lexample[["y.eset"]]
#' yi <- exprs(y.eset)
#' 
#' # get param object
#' param <- bisqueParam(y.eset=y.eset, yi=yi,
#'                      sc.data=lexample[["sc.eset"]], 
#'                      batch.variable="SubjectName", 
#'                      celltype.variable="cellType", 
#'                      use.overlap=FALSE)
#' 
#' # get predicted proportions
#' res <- deconvolution(param)
#' 
#' @references Brandon Jew and Marcus Alvarez (2021). BisqueRNA: Decomposition of Bulk 
#' Expression with Single-Cell Sequencing. CRAN, R package version 1.0.5.
#' URL: https://CRAN.R-project.org/package=BisqueRNA
#' 
#' Brandon Jew et al. Accurate estimation of cell composition in bulk 
#' expression through robust integration of single-cell information. 
#' Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
#'
#' @returns New object of class \linkS4class{bisqueParam}.
#'
#' @aliases 
#' BisqueParam-class
#'
setClass("bisqueParam", contains="independentbulkParam", 
         slots=c(y.eset="ExpressionSet", sc.eset="ExpressionSet", 
                 assay.name="character", batch.variable="character", 
                 celltype.variable="character", use.overlap="logical"))

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
#' @param sc.data SummarizedExperiment-type object of single-cell transcriptomics data. Accepts
#' ExpressionSet, SummarizedExperiment, and SingleCellExperiment object types.
#' @param assay.name Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batch.variable Name of variable identifying the batches in sc.eset pData/coldata.
#' @param celltype.variable Name of cell type labels variable in sc.eset pData/coldata.
#' @param use.overlap Whether to deconvolve samples overlapping bulk and sc 
#' esets (logical, FALSE).
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#' 
#' @examples
#' # get data
#' lexample <- get_decon_example_data_bisque()
#' y.eset <- lexample[["y.eset"]]
#' yi <- exprs(y.eset)
#' sc.data <- lexample[["sc.eset"]]
#' # get param object
#' param <- bisqueParam(
#' y.eset=y.eset, yi=yi, sc.data=sc.data, batch.variable="SubjectName", 
#' celltype.variable="cellType", use.overlap=FALSE
#' )
#' # get predicted proportions
#' res <- deconvolution(param)
#'
#' @returns New object of class \linkS4class{bisqueParam}.
#'
#' @details Takes standard inputs for the Bisque method. If user provides matrices, will convert these
#' into ExpressionSet objects compatible with the main bisque method.
#' 
#' @export
bisqueParam <- function(y=NULL, yi=NULL, z=NULL, s=NULL, 
                        y.eset=NULL, sc.data=NULL, assay.name="counts", 
                        batch.variable="batch.id", 
                        celltype.variable="celltype", 
                        use.overlap=FALSE, return.info=FALSE) {
  # check y.eset/y
  list.y <- .parse_y(y, y.eset)
  # parse sc.data
  sc.eset <- .parse_sc(sc.data, assay.name)
  # parse z data
  list.z <- .parse_z(sc.eset, z, assay.name, batch.variable, 
                     celltype.variable)
  # parse s
  s <- .parse_s(list.z[["z"]], s)
  # parse batch ids in bulk and sc
  list.batchid <- .parse_batches(batch.variable=batch.variable,
                                 y.eset=y.eset, id.sc=list.z[["id.sc"]])
  # parse independent bulk samples
  y <- .parse_independent_bulk(id.onlybulk=list.batchid[["id.onlybulk"]], 
                               y=list.y[["y"]], yi=yi, 
                               y.eset=list.y[["y.eset"]])
  new("bisqueParam", y=y, yi=yi, z=list.z[["z"]], s=s, 
      y.eset=list.y[["y.eset"]], sc.eset=sc.eset, 
      assay.name=assay.name, batch.variable=batch.variable, 
      celltype.variable=celltype.variable, 
      use.overlap=use.overlap, return.info=return.info)
}

#'
.parse_independent_bulk <- function(id.onlybulk=NULL, y=NULL,
                                    yi=NULL, y.eset=NULL){
  stop.option <- FALSE
  if(length(id.onlybulk)==0){
    if(is(yi, "NULL")){
      stop.option <- TRUE
    } else{}
  } else{
    if(is(yi, "NULL")){
      message("Making yi from provided y bulk...")
      filter.yi <- colnames(y.eset) %in% id.onlybulk
      yi <- exprs(y.eset)[,filter.yi]
      colnames(yi) <- colnames(y.eset)[filter.yi]
      rownames(yi) <- rownames(y.eset)
    } else{}
  }
  if(stop.option){stop("Error parsing independent bulk data.")}
  filter.bulk.samples.y <- colnames(y) %in% colnames(yi)
  y <- y[,!filter.bulk.samples.y]
  return(y)
}

#'
.parse_batches <- function(batch.variable=NULL, y.eset=NULL,
                           id.sc=NULL){
  stop.option <- FALSE
  message("Checking batch ids in bulk and sc eset...")
  if(batch.variable %in% colnames(pData(y.eset))){
    id.bulk <- unique(y.eset[[batch.variable]])
  } else{
    stop.option <- TRUE
  }
  id.overlap <- intersect(id.sc, id.bulk)
  id.unique <- unique(c(id.sc, id.bulk))
  id.onlybulk <- id.bulk[!id.bulk %in% id.overlap]
  id.onlysc <- id.sc[!id.sc %in% id.overlap]
  if(length(id.overlap) == 0){stop.option <- TRUE}
  if(stop.option){stop("Error parsing batches.")}
  return(list(id.sc=id.sc, 
              id.bulk=id.bulk, 
              id.overlap=id.overlap, 
              id.unique=id.unique, 
              id.onlybulk=id.onlybulk, 
              id.onlysc=id.onlysc))
}

#'
.parse_s <- function(z=NULL, s=NULL){
  unique.types <- colnames(z)
  unique.types <- unique.types[order(unique.types)]
  if(is(s, "NULL")){
    s <- rep(1, ncol(z)); names(s) <- unique.types
  }
  return(s=s)
}

#'
.parse_z <- function(sc.eset=NULL, z=NULL,
                     assay.name="counts",
                     batch.variable="group",
                     celltype.variable="celltype"){
  stop.option <- FALSE
  if(!celltype.variable %in% colnames(pData(sc.eset))){
    stop.option <- TRUE
  }
  if(is(z, "NULL")){
    sce <- eset_to_sce(sc.eset, "counts")
    z <- get_z_from_sce(sce=sce, 
                        assay.name=assay.name, 
                        celltype.variable=celltype.variable)
  }
  if(batch.variable %in% colnames(pData(sc.eset))){
    id.sc <- unique(sc.eset[[batch.variable]])
  } else{
    stop.option <- TRUE
  }
  if(stop.option){stop("Error parsing Z data.")}
  return(list(sce=sce, z=z, id.sc=id.sc))
}

#'
.parse_sc <- function(sc.data=NULL, assay.name=assay.name){
  stop.option <- FALSE
  if(is(sc.data, "SingleCellExperiment")){
    sc.eset <- sce_to_eset(sc.data, assay.name=assay.name)
  } else if(is(sc.data, "SummarizedExperiment")){
    sc.eset <- se_to_eset(sc.data, assay.name=assay.name)
  } else if(is(sc.data, "ExpressionSet")){
    sc.eset <- sc.data
  } else if(is(sc.data, "NULL")){
    stop.option <- TRUE
  } else{
    stop.option <- TRUE
  }
  if(stop.option){stop("Error parsing sc data.")}
  return(sc.eset)
}

#'
.parse_y <- function(y=NULL, y.eset=NULL){
  if(is(y, "NULL")){
    y <- as.matrix(exprs(y.eset))
  } else{
    if(is(y.eset, "NULL")){
      y.eset <- get_eset_from_matrix(mat=y, batch.variable="SubjectName")
      # need at least 2 columns/samples to pass to bisque
      if(ncol(y.eset) == 1){
        sample.name <- colnames(y.eset)
        y.eset <- cbind(y.eset, y.eset)
        colnames(y.eset) <- c(sample.name, paste0(sample.name, "_rep1"))
      }
    }
  }
  return(list(y=y, y.eset=y.eset))
}

#' Deconvolution method for bisqueParam
#'
#' Main method to access the Bisque deconvolution method from the main lute 
#' \code{deconvolution} generic.
#'
#' @param object Object of type \linkS4class{bisqueParam}.
#' @details Takes an object of class \linkS4class{bisqueParam} as input, 
#' returning a list.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#' 
#' @examples
#' # get data
#' lexample <- get_decon_example_data_bisque()
#' y.eset <- lexample[["y.eset"]]
#' yi <- exprs(y.eset)
#' 
#' # get param object
#' param <- bisqueParam(y.eset=y.eset, yi=yi,
#'                      sc.data=lexample[["sc.eset"]], 
#'                      batch.variable="SubjectName", 
#'                      celltype.variable="cellType", 
#'                      use.overlap=FALSE)
#' 
#' # get predicted proportions
#' res <- deconvolution(param)
#' 
#' @references Brandon Jew and Marcus Alvarez (2021). BisqueRNA: Decomposition of Bulk 
#' Expression with Single-Cell Sequencing. CRAN, R package version 1.0.5.
#' URL: https://CRAN.R-project.org/package=BisqueRNA
#' 
#' Brandon Jew et al. Accurate estimation of cell composition in bulk 
#' expression through robust integration of single-cell information. 
#' Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
#'
#' @export
setMethod("deconvolution", signature(object="bisqueParam"), function(object){
  lparam <- callNextMethod()
  y.eset <- object[["y.eset"]]
  sc.eset <- object[["sc.eset"]]
  use.overlap <- object[["use.overlap"]]
  result <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset=y.eset, 
                                                   sc.eset=sc.eset, 
                                                   use.overlap=use.overlap)
  predictions <- result$bulk.props
  lpred <- lapply(seq(ncol(predictions)), function(index){predictions[,index]})
  lr <- .parse_deconvolution_predictions_results(lpred, 
                                                 row.names(predictions), 
                                                 colnames(predictions))
  if(object[["return.info"]]){
    lr <- list(predictions=predictions, result.info=result, 
               metadata=list(lmd=lparam[["metadata"]], 
                y.eset=y.eset, sc.eset=sc.eset))}
  return(lr)
})

#' Show generic behavior for object of class \linkS4class{bisqueParam}
#' @param object An object of class \linkS4class{bisqueParam}.
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{bisqueParam}
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @returns Shows object summaries.
#' 
#' @export
setMethod("show", signature(object="bisqueParam"), function(object){
  show(object)
})