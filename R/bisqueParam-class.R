#!/usr/bin/env R

### Author: Sean Maden

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
#' ## get data
#' exampleList <- get_decon_example_data_bisque()
#' bulkExpressionSet <- exampleList[["bulkExpressionSet"]][,seq(10)]
#' bulkExpression <- exprs(exampleList[["bulkExpressionSet"]])
#' bulkExpression <- bulkExpression[,c(11:ncol(bulkExpression))]
#' 
#' ## get param object
#' param <- bisqueParam(bulkExpressionSet=bulkExpressionSet, bulkExpression=bulkExpression,
#'                      singleCellExperimentData=exampleList[["singleCellExpressionSet"]], 
#'                      batch.variable="SubjectName", 
#'                      celltype.variable="cellType", 
#'                      use.overlap=FALSE)
#' 
#' ## get predicted proportions
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
         slots=c(bulkExpressionSet="ExpressionSet", 
                 singleCellExpressionSet="ExpressionSet", 
                 assayName="character", 
                 batch.variable="character", 
                 celltype.variable="character", 
                 use.overlap="logical"))

#' Make new object of class bisqueParam
#'
#' Main constructor for class \linkS4class{bisqueParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param yi Bulk mixed signals matrix of independent samples, which should not overlap samples in y.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param bulkExpressionSet ExpressionSet of bulk mixed signals.
#' @param sc.data SummarizedExperiment-type object of single-cell transcriptomics data. Accepts
#' ExpressionSet, SummarizedExperiment, and SingleCellExperiment object types.
#' @param assayName Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batch.variable Name of variable identifying the batches in singleCellExpressionSet pData/coldata.
#' @param celltype.variable Name of cell type labels variable in singleCellExpressionSet pData/coldata.
#' @param use.overlap Whether to deconvolve samples overlapping bulk and sc 
#' esets (logical, FALSE).
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#' 
#' @examples
#' ## get data
#' lexample <- get_decon_example_data_bisque()
#' input_y_eset <- lexample[["bulkExpressionSet"]][,seq(10)]
#' input_yi <- exprs(lexample[["bulkExpressionSet"]])
#' input_yi <- input_yi[,c(11:ncol(input_yi))]
#' 
#' ## get param object
#' param <- bisqueParam(bulkExpressionSet=input_y_eset, yi=input_yi,
#'                      sc.data=lexample[["singleCellExpressionSet"]], 
#'                      batch.variable="SubjectName", 
#'                      celltype.variable="cellType", 
#'                      use.overlap=FALSE)
#' 
#' ## get predicted proportions
#' res <- deconvolution(param)
#'
#' @returns New object of class \linkS4class{bisqueParam}.
#'
#' @details Takes standard inputs for the Bisque method. If user provides matrices, will convert these
#' into ExpressionSet objects compatible with the main bisque method.
#' 
#' @export
bisqueParam <- function(y=NULL, yi=NULL, z=NULL, s=NULL, 
                        bulkExpressionSet=NULL, sc.data=NULL, assayName="counts", 
                        batch.variable="batch.id", 
                        celltype.variable="celltype", 
                        use.overlap=FALSE, return.info=FALSE) {
  input_y <- y; input_yi <- yi; input_z <- z; input_s <- s; 
  input_y_eset <- bulkExpressionSet
  ## check bulkExpressionSet/y
  list.y <- .parse_y(input_y, input_y_eset)
  ## parse sc.data
  singleCellExpressionSet <- .parse_sc(sc.data, assayName)
  ## parse z data
  list.z <- .parse_z(singleCellExpressionSet, input_z, assayName, batch.variable, 
                     celltype.variable)
  ## parse s
  input_s <- .parse_s(list.z[["z"]], input_s)
  ## parse batch ids in bulk and sc
  list.batchid <- .parse_batches(batch.variable=batch.variable,
                                 bulkExpressionSet=bulkExpressionSet, id.sc=list.z[["id.sc"]])
  ## parse independent bulk samples
  input_y <- .parse_independent_bulk(
    id.onlybulk=list.batchid[["id.onlybulk"]], y=list.y[["y"]], 
    yi=input_yi, bulkExpressionSet=list.y[["bulkExpressionSet"]])
  
  new("bisqueParam", y=input_y, yi=input_yi, z=list.z[["z"]], s=input_s, 
      bulkExpressionSet=list.y[["bulkExpressionSet"]], singleCellExpressionSet=singleCellExpressionSet, 
      assayName=assayName, batch.variable=batch.variable, 
      celltype.variable=celltype.variable, 
      use.overlap=use.overlap, return.info=return.info)
}

#'
.parse_independent_bulk <- function(id.onlybulk=NULL, y=NULL,
                                    yi=NULL, bulkExpressionSet=NULL){
  input_y <- y; input_yi <- yi; input_y_eset <- bulkExpressionSet
  stop.option <- FALSE
  if(length(id.onlybulk) == 0){
    if(is(input_yi, "NULL")){
      stop.option <- TRUE
    } else{}
  } else{
    if(is(input_yi, "NULL")){
      message("Making yi from provided y bulk...")
      filter.yi <- colnames(input_y_eset) %in% id.onlybulk
      input_yi <- exprs(input_y_eset)[,filter.yi]
      colnames(input_yi) <- colnames(input_y_eset)[filter.yi]
      rownames(input_yi) <- rownames(input_y_eset)
    } else{}
  }
  if(stop.option){stop("Error parsing independent bulk data.")}
  filter.bulk.samples.y <- colnames(input_y) %in% colnames(input_yi)
  input_y <- input_y[,!filter.bulk.samples.y]
  return(input_y)
}

#'
.parse_batches <- function(batch.variable=NULL, bulkExpressionSet=NULL,
                           id.sc=NULL){
  stop.option <- FALSE; input_y_eset <- bulkExpressionSet
  message("Checking batch ids in bulk and sc eset...")
  if(batch.variable %in% colnames(pData(input_y_eset))){
    id.bulk <- unique(input_y_eset[[batch.variable]])
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
  input_z <- z; input_s <- s
  unique.types <- colnames(input_z)
  unique.types <- unique.types[order(unique.types)]
  if(is(input_s, "NULL")){
    input_s <- rep(1, ncol(input_z)); names(input_s) <- unique.types
  }
  return(s=input_s)
}

#'
.parse_z <- function(singleCellExpressionSet=NULL, z=NULL,
                     assayName="counts",
                     batch.variable="group",
                     celltype.variable="celltype"){
  stop.option <- FALSE; input_z <- z
  if(!celltype.variable %in% colnames(pData(singleCellExpressionSet))){
    stop.option <- TRUE
  }
  if(is(z, "NULL")){
    sce <- eset_to_sce(singleCellExpressionSet, "counts")
    input_z <- get_z_from_sce(sce=sce, 
                        assayName=assayName, 
                        celltype.variable=celltype.variable)
  }
  if(batch.variable %in% colnames(pData(singleCellExpressionSet))){
    id.sc <- unique(singleCellExpressionSet[[batch.variable]])
  } else{
    stop.option <- TRUE
  }
  if(stop.option){stop("Error parsing Z data.")}
  return(list(sce=sce, z=input_z, id.sc=id.sc))
}

#'
.parse_sc <- function(sc.data=NULL, assayName=assayName){
  stop.option <- FALSE
  if(is(sc.data, "SingleCellExperiment")){
    singleCellExpressionSet <- sce_to_eset(sc.data, assayName=assayName)
  } else if(is(sc.data, "SummarizedExperiment")){
    singleCellExpressionSet <- se_to_eset(sc.data, assayName=assayName)
  } else if(is(sc.data, "ExpressionSet")){
    singleCellExpressionSet <- sc.data
  } else if(is(sc.data, "NULL")){
    stop.option <- TRUE
  } else{
    stop.option <- TRUE
  }
  if(stop.option){stop("Error parsing sc data.")}
  return(singleCellExpressionSet)
}

#'
.parse_y <- function(y=NULL, bulkExpressionSet=NULL){
  input_y <- y; input_y_eset <- bulkExpressionSet
  if(is(input_y, "NULL")){
    input_y <- as.matrix(exprs(input_y_eset))
  } else{
    if(is(input_y_eset, "NULL")){
      input_y_eset <- get_eset_from_matrix(
        mat=input_y, batch.variable="SubjectName")
      ## need at least 2 columns/samples to pass to bisque
      if(ncol(input_y_eset) == 1){
        sample.name <- colnames(input_y_eset)
        input_y_eset <- cbind(input_y_eset, input_y_eset)
        colnames(input_y_eset) <- c(sample.name, paste0(sample.name, "_rep1"))
      }
    }
  }
  return(list(y=input_y, bulkExpressionSet=input_y_eset))
}

#' Deconvolution method for bisqueParam
#'
#' Main method to access the Bisque deconvolution method from the main lute 
#' \code{deconvolution} generic.
#'
#' @param object Object of type \linkS4class{bisqueParam} (see 
#' \code{?bisqueParam}).
#' @details Takes an object of class \linkS4class{bisqueParam} as input, 
#' returning a list.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#' 
#' @examples
#' ## get data
#' lexample <- get_decon_example_data_bisque()
#' input_y_eset <- lexample[["bulkExpressionSet"]][,seq(10)]
#' input_yi <- exprs(lexample[["bulkExpressionSet"]])
#' input_yi <- input_yi[,c(11:ncol(input_yi))]
#' 
#' ## get param object
#' param <- bisqueParam(bulkExpressionSet=input_y_eset, yi=input_yi,
#'                      sc.data=lexample[["singleCellExpressionSet"]], 
#'                      batch.variable="SubjectName", 
#'                      celltype.variable="cellType", 
#'                      use.overlap=FALSE)
#' 
#' ## get predicted proportions
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
  input_y_eset <- object[["bulkExpressionSet"]]
  singleCellExpressionSet <- object[["singleCellExpressionSet"]]
  use.overlap <- object[["use.overlap"]]
  result <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset=input_y_eset, 
                                                   singleCellExpressionSet=singleCellExpressionSet, 
                                                   use.overlap=use.overlap)
  predictions <- result$bulk.props
  lpred <- lapply(seq(ncol(predictions)), function(index){predictions[,index]})
  return_list <- .parse_deconvolution_predictions_results(lpred, 
                                                 row.names(predictions), 
                                                 colnames(predictions))
  if(object[["return.info"]]){
    return_list <- list(
      predictions=predictions, result.info=result, 
      metadata=
        list(lmd=lparam[["metadata"]], bulkExpressionSet=input_y_eset, singleCellExpressionSet=singleCellExpressionSet))}
  return(return_list)
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