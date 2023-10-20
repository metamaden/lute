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
#' @seealso \linkS4class{deconvolutionParam}, 
#' \linkS4class{referencebasedParam}, 
#' \linkS4class{independentbulkParam}
#' 
#' @examples
#' ## get data
#' exampleList <- getDeconvolutionExampleData()
#' bulkExpressionSet <- exampleList[["bulkExpressionSet"]][,seq(10)]
#' bulkExpression <- exprs(exampleList[["bulkExpressionSet"]])
#' bulkExpression <- bulkExpression[,c(11:ncol(bulkExpression))]
#' 
#' ## get param object
#' param <- bisqueParam(bulkExpressionSet=bulkExpressionSet, 
#'                      bulkExpression=bulkExpression,
#'                      singleCellExperimentData=
#'                       exampleList[["singleCellExpressionSet"]], 
#'                      batchVariable="SubjectName", 
#'                      cellTypeVariable="cellType", 
#'                      useOverlap=FALSE)
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
setClass("bisqueParam", 
         contains="independentbulkParam", 
         slots=c(
           bulkExpressionSet="ExpressionSet",
           singleCellExpressionSet="ExpressionSet",
           assayName="character",
           batchVariable="character",
           cellTypeVariable="character",
           useOverlap="logical"))

#' Make new object of class bisqueParam
#'
#' Main constructor for class \linkS4class{bisqueParam}.
#'
#' @param bulkExpressionBulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param bulkExpressionIndependent Bulk mixed signals matrix of independent 
#' samples, which should not overlap samples in y.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a provided ExpressionSet containing 
#' single-cell data.
#' @param cellSizeFactor size factor transformations of length equal to the K 
#' cell types to deconvolve.
#' @param bulkExpressionSet ExpressionSet of bulk mixed signals.
#' @param scData SummarizedExperiment-type object of single-cell transcriptomics 
#' data. Accepts ExpressionSet, SummarizedExperiment, and SingleCellExperiment 
#' object types.
#' @param assayName Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batchVariable Name of variable identifbulkExpressionIndependentng the 
#' batches in singleCellExpressionSet pData/coldata.
#' @param cellTypeVariable Name of cell type labels variable in 
#' singleCellExpressionSet pData/coldata.
#' @param useOverlap Whether to deconvolve samples overlapping bulk and sc 
#' esets (logical, FALSE).
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples
#' ## get data
#' exampleList <- getDeconvolutionExampleData()
#' bulkExpressionSet <- exampleList[["bulkExpressionSet"]][,seq(10)]
#' bulkExpression <- exprs(exampleList[["bulkExpressionSet"]])
#' bulkExpression <- bulkExpression[,c(11:ncol(bulkExpression))]
#' 
#' ## get param object
#' param <- bisqueParam(bulkExpressionSet=bulkExpressionSet, 
#'                      bulkExpressionIndependent=bulkExpression,
#'                      scData=exampleList[["singleCellExpressionSet"]], 
#'                      batchVariable="SubjectName", 
#'                      cellTypeVariable="cellType", 
#'                      useOverlap=FALSE)
#' 
#' ## get predicted proportions
#' res <- deconvolution(param)
#'
#' @returns New object of class \linkS4class{bisqueParam}.
#'
#' @details Takes standard inputs for the Bisque method. If user provides 
#' matrices, will convert these into ExpressionSet objects compatible with the 
#' main bisque method.
#' 
#' @export
bisqueParam <- function(bulkExpression=NULL, 
                        bulkExpressionIndependent=NULL, 
                        referenceExpression=NULL, 
                        cellSizeFactor=NULL, 
                        bulkExpressionSet=NULL, 
                        scData=NULL, 
                        assayName="counts", 
                        batchVariable="batch.id", 
                        cellTypeVariable="celltype", 
                        useOverlap=FALSE, 
                        return.info=FALSE) {
  ## check bulkExpressionSet/y
  list.bulkExpression<- 
    .parseBulkExpression(bulkExpression, bulkExpressionSet)
  ## parse scData
  singleCellExpressionSet <- .parseSingleCellData(scData, assayName)
  ## parse z data
  listReferenceExpression <- .parseReferenceExpression(
    singleCellExpressionSet, referenceExpression, assayName, batchVariable, 
    cellTypeVariable)
  ## parse s
  cellSizeFactor <- .parseCellSize(
    listReferenceExpression[["referenceExpression"]], cellSizeFactor)
  ## parse batch ids in bulk and sc
  listBatchID <- .parseBatches(batchVariable=batchVariable,
                                bulkExpressionSet=bulkExpressionSet, 
                                idSC=listReferenceExpression[["idSC"]])
  ## parse independent bulk samples
  bulkExpression <- .parseBulkExpressionIndependent(
    idOnlyBulk=listBatchID[["idOnlyBulk"]], 
    bulkExpression=listBulk[["bulkExpression"]], 
    bulkExpressionIndependent=bulkExpression, 
    bulkExpressionSet=listBulk[["bulkExpressionSet"]])
  
  new("bisqueParam", 
      bulkExpression=bulkExpression, 
      bulkExpressionIndependent=bulkExpressionIndependent, 
      referenceExpression=listReferenceExpression[["referenceExpression"]], 
      cellSizeFactor=cellSizeFactor, 
      bulkExpressionSet=listBulk[["bulkExpressionSet"]], 
      singleCellExpressionSet=singleCellExpressionSet, 
      assayName=assayName, 
      batchVariable=batchVariable, 
      cellTypeVariable=cellTypeVariable, 
      useOverlap=useOverlap, 
      return.info=return.info)
}

#'
.parseBulkExpressionIndependent <- function(idOnlyBulk=NULL, 
                                            bulkExpression=NULL,
                                            bulkExpressionIndependent=NULL, 
                                            bulkExpressionSet=NULL){
  stopOption <- FALSE
  if(length(idOnlyBulk) == 0){
    if(is(bulkExpression, "NULL")){
      stopOption <- TRUE
    } else{}
  } else{
    if(is(bulkExpressionIndependent, "NULL")){
      message(
        "Making bulkExpressionIndependent from provided bulkExpression...")
      filterBulkExpressionIndependent <- 
        colnames(bulkExpressionSet) %in% idOnlyBulk
      bulkExpressionIndependent <- 
        exprs(bulkExpressionSet)[,filterBulkExpressionIndependent]
      colnames(bulkExpressionIndependent) <- 
        colnames(bulkExpressionSet)[filterBulkExpressionIndependent]
      rownames(bulkExpressionIndependent) <- rownames(bulkExpressionSet)
    } else{}
  }
  if(stopOption){stop("Error parsing independent bulk data.")}
  filterBulkSamples<- 
    colnames(bulkExpression) %in% colnames(bulkExpressionIndependent)
  bulkExpression<- bulkExpression[,!filterBulkSamples]
  return(bulkExpression)
}

#'
.parseBatches <- function(
    batchVariable=NULL, bulkExpressionSet=NULL,idSC=NULL){
  stopOption <- FALSE
  message("Checking batch ids in bulk and sc esets...")
  if(batchVariable %in% colnames(pData(bulkExpressionSet))){
    idBulk <- unique(bulkExpressionSet[[batchVariable]])
  } else{
    stopOption <- TRUE
  }
  idOverlap <- intersect(idSC, idBulk)
  idUnique <- unique(c(idSC, idBulk))
  idOnlyBulk <- idBulk[!idBulk %in% idOverlap]
  idOnlySC <-idSC[!idSC %in% idOverlap]
  if(length(idOverlap) == 0){stopOption <- TRUE}
  if(stopOption){stop("Error parsing batches.")}
  return(
    list(idSC=idSC, idBulk=idBulk, idOverlap=idOverlap,
         idUnique=idUnique, idOnlyBulk=idOnlyBulk, idOnlySC=idOnlySC)
  )
}

#'
.parseCellSize <- function(referenceExpression=NULL, cellSizeFactor=NULL){
  uniqueTypes <- colnames(referenceExpression)
  uniqueTypes <- uniqueTypes[order(uniqueTypes)]
  if(is(cellSizeFactor, "NULL")){
    cellSizeFactor <- rep(1, ncol(referenceExpression))
    names(cellSizeFactor) <- uniqueTypes
  }
  return(cellSizeFactor=cellSizeFactor)
}

#'
.parseReferenceExpression <- function(singleCellExpressionSet=NULL, 
                                      referenceExpression=NULL, 
                                      assayName="counts",
                                      batchVariable="group",
                                      cellTypeVariable="celltype"){
  stopOption <- FALSE
  if(!cellTypeVariable %in% colnames(pData(singleCellExpressionSet))){
    stopOption <- TRUE
  }
  if(is(referenceExpression, "NULL")){
    singleCellExperiment <- eset_to_sce(singleCellExpressionSet, "counts")
    referenceExpression <- get_z_from_sce(sce=sce, assayName=assayName, 
                                          cellTypeVariable=cellTypeVariable)
  }
  if(batchVariable %in% colnames(pData(singleCellExpressionSet))){
   idSC <- unique(singleCellExpressionSet[[batchVariable]])
  } else{
    stopOption <- TRUE
  }
  if(stopOption){stop("Error parsing Z data.")}
  return(
    list(singleCellExperiment=singleCellExperiment,
         referenceExpression=referenceExpression, idSC=idSC)
    )
}

#'
.parseSingleCellData <- function(scData=NULL, assayName=assayName){
  stopOption <- FALSE
  if(is(scData, "SingleCellExperiment")){
    singleCellExpressionSet <- sce_to_eset(scData, assayName=assayName)
  } else if(is(scData, "SummarizedExperiment")){
    singleCellExpressionSet <- se_to_eset(scData, assayName=assayName)
  } else if(is(scData, "ExpressionSet")){
    singleCellExpressionSet <- scData
  } else if(is(scData, "NULL")){
    stopOption <- TRUE
  } else{
    stopOption <- TRUE
  }
  if(stopOption){stop("Error parsing sc data.")}
  return(singleCellExpressionSet)
}

#'
.parseBulkExpression<- function(bulkExpression=NULL, bulkExpressionSet=NULL){
  if(is(bulkExpression, "NULL")){
    bulkExpression <- as.matrix(exprs(bulkExpressionSet))
  } else{
    if(is(bulkExpressionSet, "NULL")){
      bulkExpressionSet <- get_eset_from_matrix(
        bulkExpression=bulkExpression, batchVariable="SubjectName")
      ## need at least 2 columns/samples to pass to bisque
      if(ncol(bulkExpressionSet) == 1){
        sampleName <- colnames(bulkExpressionSet)
        bulkExpressionSet <- 
          cbind(bulkExpressionSet, bulkExpressionSet)
        colnames(bulkExpressionSet) <- 
          c(sampleName, paste0(sampleName, "_rep1"))
      }
    }
  }
  return(
    list(bulkExpression=bulkExpression, bulkExpressionSet=bulkExpressionSet))
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
#' exampleList <- getDeconvolutionExampleData()
#' bulkExpressionSet <- exampleList[["bulkExpressionSet"]][,seq(10)]
#' bulkExpression <- exprs(exampleList[["bulkExpressionSet"]])
#' bulkExpression <- bulkExpression[,c(11:ncol(bulkExpression))]
#' 
#' ## get param object
#' param <- bisqueParam(bulkExpressionSet=bulkExpressionSet, 
#'                      bulkExpressionIndependent=bulkExpression,
#'                      scData=exampleList[["singleCellExpressionSet"]], 
#'                      batchVariable="SubjectName", 
#'                      cellTypeVariable="cellType", 
#'                      useOverlap=FALSE)
#' 
#' ## get predicted proportions
#' res <- deconvolution(param)
#'
#' @references Brandon Jew and Marcus Alvarez (2021). BisqueRNA: Decomposition 
#' of Bulk Expression with Single-Cell Sequencing. CRAN, R package version 
#' 1.0.5. URL: https://CRAN.R-project.org/package=BisqueRNA
#' 
#' Brandon Jew et al. Accurate estimation of cell composition in bulk 
#' expression through robust integration of single-cell information. 
#' Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
#'
#' @export
setMethod("deconvolution", signature(object="bisqueParam"), function(object){
  parametersList <- callNextMethod()
  bulkExpressionSet <- object[["bulkExpressionSet"]]
  singleCellExpressionSet <- object[["singleCellExpressionSet"]]
  useOverlap <- object[["useOverlap"]]
  result <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset=bulkExpressionSet, 
    singleCellExpressionSet=singleCellExpressionSet, 
    useOverlap=useOverlap
  )
  predictions <- result$bulk.props
  predictionsList <- lapply(seq(ncol(predictions)), 
                            function(index){predictions[,index]})
  returnList <- .parsedeconvolution_predictions_results(predictionsList, 
                                                 row.names(predictions), 
                                                 colnames(predictions))
  if(object[["return.info"]]){
    returnList <- list(
      predictioncellSizeFactor=predictions, 
      resultInfo=result, 
      metadata=
        list(metadataList=parametersList[["metadata"]], 
             bulkExpressionSet=bulkExpressionSet, 
             singleCellExpressionSet=singleCellExpressionSet))
  }
  return(returnList)
})

#' Show generic behavior for object of class bisqueParam
#' @param object Object of class \linkS4class{bisqueParam} (see 
#' \code{?bisqueParam}).
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' referencebasedParam(
#' bulkExpression=exampleList$bulkExpression, 
#' referenceExpression=exampleList$referenceExpression, 
#' cellScaleFactors=exampleList$cellScaleFactors)
#' 
#' @returns Prints data summary messages to console.
#' @export
setMethod("show", "bisqueParam", function(object) {
  # nnlsParam inherits from deconvolutionParam -> referencebasedParam
  # needs to show standard properties for each parent class
  
  ## deconvolutionParam -- show properties
  bulkExpression <- object[["bulkExpression"]]
  message("Object of class deconvolutionParam")
  message("\nData summaries:")
  message("\tNumber of bulk markers: ", nrow(bulkExpression))
  message("\tNumber of bulk samples: ", ncol(bulkExpression))
  markers <- rownames(bulkExpression)
  if(length(markers) > 10){markers <- markers[seq(10)]}
  message("\tFirst bulk marker labels:\n", 
          paste0(rownames(bulkExpression), collapse="; "))
  samples <- colnames(bulkExpression)
  if(length(samples) > 10){samples <- samples[seq(10)]}
  message("\tFirst sample labels:\n", paste0(samples, collapse="; "), "\n\n")
  
  ## referencebasedParam -- show properties
  ## get metadata
  cellScaleFactors <- object[["cellScaleFactors"]]
  bulkExpression <- object[["bulkExpression"]]
  referenceExpression <- object[["referenceExpression"]]
  uniqueTypes <- try(colnames(object[["referenceExpression"]]))
  markersBulkExpression <- rownames(bulkExpression)
  markersReferenceExpression <- rownames(referenceExpression)
  uniqueMarkers <- unique(
    c(markersBulkExpression, markersReferenceExpression))
  overlappingMarkers <- 
    intersect(markersBulkExpression, markersReferenceExpression)
  markerGenes <- nrow(referenceExpression)
  bulkSamples <- ncol(bulkExpression)
  numberCellTypesK <- ncol(referenceExpression)
  metadataList <- list(
    markerGenes = markerGenes, bulkSamples = bulkSamples, 
    numberCellTypesK = numberCellTypesK, cellScaleFactors = cellScaleFactors, 
    uniqueTypes = uniqueTypes, markersBulkExpression = markersBulkExpression, 
    markersReferenceExpression = markersReferenceExpression)
  ## post console messages
  cat(paste0("class: ", class(object)[1], "\n\n"))
  cat("key deconvolution run info:\n")
  cat("\tmarker info:\n")
  cat("\tsignature markers (Gz): ", markerGenes, "\n")
  cat("\tunique marker labels (Gy | Gz): ", length(uniqueMarkers), "\n")
  cat("\toverlapping marker labels (Gy & Gz): ", 
      length(overlappingMarkers), "\n\n")
  cat("\tsamples info:\n")
  cat("\tnumber of bulk samples (J): ", ncol(object[["bulkExpression"]]), "\n")
  cat("\tsample labels: ", 
      paste0(colnames(bulkExpression), collapse = "; "), "\n")
  cat("\n")
  cat("\tcell size factor properties:\n")
  if(!is(cellScaleFactors, "NULL")){
    for(type in names(cellScaleFactors)){
      cat("\tscale factor for type ", 
          type, ": ", cellScaleFactors["type"], "\n")}
    if(length(cellScaleFactors) == ncol(referenceExpression)){
      referenceExpression <- 
        .zstransform(referenceExpression, cellScaleFactors)}
  }; cat("\n")
  cat("\ttypes info:\n")
  cat("\tnumber of types (K): ", ncol(object[["referenceExpression"]]), "\n")
  if(!(is(uniqueTypes, "NULL")|is(uniqueTypes, "try-error"))){
    uniqueTypes <- uniqueTypes[order(uniqueTypes)]
    cat("\tunique type labels: ", paste0(uniqueTypes, collapse = ";"), "\n")
  } else{
    cat(
      "\nWarning, object 'referenceExpression' has no type labels (colnames)\n")
  }; cat("\n")
  ## parse additional warnings
  if(is(markersBulkExpression, "NULL")){
    cat("Warning, object 'bulkExpression' has no marker labels (rownames)\n\n")}
  if(is(markersReferenceExpression, "NULL")){
    cat(paste0("Warning, object 'referenceExpression' has no marker labels",
               " (rownames)\n\n"))}
  
})