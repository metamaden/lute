#!/usr/bin/env R

### Author: Sean Maden
###
### Utilities and misingleCellExperimentllaneous functions supporting the lute package for deconvolution experiments.
###

#' get_celltypes_from_singleCellExperiment
#' 
#' Extract cell type values from SingleCellExperiment.
#' 
#' @param singleCellExperiment A SingleCellExperiment object.
#' @param cellTypeVariable Variable containing cell type labels (e.g. "type1", 
#' "type2", etc.).
#' @returns List of cell type variable metadata and values.
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame pData exprs
#' @importFrom SingleCellExperiment counts SingleCellExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' 
#' @examples
#' lexample <- getDeconvolutionExampleData()
#'
#' @export
get_celltypes_from_singleCellExperiment <- function(singleCellExperiment, cellTypeVariable="celltype"){
  celltype.vector <- as.data.frame(SummarizedExperiment::colData(singleCellExperiment))[,cellTypeVariable]
  celltype.char <- as.character(celltype.vector)
  unique.types <- unique(celltype.char)
  unique.types <- unique.types[order(unique.types)]
  celltype.fact <- factor(celltype.vector, levels=unique.types)
  lr <- list(variable=cellTypeVariable, 
             unique.types=unique.types, 
             character=celltype.char, 
             factor=celltype.fact)
  return(lr)
}

#' ypb_from_singleCellExperiment
#'
#' Get pseudobulk from a SingleCellExperiment object.
#'
#' @param singleCellExperiment An object of type \linkS4class{SingleCellExperiment}.
#' @param assayName Name of expression matrix in \code{singleCellExperiment} assays.
#' @param cellTypeVariable Variable name for cell type labels in \code{singleCellExperiment} 
#' coldata.
#' @param sample.id.variable Variable name for sample/group ID labels in 
#' \code{singleCellExperiment} coldata.
#' @param S Vector of cell type size scale factors. Optional.
#' @returns Matrix of simulated bulk convoluted signals.
#' 
#' @examples
#' singleCellExperiment.example <- random_singleCellExperiment()
#' ypb_from_singleCellExperiment(singleCellExperiment.example)
#' 
#' @export
ypb_from_singleCellExperiment <- function(singleCellExperiment, assayName="counts", 
                         cellTypeVariable="celltype", 
                         sample.id.variable=NULL, S=NULL){
  num.groups <- 1; unique.group.id.vector <- ""
  if(!is(sample.id.variable, "NULL")){
    group.id.vector <- singleCellExperiment[[sample.id.variable]]
    unique.group.id.vector <- group.id.vector 
    unique.group.id.vector <- unique(unique.group.id.vector)
    unique.group.id.vector <- as.character(unique.group.id.vector)
    num.groups <- length(unique.group.id.vector)
  }
  list.cell.types <- get_celltypes_from_singleCellExperiment(
    singleCellExperiment=singleCellExperiment, cellTypeVariable=cellTypeVariable)
  numberType <- length(list.cell.types[["unique.types"]])
  ypb.list <- lapply(unique.group.id.vector, function(group.id){
    singleCellExperiment.filter <- singleCellExperiment
    if(num.groups > 1){
      filter.group <- singleCellExperiment[[sample.id.variable]]==group.id
      singleCellExperiment.filter <- singleCellExperiment[,filter.group]
    }
    
    cellScaleFactors <- S
    if(is(cellScaleFactors, "NULL")){
      cellScaleFactors <- rep(1, numberType)
      names(cellScaleFactors) <- list.cell.types[["unique.types"]]
    }
    
    referenceExpressionZnew <- get_z_from_singleCellExperiment(singleCellExperiment.filter, assayName, cellTypeVariable)
    input_P <- table(list.cell.types[["character"]])
    input_P <- prop.table(input_P)
    order.p <- match(names(input_P), list.cell.types[["unique.types"]])
    order.p <- order(order.p)
    input_P <- input_P[order.p]
    referenceExpressionZSnew <- .zstransform(referenceExpressionZnew, cellScaleFactors)
    bulkExpressionPseudobulk <- t(t(input_P) %*% t(referenceExpressionZSnew))
    return(bulkExpressionPseudobulk)
  })
  ypb.table <- do.call(cbind, ypb.list)
  ypb.table <- as.data.frame(ypb.table)
  if(num.groups > 1){
    colnames(ypb.table) <- unique.group.id.vector
  } else{
    colnames(ypb.table) <- "singleCellExperiment.pseudobulk"
  }
  return(ypb.table)
}

#' signature_matrix_from_singleCellExperiment
#' 
#' Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' 
#' @param singleCellExperiment An object of type \linkS4class{SingleCellExperiment}.
#' @param assayName Name of expression matrix in \code{singleCellExperiment} assays (e.g. 
#' "counts").
#' @param cellTypeVariable Variable name for cell type labels in \code{singleCellExperiment} 
#' coldata (e.g. "type1", "type2", etc.). 
#' @param summary.method Summary statistic function to use.
#' @details Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' @returns New Z signature matrix.
#' 
#' @examples
#' singleCellExperiment.example <- random_singleCellExperiment()
#' signature_matrix_from_singleCellExperiment(singleCellExperiment.example)
#' 
#' @export
signature_matrix_from_singleCellExperiment <- function(singleCellExperiment, 
                                      cellTypeVariable="celltype", 
                                      summary.method="mean", 
                                      assayName="counts"){
  ## gets the z signature matrix from an singleCellExperiment object
  expression.matrix <- assays(singleCellExperiment)[[assayName]]
  expression.matrix <- as.matrix(expression.matrix)
  cd <- colData(singleCellExperiment)
  unique.cell.types <- unique(cd[,cellTypeVariable])
  unique.cell.types <- unique.cell.types[order(unique.cell.types)]
  referenceExpression <- do.call(cbind, lapply(unique.cell.types, function(cell.type.index){
    filter.index <- cd[,cellTypeVariable]==cell.type.index
    if(summary.method == "mean"){
      DelayedArray::rowMeans(expression.matrix[,filter.index])
    } else{
      Biobase::rowMedians(expression.matrix[,filter.index])
    }
  }))
  colnames(referenceExpression) <- unique.cell.types
  return(referenceExpression)
}

#' get_z_from_singleCellExperiment
#' 
#' Makes the Z cell atlas reference from a SingleCellExperiment.
#' 
#' @param singleCellExperiment A SingleCellExperiment object.
#' @param assayName Name of expression assay type (e.g. "counts").
#' @param cellTypeVariable Name of variable containing cell type labels (e.g. 
#' "type1", "type2", etc.).
#' @returns Matrix of cell summary values (Z reference atlas).
#'
#' @importFrom SummarizedExperiment assays
#' 
#' @examples
#' lexample <- getDeconvolutionExampleData()
#' @export
referenceFromSingleCellExperiment <- function(
    singleCellExperiment, assayName="counts", cellTypeVariable="celltype"){
  typesList<- get_celltypes_from_singleCellExperiment(
    singleCellExperiment=singleCellExperiment, 
    cellTypeVariable=cellTypeVariable)
  mexpr <- as.matrix(
    assays(
      singleCellExperiment)[[assayName]])
  referenceExpressionZnew <- 
    do.call(cbind, lapply(typesList[["unique.types"]], function(typei){
    datav <- mexpr[,typesList[["character"]]==typei]; .z_operator(datav)
  }))
  colnames(referenceExpressionZnew) <- typesList[["unique.types"]]
  rownames(referenceExpressionZnew) <- rownames(singleCellExperiment)
  return(referenceExpressionZnew)
}

.referenceExpressionCellScaleTransform <- function(
    referenceExpression, cellScaleFactor){
  sweep(referenceExpression, 2, cellScaleFactor, FUN="*")
}

.referenceExpressionOperator <- function(datav){
  rowMeans(datav)
}

#' getDeconvolutionExampleData
#' 
#' Make example data for deconvolution.
#' 
#' @param numberBulkSamples Number of bulk samples.
#' @param numberMarkers Number of cell type markers.
#' @param numberType Number of cell types.
#' @returns Example data as list.
#' 
#' @importFrom stats rpois
#' 
#' @examples
#' exampleData <- getDeconvolutionExampleData()
#' 
#' @export
getDeconvolutionExampleData <- function(
    numberBulkSamples=2, numberMarkers=10, numberTypes=2){
  bulkExpression <- matrix(
    rpois(n=numberMarkers*numberBulkSamples, lambda=seq(0, 50, 5)), 
    ncol=numberBulkSamples)
  referenceExpression <- matrix(
    rpois(n=numberType*numberMarkers, lambda=seq(0, 50, 5)), 
    ncol=numberType)
  rownames(bulkExpression) <- rownames(referenceExpression) <- 
    paste0("marker", seq(numberMarkers))
  colnames(referenceExpression) <- paste0("type", seq(numberType))
  colnames(bulkExpression) <- paste0("sample", seq(numberBulkSamples))
  cellScaleFactors <- c(1, 10)
  names(cellScaleFactors) <- colnames(referenceExpression)
  return(list(
    referenceExpression=referenceExpression, 
    bulkExpression=bulkExpression, 
    cellScaleFactors=cellScaleFactors))
}

.getReferenceExpressionVariable <- function(referenceExpression){
  referenceExpressionVariable <- 
    matrix(0, nrow=nrow(referenceExpression), ncol=ncol(referenceExpression))
  rownames(referenceExpressionVariable) <- rownames(referenceExpression)
  colnames(referenceExpressionVariable) <- colnames(referenceExpression)
  referenceExpressionVariable
}

#' getDeconvolutionExampleData_bisque
#'
#' Get example data for Bisque algorithm.
#'
#' @param numberBulkSamples Number of bulk samples.
#' @param numberMarkers Number of cell type markers.
#' @param numberCells Number of cells.
#' @param numberType Number of cell types.
#' @returns Example data as list.
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#' 
#' @examples
#' exampleData <- getDeconvolutionExampleDataBisque()
#'
#' @export
getDeconvolutionExampleDataBisque <- function(numberBulkSamples=100,
                                               numberMarkers=1000, 
                                               numberCells=1000, 
                                               numberType=2){
  exampleList <- getDeconvolutionExampleData(
    numberBulkSamples=numberBulkSamples,
    numberMarkers=numberMarkers,
    numberType=numberType)
  bulkExpression <- exampleList[["bulkExpression"]]
  colnames(bulkExpression) <- c(paste0("sample", seq(numberBulkSamples/2)), 
                   paste0("bulk", seq(numberBulkSamples/2)))
  dfBulkPheno <- data.frame(SubjectName=colnames(bulkExpression))
  rownames(dfBulkPheno) <- colnames(bulkExpression)
  bulkExpressionSet <- ExpressionSet(assayData=bulkExpression, 
                          phenoData=AnnotatedDataFrame(dfBulkPheno))
  singleCellExperiment <- random_singleCellExperiment(num.genes=numberMarkers, 
                    numberCells=numberCells, 
                    numberType=numberType)
  dfReferenceExpressionPheno <- data.frame(
    cellType=singleCellExperiment[["celltype"]], 
                           SubjectName=
                             paste0("sample", seq(numberCells)))
  rownames(dfReferenceExpressionPheno) <- colnames(singleCellExperiment)
  referenceExpressionSet <- ExpressionSet(
    assayData=counts(singleCellExperiment), 
    phenoData=AnnotatedDataFrame(dfReferenceExpressionPheno))
  rownames(referenceExpressionSet) <- rownames(bulkExpressionSet)
  returnList <- list(bulkExpressionSet=bulkExpressionSet,
                     singleCellExpressionSet=referenceExpressionSet)
  return(returnList)
}

#' getDeconvolutionExampleDataSCDC
#' 
#' Get example data for SCDC
#'
#' @returns Example data as list.
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#' 
#' @examples
#' exampleData <- getDeconvolutionExampleDataSCDC()
#' 
#' @export
getDeconvolutionExampleDataSCDC <- function(){
  ## get bulkExpressionSet
  bulkExpression <- getDeconvolutionExampleData()[["y"]]
  bulkExpression <- cbind(bulkExpression, bulkExpression, bulkExpression, 
                          bulkExpression, bulkExpression, bulkExpression)
  colnames(bulkExpression) <- c(paste0("sample", seq(2)), paste0("bulk",seq(4)))
  dfBulkPheno <- data.frame(SubjectName=colnames(bulkExpression))
  rownames(dfBulkPheno) <- colnames(bulkExpression)
  bulkExpressionSet <- ExpressionSet(assayData=bulkExpression, 
                          phenoData=AnnotatedDataFrame(dfBulkPheno))

  ## get referenceExpressionSet
  singleCellExperiment <- random_singleCellExperiment(num.genes=10, 
                                                      numberCells=300, 
                                                      numberType=4)
  dfReferenceExpressionPheno <- data.frame(
    cellType=singleCellExperiment[["celltype"]], 
                           SubjectName=
      paste0("sample", seq(ncol(singleCellExperiment))))
  rownames(dfReferenceExpressionPheno) <- colnames(singleCellExperiment)
  referenceExpressionSet <- ExpressionSet(
    assayData=counts(singleCellExperiment),
    phenoData=AnnotatedDataFrame(dfReferenceExpressionPheno))
  rownames(referenceExpressionSet) <- rownames(bulkExpressionSet)
  
  ## return
  returnList <- list(
    bulkExpressionSet=bulkExpressionSet, 
    singleCellExpressionSet=referenceExpressionSet)
  return(returnList)
}

#'
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#'
getDeconvolutionExampleData_music2 <- function(){
  ## get bulkExpressionSet
  bulkExpression <- getDeconvolutionExampleData()[["y"]]
  bulkExpression <- cbind(bulkExpression, bulkExpression, bulkExpression, 
                          bulkExpression, bulkExpression, bulkExpression)
  colnames(bulkExpression) <- c(paste0("sample", seq(2)), 
                   paste0("bulk",seq(ncol(bulkExpression)-2)))
  
  dfBulkPheno <- data.frame(SubjectName=colnames(bulkExpression))
  rownames(dfBulkPheno) <- colnames(bulkExpression)
  bulkExpressionSet <- ExpressionSet(assayData=bulkExpression, 
                          phenoData=AnnotatedDataFrame(dfBulkPheno))

  ## get referenceExpressionSet
  singleCellExperiment <- random_singleCellExperiment(num.genes=10, 
                                                      numberCells=300, 
                                                      numberType=2)
  dfReferenceExpressionPheno <- data.frame(
    cellType=singleCellExperiment[["celltype"]], 
                           SubjectName=paste0("sample", 
                                              seq(ncol(singleCellExperiment))))
  rownames(dfReferenceExpressionPheno) <- colnames(singleCellExperiment)
  referenceExpressionSet <- ExpressionSet(
    assayData=counts(singleCellExperiment),
    phenoData=AnnotatedDataFrame(dfReferenceExpressionPheno))
  rownames(referenceExpressionSet) <- rownames(bulkExpressionSet)
  
  ## return
  returnList <- list(bulkExpressionSet=bulkExpressionSet, 
                      singleCellExpressionSet=referenceExpressionSet)
  return(returnList)
}


#' getDeconvolutionExampleDataSCDC
#' 
#' Get example data for SCDC
#'

#' parseDeconvolutionPredictionsResults
#' 
#' Gets formatted predicted cell type proportions table from deconvolution 
#' results list.
#' 
#' @param listPred List of cell type proportions predictions.
#' @param columnLabels Vector of cell type labels 
#' (e.g. "type1", "type2", etc.).
#' @param rowLabels Vector of sample id labels 
#' (e.g. "sample1", "sample2", etc.).
#' 
#' @examples
#' exampleData <- getDeconvolutionExampleData()
#' 
#' @returns Example data as list.
#' 
#' @export
parseDeconvolutionPredictionsResults <- function(listPred, columnLabels, 
                                                  rowLabels){
  if(is(columnLabels, "NULL")){
    columnLabels <- seq(length(listPred[[1]]))}
  if(is(rowLabels, "NULL")){
    rowLabels <- seq(length(listPred))}
  tablePred <- do.call(rbind, listPred)
  tablePred <- apply(tablePred, 1, function(ri){ri/sum(ri)}) 
  tablePred <- t(tablePred)
  colnames(tablePred) <- columnLabels
  rownames(tablePred) <- rowLabels
  ## convert
  tablePred <- cellProportionsPredictions(tablePred, 
                                          columnLabels,
                                          rowLabels)
  return(tablePred)
}
