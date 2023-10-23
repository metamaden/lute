#!/usr/bin/env R

### Author: Sean Maden
###
### Utilities and misingleCellExperimentllaneous functions supporting the lute package for deconvolution experiments.
###

#' get_celltypes_from_sce
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
#' exampleList <- getDeconvolutionExampleData()
#'
#' @export
get_celltypes_from_sce <- function(
    singleCellExperiment, cellTypeVariable="celltype"){
  celltype.vector <- as.data.frame(
    SummarizedExperiment::colData(singleCellExperiment))[,cellTypeVariable]
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

#' ypb_from_sce
#'
#' Get pseudobulk from a SingleCellExperiment object.
#'
#' @param singleCellExperiment An object of type 
#' \linkS4class{SingleCellExperiment}.
#' @param assayName Name of expression matrix in \code{singleCellExperiment} 
#' assays.
#' @param cellTypeVariable Variable name for cell type labels in 
#' \code{singleCellExperiment} 
#' coldata.
#' @param sampleIdVariable Variable name for sample/group ID labels in 
#' \code{singleCellExperiment} coldata.
#' @param cellScaleFactors Vector of cell type size scale factors. Optional.
#' @returns Matrix of simulated bulk convoluted signals.
#' 
#' @examples
#' singleCellExperimentExample <- randomSingleCellExperiment()
#' ypb_from_sce(singleCellExperimentExample)
#' 
#' @export
ypb_from_sce <- function(singleCellExperiment, assayName="counts", 
                         cellTypeVariable="celltype", 
                         sampleIdVariable=NULL, cellScaleFactors=NULL){
  groupNumber <- 1
  uniqueGroupIdVector <- ""
  if(!is(sampleIdVariable, "NULL")){
    group.id.vector <- singleCellExperiment[[sampleIdVariable]]
    uniqueGroupIdVector <- group.id.vector 
    uniqueGroupIdVector <- unique(uniqueGroupIdVector)
    uniqueGroupIdVector <- as.character(uniqueGroupIdVector)
    groupNumber <- length(uniqueGroupIdVector)
  }
  cellTypesList <- get_celltypes_from_sce(
    singleCellExperiment=singleCellExperiment, 
    cellTypeVariable=cellTypeVariable)
  numberTypes <- length(cellTypesList[["unique.types"]])
  pseudobulkList <- lapply(uniqueGroupIdVector, function(group.id){
    filterSingleCellExperiment <- singleCellExperiment
    if(groupNumber > 1){
      filter.group <- singleCellExperiment[[sampleIdVariable]]==group.id
      filterSingleCellExperiment <- singleCellExperiment[,filter.group]
    }
    if(is(cellScaleFactors, "NULL")){
      cellScaleFactors <- rep(1, numberTypes)
      names(cellScaleFactors) <- cellTypesList[["unique.types"]]
    }
    
    referenceExpressionZnew <- 
      referenceFromSingleCellExperiment(
        filterSingleCellExperiment, assayName, cellTypeVariable)
    cellTypeProportions <- table(cellTypesList[["character"]])
    cellTypeProportions <- prop.table(cellTypeProportions)
    orderProportions <- match(
      names(cellTypeProportions), cellTypesList[["unique.types"]])
    orderProportions <- order(orderProportions)
    cellTypeProportions <- cellTypeProportions[orderProportions]
    referenceExpressionZSnew <- .zstransform(
      referenceExpressionZnew, cellScaleFactors)
    bulkExpressionPseudobulk <- t(t(cellTypeProportions) %*% 
                                    t(referenceExpressionZSnew))
    return(bulkExpressionPseudobulk)
  })
  pseudobulkTable <- do.call(cbind, pseudobulkList)
  pseudobulkTable <- as.data.frame(pseudobulkTable)
  if(groupNumber > 1){
    colnames(pseudobulkTable) <- uniqueGroupIdVector
  } else{
    colnames(pseudobulkTable) <- "singleCellExperiment.pseudobulk"
  }
  return(pseudobulkTable)
}

#' z_matrix_from_sce
#' 
#' Calculate a Z signature matrix (referenceExpression) from object of type 
#' \linkS4class{SingleCellExperiment}.
#' 
#' @param singleCellExperiment An object of type \linkS4class{SingleCellExperiment}.
#' @param assayName Name of expression matrix in \code{singleCellExperiment} assays (e.g. 
#' "counts").
#' @param cellTypeVariable Variable name for cell type labels in \code{singleCellExperiment} 
#' coldata (e.g. "type1", "type2", etc.). 
#' @param summaryMethod Summary statistic function to use.
#' @details Calculate a Z signature matrix from object of type 
#' \linkS4class{SingleCellExperiment}.
#' @returns New Z signature matrix.
#' 
#' @examples
#' singleCellExperiment.example <- randomSingleCellExperiment()
#' signature_matrix_from_singleCellExperiment(singleCellExperiment.example)
#' 
#' @export
z_matrix_from_sce <- function(singleCellExperiment,
                              cellTypeVariable="celltype",
                              summaryMethod="mean",
                              assayName="counts"){
  ## gets the z signature matrix from an singleCellExperiment object
  expressionMatrix <- assays(singleCellExperiment)[[assayName]]
  expressionMatrix <- as.matrix(expressionMatrix)
  newColData <- colData(singleCellExperiment)
  uniqueCellTypes <- unique(newColData[,cellTypeVariable])
  uniqueCellTypes <- uniqueCellTypes[order(uniqueCellTypes)]
  referenceExpression <- do.call(cbind, 
                                 lapply(uniqueCellTypes, 
                                        function(cellTypeIndex){
    filterIndex <- newColData[,cellTypeVariable]==cellTypeIndex
    if(summaryMethod == "mean"){
      DelayedArray::rowMeans(expressionMatrix[,filterIndex])
    } else{
      Biobase::rowMedians(expressionMatrix[,filterIndex])
    }
  }))
  colnames(referenceExpression) <- uniqueCellTypes
  return(referenceExpression)
}

#' referenceFromSingleCellExperiment
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
#' exampleList <- getDeconvolutionExampleData()
#' 
#' @export
referenceFromSingleCellExperiment <- function(
    singleCellExperiment, assayName="counts", cellTypeVariable="celltype"){
  typesList<- get_celltypes_from_sce(
    singleCellExperiment=singleCellExperiment, 
    cellTypeVariable=cellTypeVariable)
  expressionMatrix <- as.matrix(
    assays(
      singleCellExperiment)[[assayName]])
  referenceExpressionZnew <- 
    do.call(cbind, lapply(typesList[["unique.types"]], function(typei){
    datav <- expressionMatrix[,typesList[["character"]]==typei]
    .z_operator(datav)
  }))
  colnames(referenceExpressionZnew) <- typesList[["unique.types"]]
  rownames(referenceExpressionZnew) <- rownames(singleCellExperiment)
  return(referenceExpressionZnew)
}

.zstransform <- function(referenceExpression, cellScaleFactors){
  markersReferenceExpression <- colnames(referenceExpression)
  markersFactors <- names(cellScaleFactors)
  reorderFactors <- order(match(markersFactors, markersReferenceExpression))
  cellScaleFactors <- cellScaleFactors[reorderFactors]
  if(identical(colnames(referenceExpression),names(cellScaleFactors))){
    referenceTransformed <- sweep(
      referenceExpression, 2, cellScaleFactors, FUN="*")
  } else{
    allTypes <- unique(c(markersReferenceExpression, markersFactors))
    sharedTypes <- intersect(markersReferenceExpression, markersFactors)
    outersectTypes <- c(markersReferenceExpression[
      !markersReferenceExpression %in% markersFactors],
      markersFactors[!markersFactors %in% markersReferenceExpression])
    message(
      "Cell types not shared in referenceExpression and cellScaleFactors:", 
      outersectTypes)
    stop("Error matching types in referenceExpression and cellScaleFactors.")
  }
  return(referenceTransformed)
}

.z_operator <- function(expressionDataVector){
  rowMeans(expressionDataVector)
}

#' getDeconvolutionExampleData
#' 
#' Make example data for deconvolution.
#' 
#' @param cellScaleFactors Vector of cell scale factors
#' @param numberBulkSamples Number of bulk samples.
#' @param numberMarkers Number of cell type markers.
#' @param numberTypes Number of cell types.
#' @returns Example data as list.
#' 
#' @importFrom stats rpois
#' 
#' @examples
#' exampleData <- getDeconvolutionExampleData()
#' 
#' @export
getDeconvolutionExampleData <- function(
    cellScaleFactors=c(1,10),numberBulkSamples=2, numberMarkers=10, numberTypes=2
  ){
  if(!length(cellScaleFactors)==numberTypes){
    stop("Error, cellScaleFactors length should equal numberTypes.")}
  bulkExpression <- matrix(
    rpois(n=numberMarkers*numberBulkSamples, lambda=seq(0, 50, 5)), 
    ncol=numberBulkSamples)
  referenceExpression <- matrix(
    rpois(n=numberTypes*numberMarkers, lambda=seq(0, 50, 5)), 
    ncol=numberTypes)
  rownames(bulkExpression) <- rownames(referenceExpression) <- 
    paste0("marker", seq(numberMarkers))
  colnames(referenceExpression) <- paste0("type", seq(numberTypes))
  colnames(bulkExpression) <- paste0("sample", seq(numberBulkSamples))
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

#' getDeconvolutionExampleDataBisque
#'
#' Get example data for Bisque algorithm.
#'
#' @param numberBulkSamples Number of bulk samples.
#' @param numberMarkers Number of cell type markers.
#' @param numberCells Number of cells.
#' @param numberTypes Number of cell types.
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
                                               numberTypes=2){
  exampleList <- getDeconvolutionExampleData(
    numberBulkSamples=numberBulkSamples,
    numberMarkers=numberMarkers,
    numberTypes=numberTypes)
  bulkExpression <- exampleList[["bulkExpression"]]
  colnames(bulkExpression) <- c(paste0("sample", seq(numberBulkSamples/2)), 
                   paste0("bulk", seq(numberBulkSamples/2)))
  dfBulkPheno <- data.frame(SubjectName=colnames(bulkExpression))
  rownames(dfBulkPheno) <- colnames(bulkExpression)
  bulkExpressionSet <- ExpressionSet(assayData=bulkExpression, 
                          phenoData=AnnotatedDataFrame(dfBulkPheno))
  singleCellExperiment <- randomSingleCellExperiment(numberGenes=numberMarkers, 
                    numberCells=numberCells, 
                    numberTypes=numberTypes)
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

#' getDeconvolutionExampleDataSnewColDataC
#' 
#' Get example data for SnewColDataC
#'
#' @returns Example data as list.
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom BiocGenerics counts
#' 
#' @examples
#' exampleData <- getDeconvolutionExampleDataSnewColDataC()
#' 
#' @export
getDeconvolutionExampleDataSnewColDataC <- function(){
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
  singleCellExperiment <- randomSingleCellExperiment(numberGenes=10, 
                                                     numberCells=300, 
                                                     numberTypes=4)
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
  singleCellExperiment <- randomSingleCellExperiment(numberGenes=10, 
                                                     numberCells=300, 
                                                     numberTypes=2)
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


#' getDeconvolutionExampleDataSnewColDataC
#' 
#' Get example data for SnewColDataC
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
