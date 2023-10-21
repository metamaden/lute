#!/usr/bin/env R

### Author: Sean Maden
###
### Functions to get randomized or synthetic datasets.
###
###
###

#' proportionsVectorsList
#'
#' Get complementary proportions for k types. The first type k1 is the vector of
#' proportions for the first type. The remaining types up to totalCellTypesK are based on
#' the reverse of k1. Types k > 1 are assumed to have equal proportions 
#' complementary to k1.
#' 
#' For k1=c(0, 0.5, 1), totalCellTypesK=2 will generate an additional type with 
#' proportions c(1, 0.5, 0).
#' 
#' For the same k1 above, totalCellTypesK=3, will generate 2 types with the same 
#' proportions as c(0.5, 0.25, 0).
#'
#' @param totalCellTypesK Total number of cell types to simulate.
#' @param firstCellTypeProportions Vector of first cell type proportions. If NULL, uses 
#' seq(1e-3, 1-1e-3, 1e-3).
#'
#' @returns lpv, a list of proportions vectors for simulation iterations.
#'
#' @examples
#' proportionsVectorsList(firstCellTypeProportions=c(0, 0.5, 1))
#'
#' @export
proportionsVectorsList <- function(
    totalCellTypesK=2, firstCellTypeProportions=NULL){
  if(is(firstCellTypeProportions, "NULL")){
    firstCellTypeProportions <- seq(0, 1, 1e-3)}
  numberOfIterations <- length(firstCellTypeProportions)
  secondCellTypeProportions <- rev(firstCellTypeProportions)/(totalCellTypesK-1)
  proportionsVectorsList <- lapply(seq(numberOfIterations), function(ii){
    c(firstCellTypeProportions[ii], 
      rep(secondCellTypeProportions[ii], totalCellTypesK-1))})
  return(proportionsVectorsList)
}

#' randomMarkersVectorsList
#'
#' Get randomized markers using Poisson distribution sampling. For a given K,
#' we assume "positive" markers have higher values than for non-K types, and 
#' thus we sample from 2 different Poisson distributions defined by different
#' lambda values (e.g. arguments lambdaMean, lambdaMeanNegative). WE also use argument 
#' markerIndexVector to define total markers as length(markerIndexVector) and the marker balance as
#' relative counts of each type index.
#' 
#' For example, if gindex is c(1, 1, 2), we define 3 total markers, 2 positive
#' markers for type 1 (negative for type 2) and a single positive marker for 
#' type 2 (negative for type 1).
#' 
#' @param markerIndexVector Vector of marker indices. Index values correspond to the k types,
#' and each index position represents a marker (e.g. c(1,2,2) means two markers 
#' for the second type, etc.).
#' @param numberIterations Total simulation iterations.
#' @param lambdaMean Value of lambda (Poisson dist. mean) for "positive" marker 
#' status (e.g. mean of dist. for k when marker is positive for k, negative for 
#' not-k). This is passed to the argument mu when method is "nbinom".
#' @param lambdaMeanNegative Value of lambda (Poisson dist. mean) for "negative" marker 
#' status (e.g. mean of dist. for k when marker is positive for not-k, negative 
#' for k). This is passed to the argument mu when method is "nbinom".
#' @param method Type of randomization method to use. Accepts either "poisson"
#' for poisson distribution (see `?rpois` for details), or "nbinom" for the 
#' negative binomial (a.k.a. gamm poisson) distribution (see `?rnbinom` for 
#' details).
#' @param gammaSize The gamma distribution magnitude for "positive" markers. 
#' This is applied when the "nbinom" method is used.
#' @param gammaSizeNegative The gamma distribution magnitude for "negative" markers. 
#' This is applied when the "nbinom" method is used.
#' @returns Listed lgv object containing the randomized marker values across 
#' types.
#' 
#' @importFrom stats rnbinom
#' 
#' @examples 
#' randomMarkersVectorsList(markerIndexVector=c(rep(1, 10), rep(2, 5)))
#'
#' @export
randomMarkersVectorsList <- function(markerIndexVector, numberIterations=1, 
                                     lambdaMean=25, lambdaMeanNegative=2,
                                     method="nbinom", gammaSize=10,
                                     gammaSizeNegative=10){
  totalCellTypesK <- length(unique(markerIndexVector))
  markersVectorList <- lapply(seq(totalCellTypesK), function(ki){
    markersVectorG <- markerIndexVector
    whichPositive <- which(markerIndexVector == ki)
    whichNegative <- which(!markerIndexVector == ki)
    if(method == "poisson"){
      markersVectorG[whichPositive] <- 
        rpois(lambda=lambdaMean, n=length(whichPositive))
      markersVectorG[whichNegative] <- 
        rpois(lambda=lambdaMeanNegative, n=length(whichNegative))      
    } else if(method == "nbinom"){
      markersVectorG[whichPositive] <- 
        rnbinom(size=gammaSize, mu=lambdaMean, n=length(whichPositive))
      markersVectorG[whichNegative] <- rnbinom(
        size=gammaSizeNegative, mu=lambdaMeanNegative, n=length(whichNegative))
    } else{
      stop("Error, invalid method.")
    }
    markersVectorG
  })
  returnList <- lapply(seq(numberIterations), function(ii){markersVectorList})
  return(returnList)
}

#' randomSingleCellExperiment
#'
#' Make a random object of type SingleCellExperiment. Uses the negative binomial 
#' distribution to randomly generate gene expression data for simulated cells.
#'
#' @param numberGenes Number of genes to randomize.
#' @param numberCells Numnber of cells to randomize.
#' @param numberTypes Number of cell types to annotate.
#' @param fractionTypes Vector of fractions by type.
#' @param dispersion Disperison of gene expression. If NULL, uses the mean from 
#' expressionMean
#' @param expressionMean Poisson dist mean for random expression data.
#' @param seedNumber Seed value for randomization of expression data.
#' @param naInclude Whether to include random NA values.
#' @param naFraction Fraction of NA values to include.
#' @param zeroInclude Whether to include random zero-count values.
#' @param zeroFraction Fraction of zero-count values to include.
#' @param verbose Whether to show verbose status messages.
#'
#' @return New randomized SingleCellExperiment object.
#' 
#' @importFrom stats rnbinom
#' @importFrom S4Vectors metadata
#' 
#' @examples 
#' singleCellExperiment <- randomSingleCellExperiment()
#' 
#' @export
randomSingleCellExperiment <- function(
    numberGenes=20, numberCells=12, numberTypes=2, fractionTypes=NULL, 
    dispersion=NULL, expressionMean=10, naInclude=FALSE, naFraction=0.2, 
    zeroInclude=FALSE, zeroFraction=0.2, verbose=FALSE, seedNumber=0
    ){
  if(verbose){message("Getting random expression data...")}
  if(is(dispersion, "NULL")){dispersion <- expressionMean}
  matrixData <- rnbinom(n=(numberCells*numberGenes), 
                  size=dispersion, mu=expressionMean)
  if(naInclude){ # manually add NAs
    if(verbose){message("Including NA values...")}
    numberNAs <- round(length(matrixData)*naFraction, digits=0)
    indexNAs <- sample(seq(length(matrixData)), numberNAs)
    matrixData[indexNAs] <- NA
  }
  if(zeroInclude){ # manually add zero counts
    if(verbose){message("Including NA values...")}
    numberZeroes <- round(length(matrixData)*zeroFraction, digits=0)
    zeroIndex <- sample(seq(length(matrixData)), numberZeroes)
    matrixData[zeroIndex] <- 0
  }
  expressionCounts <- matrix(matrixData, ncol=numberCells, nrow=numberGenes)
  if(verbose){message("Getting new colData...")}
  cellVector <- paste0("cell.barcode.", seq(numberCells))
  cellsPerType <- round(numberCells/numberTypes, 0)
  
  if(is(fractionTypes, "NULL")){
    fractionTypes <- rep((1/numberTypes), numberTypes)}
  typeVector <- paste0("type", seq(numberTypes))
  typeVector <- unlist(lapply(seq(length(typeVector)), function(ti){
    num <- fractionTypes[ti]*numberCells; rep(typeVector[ti], num)
  }))
  
  newColData <- data.frame(cell.id=cellVector, celltype=typeVector)
  colnames(expressionCounts) <- cellVector
  if(verbose){message("Getting new rowData...")}
  genev <- paste0("gene", seq(nrow(expressionCounts)))
  newRowData <- data.frame(gene.id=genev)
  rownames(expressionCounts) <- genev
  ## manage new metadata
  descriptionString <- "random SingleCellExperiment made using randomSingleCellExperiment()"
  metadataList <- list(description=descriptionString)
  if(verbose){message("Making new sce object...")}
  singleCellExperiment <- SingleCellExperiment::SingleCellExperiment(
    assays=list(counts=expressionCounts), colData=newColData, 
    rowData=newRowData, metadata=metadataList
  )
  return(singleCellExperiment)
}
