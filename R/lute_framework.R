#!/usr/bin/env r

### Author: Sean Maden
###
### Framework functions and utilities for the lute package.

#' lute framework
#'
#' Obtain cell type markers and proportion predictions from various algorithms. 
#' Allows flexible data types and standard application of cell size scale 
#' factors.
#' 
#' @param singleCellExperiment Object of type SingleCellExperiment. Optional (see argument z).
#' @param referenceExpression Signature matrix of cell type-specific signals. Optional (see 
#' argument singleCellExperiment).
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples. Optional (see argument y.se).
#' @param bulkSummarizedExperiment SummarizedExperiment or similar data type containing the bulk signals
#' matrix in its assays (e.g. accessible with assays(y.se)[[assayName]] using
#' the provided assayName argument). Optional (see argument y).
#' @param cellScaleFactors Cell size factor transformations of length equal to the K cell types 
#' to deconvolve. Optional, if not provided, uses equal weights for types.
#' @param assayName Name of expression matrix in singleCellExperiment, and optionally y.se, 
#' assays. Optional (e.g. "counts"; see arguments singleCellExperiment, y.se).
#' @param cellTypeVariable Name of cell type variable in singleCellExperiment coldata.
#' @param markersPerType Number of top markers to get per cell type.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' @param typemarkerAlgorithm Which type-specific marker selection algorithm to 
#' use. If NULL, skips type marker analyses.
#' @param deconvolutionAlgorithm Where deconvolution algorithm to use. If NULL, 
#' skips deconvolution.
#' @param verbose Whether to show verbose status messages.
#' 
#' @returns A list containing results returned from type marker selection and
#' deconvolution runs, with additional information returned if
#'  \code{returnInfo == TRUE}.
#' 
#' @details Main function to use the lute deconvolution framework. Manages data
#' conversions and mappings to deconvolution experiment steps, including setup,
#' gene marker identification, and main deconvolution runs.
#' 
#' Support is provided for \linkS4class{SummarizedExperiment}-type or 
#' matrix-type inputs for the Z signature matrix (see referenceExpression 
#' argument) and Y bulk signals matrix (see bulkExpression arguments). Note,
#' both Z and Y need to be provided or derivable in order to run deconvolution.
#' 
#' @examples 
#' # get example bulk data
#' bulkExpression <- getDeconvolutionExampleData()$reference
#' 
#' # get example singleCellExperiment
#' singleCellExperiment <- randomSingleCellExperiment()[seq(10),]
#' 
#' # get framework results
#' experiment.results <- lute(
#' singleCellExperiment=singleCellExperiment, 
#' referenceExpression=referenceExpression, typemarkerAlgorithm=NULL
#' )
#' 
#' @export
lute <- function(singleCellExperiment=NULL, 
                 referenceExpression=NULL, 
                 bulkExpression=NULL, 
                 bulkSummarizedExperiment=NULL, 
                 cellScaleFactors=NULL, 
                 returnInfo=FALSE, 
                 markersPerType=20,
                 assayName="counts", 
                 cellTypeVariable="celltype",
                 typemarkerAlgorithm="findmarkers", 
                 deconvolutionAlgorithm="nnls",
                 verbose=TRUE){
  resultsList <- list()
  if(!is(typemarkerAlgorithm, "NULL")){
    if(verbose){message("Parsing marker gene arguments...")}
    typemarkerResults <- markerVector <- map_typemarker_algorithm(
      algorithm=typemarkerAlgorithm,
      singleCellExperiment=singleCellExperiment,
      assayName=assayName,
      cellTypeVariable=cellTypeVariable,
      markersPerType=markersPerType,
      returnInfo=returnInfo)
    if(is(typemarkerResults, "list")){
      markerVector <- typemarkerResults[["markers"]]}
    if(verbose){message("Filtering singleCellExperiment...")}
    filterSingleCellExperiment <- 
      rownames(singleCellExperiment) %in% markerVector
    singleCellExperiment <- singleCellExperiment[filterSingleCellExperiment,]
    referenceExpression <- referenceFromSingleCellExperiment(
      singleCellExperiment, assayName, cellTypeVariable)
    resultsList[["typemarkerResults"]] <- typemarkerResults
  }
  bulkCondition <- is(bulkExpression, "NULL") & !is(bulkSummarizedExperiment, "NULL")
  if(bulkCondition){bulkExpression <- assays(bulkSummarizedExperiment)[[assayName]]}
  if(!is(deconvolutionAlgorithm, "NULL")){
    if(is(referenceExpression, "NULL")){
      referenceExpression <- referenceFromSingleCellExperiment(
        singleCellExperiment, assayName, cellTypeVariable)
    }
    if(is(cellScaleFactors, "NULL")){
      cellScaleFactors <- rep(1, ncol(referenceExpression))
      names(cellScaleFactors) <- colnames(referenceExpression)
    }
    if(verbose){message("Parsing deconvolution arguments...")}
    deconvolutionResults <- map_deconvolution_algorithm(
      algorithm=deconvolutionAlgorithm,
      referenceExpression=referenceExpression, 
      bulkExpression=bulkExpression, cellScaleFactors=cellScaleFactors, 
      returnInfo=returnInfo)
    resultsList[["deconvolutionResults"]] <- deconvolutionResults
  }
  return(resultsList)
}

#'
#'
map_typemarker_algorithm <- function(algorithm, singleCellExperiment, assayName, 
                                     cellTypeVariable, markersPerType, 
                                     returnInfo){
  if(algorithm %in% c("mr", "MR", "Mr", "deconvobuddies", "DeconvoBuddies", 
                         "meanratio", "meanRatio", "Meanratio", "MeanRatio",
                         "meanratios", "meanRatios", "Meanratio", "MeanRatio")){
    message("Using meanratiosParam...")
    typeMarkerString <- "meanratiosParam"
  } else if(algorithm %in% c("findmarkers", "findmarker", "Findmarker", 
                             "Findmarkers", "findMarkers", "FindMarkers")){
    message("Using meanratiosParam...")
    typeMarkerString <- "findmarkersParam"
  } else{
    message("Warning, unidentified marker selection algorithm provided. ",
            "Skipping marker selection")
    returnString <- "FALSE"
  }
  typeMarkerString <- paste0(typeMarkerString, 
      "(singleCellExperiment=singleCellExperiment, ",
      "assayName='",assayName,"', ",
      "markersPerType=", markersPerType, ", ",
      "cellTypeVariable='",cellTypeVariable,"',",
      "returnInfo=", returnInfo, ")")
  newParam <- eval(parse(text=typeMarkerString))
  return(typemarkers(newParam))
}

#'
#'
map_deconvolution_algorithm <- function(algorithm, bulkExpression, 
                                        referenceExpression, cellScaleFactors, 
                                        returnInfo){
  if(algorithm %in% 
     c("nnls", "Nnls", "NNLS", "nnlsParam", "NNLSParam", "NnlsParam")){
    message("Using NNLS...")
    deconvolutionString <- "nnlsParam"
  } else if(algorithm %in% 
        c("music", "Music", "MuSiC", "musicParam", "MusicParam", "MuSiCParam")){
    message("Using MuSiC...")
    deconvolutionString <- "musicParam"
  } else if(algorithm %in% 
  c("music2", "Music2", "MuSiC2", "musicParam2", "MusicParam2", "MuSiCParam2")){
    message("Using MuSiC2...")
    deconvolutionString <- "music2Param"
  } else if(algorithm %in% c("epic", "Epic", "EPIC", "epicParam", "EpicParam", "EPICParam")){
    message("Using EPIC...")
    deconvolutionString <- "epicParam"
  } else if(algorithm %in% c("bisque", "Bisque", "BISQUE", "bisqueParam", "BisqueParam", "BISQUEParam")){
    message("Using Bisque...")
    deconvolutionString <- "bisqueParam"
  } else if(algorithm %in% c("deconrnaseq", "Deconrnaseq", "DECONRNASEQ", 
                                "deconrnaseqParam", "DeconrnaseqParam", "DECONRNASEQParam")){
    message("Using DeconRNASeq...")
    deconvolutionString <- "deconrnaseqParam"
  } else if(algorithm %in% c("scdc", "Scdc", "SCDC", 
                                "scdcParam", "ScdcParam", "SCDCParam")){
    message("Using SCDC...")
    deconvolutionString <- "scdcParam"
  } else{
    stop("Error, unidentified deconvolution algorithm provided. ")
  }
  deconvolutionString <- paste0(deconvolutionString, 
                                "(bulkExpression=bulkExpression, ",
                                "referenceExpression=referenceExpression, ",
                                "cellScaleFactors=cellScaleFactors, ",
                                "returnInfo=", returnInfo, ")")
  newParameter <- eval(parse(text=deconvolutionString))
  return(deconvolution(newParameter))
}
