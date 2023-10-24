#!/usr/bin/env R

### Author: Sean Maden

#' findmarkersParam-class
#'
#' class definition for findmarkersParam, which uses scran::findMarkers()
#' 
#' @include lute_generics.R
#' @include typemarkersParam-class.R
#' 
#' @details Main constructor for class \linkS4class{findmarkersParam}.
#' @rdname findmarkersParam-class
#' @seealso \linkS4class{typemarkersParam}
#' 
#' @param assayName Name of expression matrix in SingleCellExperiment assays 
#' (e.g. "counts").
#' @param singleCellExperiment Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @param cellTypeVariable Name of cell type variable in SingleCellExperiment 
#' coldata.
#' @param testType Test type (see \code{?findMarkers} for options).
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- random_singleCellExperiment()
#' newParam <- findmarkersParam(singleCellExperiment=singleCellExperimentExample, 
#' cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' 
#' @aliases 
#' FindmarkersParam-class, findMarkersParam-class
#' 
#' @returns New object.
#'
setClass("findmarkersParam", contains="typemarkersParam", 
         slots=c(assayName="character", 
                 singleCellExperiment="SingleCellExperiment",
                 cellTypeVariable="character",
                 testType="character"))

#' Make new object of class findmarkersParam
#'
#' Main constructor for class \linkS4class{findmarkersParam}.
#'
#' @param assayName Name of expression matrix in SingleCellExperiment assays 
#' (e.g. "counts").
#' @param singleCellExperiment Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @param cellTypeVariable Name of cell type variable in SingleCellExperiment 
#' coldata.
#' @param markersPerType Number of top markers to get per cell type.
#' @param testType Test type (see \code{?findMarkers} for options).
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @returns Object of class \linkS4class{findmarkersParam}
#' 
#' @seealso \linkS4class{typemarkersParam}
#'
#' @details Main class for mapping arguments to the findMarkers method 
#' implemented as \code{scran::findMarkers()}.
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- random_singleCellExperiment()
#' newParam <- findmarkersParam(singleCellExperiment=singleCellExperimentExample, 
#' cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' 
#' @export
findmarkersParam <- function(singleCellExperiment, 
                             assayName="counts", 
                             cellTypeVariable="cellType",
                             testType="wilcox",
                             markersPerType=20, 
                             returnInfo=FALSE) {
  new("findmarkersParam", singleCellExperiment=singleCellExperiment, 
      assayName=assayName, cellTypeVariable=cellTypeVariable, testType=testType, 
      markersPerType=markersPerType, returnInfo=returnInfo)
}

#' Cell type markers method for findmarkersParam
#'
#' Defines the typemarkers method for \linkS4class{findmarkersParam}.
#'
#' @param object An object of class \linkS4class{findmarkersParam} (see 
#' \code{?findmarkersParam}).
#'
#' @details Takes an object of class \linkS4class{findmarkersParam} as input, 
#' returning either a vector of cell type gene markers, or (if 
#' \code{returnInfo == TRUE}) a list containing such a vector along with 
#' original function outputs.
#' 
#' @importFrom scran findMarkers
#' @importFrom dplyr %>%
#' 
#' @examples 
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- random_singleCellExperiment()
#' newParam <- findmarkersParam(singleCellExperiment=singleCellExperimentExample, 
#' cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#'
#' @returns Returns the top available markers, with type-specific marker filters,
#' as either a vector of marker IDs or a results list.
#'
#' @export
setMethod("typemarkers", signature(object="findmarkersParam"), function(object){
  singleCellExperiment <- object[["singleCellExperiment"]]
  cellTypeVariable <- object[["cellTypeVariable"]]
  assayName <- object[["assayName"]]
  markersPerType <- object[["markersPerType"]]
  testType <- object[["testType"]]
  ## get marker results
  resultList <- markerTableList <- list(); markersFilter <- c()
  uniqueCellTypes <- singleCellExperiment[[cellTypeVariable]] %>% 
    as.character() %>% unique()
  for(type in uniqueCellTypes){
    singleCellExperiment <- 
      singleCellExperiment[!rownames(singleCellExperiment) %in% markersFilter,]
    message("selecting among ",nrow(singleCellExperiment),
            " genes for markers of type: ", type, "...")
    resultList[[type]] <- findMarkers(x=singleCellExperiment, 
                                      group=singleCellExperiment[[cellTypeVariable]],
                                      assay.type=assayName,
                                      test.type=testType)[[type]]
    dfIterate <- resultList[[type]][,seq(4)]
    summaryColname <- 
      colnames(dfIterate)[grepl("summary\\..*", colnames(dfIterate))]
    newColname <- paste0("abs.", summaryColname)
    dfIterate$abs.summary <- 
      dfIterate[,summaryColname] %>% as.numeric() %>% abs()
    dfIterate$cellType.target <- type
    dfIterate$gene <- rownames(dfIterate)
    dfIterate <- as.data.frame(dfIterate)
    ## get filtered markers
    markerTableList[[type]] <- dfIterate %>% dplyr::arrange(abs.summary) %>% 
      dplyr::top_n(n=markersPerType) %>% as.data.frame()
    markersFilter <- markerTableList[[type]]$gene
  }
  topMarkerTable <- do.call(rbind, markerTableList) %>% as.data.frame()
  topMarkersVector <- topMarkerTable$gene
  ## parse returnInfo
  returnList <- topMarkersVector %>% unique()
  if(object[["returnInfo"]]){
    returnList <- list(
      markers=topMarkersVector, result.info=markerTable, metadata=object)}
  return(returnList)
})

#' Show generic behavior for object of class \linkS4class{findmarkersParam}
#' @param object An object of class \linkS4class{findmarkersParam} (see 
#' \code{?findmarkersParam}).
#' @details Method for behavior of show generic when called for object of class 
#' \linkS4class{findmarkersParam}
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#' singleCellExperimentExample <- random_singleCellExperiment()
#' newParam <- findmarkersParam(singleCellExperiment=singleCellExperimentExample, 
#' cellTypeVariable="celltype", markersPerType=5)
#' markers <- typemarkers(newParam)
#' 
#' @returns Shows object summaries.
#' 
#' @export
setMethod("show", signature(object="findmarkersParam"), function(object){
  show(object)
})
