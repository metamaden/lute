#!/usr/bin/env R

###
### Functions to convert between data object classes.
###

#' sce_to_eset
#' Convert SingleCellExperiment to ExpressionSet.
#' @param singleCellExperiment Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @param assayName Name of assay to store in new eset.
#' @returns ExpressionSet.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase pData
#' @importFrom Biobase AnnotatedDataFrame
#'
#' @examples
#' sce <- randomSingleCellExperiment()
#' sce_to_eset(sce, "counts")
#' 
#' @export
sce_to_eset <- function(singleCellExperiment, assayName="counts"){
	ExpressionSet <- ExpressionSet(
	  assayData=assays(singleCellExperiment)[[assayName]],
	  phenoData=AnnotatedDataFrame(
	    as.data.frame(
	      colData(singleCellExperiment))))
	return(ExpressionSet)
}

#' eset_to_sce
#' Convert ExpressionSet to SingleCellExperiment.
#' @param expressionSet Object of type ExpressionSet (see \code{?ExpressionSet}).
#' @param assayName Name of new assay in new SingleCellExperiment object.
#' @returns ExpressionSet.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase pData
#' @importFrom Biobase exprs
#' @importFrom SingleCellExperiment SingleCellExperiment
#' 
#' @examples
#' expressionSet <- getDeconvolutionExampleDataBisque()$singleCellExpressionSet
#' eset_to_sce(expressionSet)
#' 
#' @export
eset_to_sce <- function(expressionSet, assayName="counts"){
  assaysList <- list(assayName=exprs(expressionSet))
  names(assaysList) <- assayName
	 singleCellExperimentNew <- SingleCellExperiment(assays=assaysList,
	                                colData=DataFrame(
	                                  as.matrix(pData(expressionSet))))
	return(singleCellExperimentNew)
}

#' sce_to_se
#' Convert SingleCellExperiment to SummarizedExperiment.
#'
#' @param singleCellExperiment Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @returns SummarizedExperiment.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors metadata
#' 
#' @examples
#' sce <- randomSingleCellExperiment()
#' sce_to_se(sce)
#' @export
sce_to_se <- function(singleCellExperiment){
  summarizedExperiment <- SummarizedExperiment(
    assays=assays(singleCellExperiment),
    colData=colData(singleCellExperiment),
    metadata=metadata(singleCellExperiment))
	return(summarizedExperiment)
}

#' se_to_sce
#'
#' Convert SummarizedExperiment to SingleCellExperiment.
#'
#' @param summarizedExperiment Object of type SummarizedExperiment (see 
#' \code{?SummarizedExperiment}).
#' @returns New SingleCellExperiment object.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment SingleCellExperiment
#' 
#' @examples
#' se_to_sce(SummarizedExperiment())
#' 
#' @export
se_to_sce <- function(summarizedExperiment){
	singleCellExperiment <- SingleCellExperiment(
	  assays=assays(summarizedExperiment),
	  colData=colData(summarizedExperiment),
	  metadata=metadata(summarizedExperiment))
	return(singleCellExperiment)
}

#' eset_to_se
#'
#' Convert ExpressionSet to SummarizedExperiment.
#'
#' @param expressionSet Object of type ExpressionSet (see \code{?ExpressionSet}).
#' @param assayName Name of assay to store in new SummarizedExperiment object.
#' @returns New object of type SummarizedExperiment.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase pData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' 
#' @examples
#' expressionSet <- getDeconvolutionExampleDataBisque()$singleCellExpressionSet
#' eset_to_se(expressionSet, "counts")
#' 
#' @export
eset_to_se <- function(expressionSet, assayName="counts"){
  assaysList <- list(assayName=exprs(expressionSet))
  names(assaysList) <- assayName
	summarizedExperimentNew <- SummarizedExperiment(
	  assays=assaysList, colData=DataFrame(as.matrix(pData(expressionSet))))
	return(summarizedExperimentNew)
}

#' se_to_eset
#' 
#' Convert SummarizedExperiment to ExpressionSet.
#' 
#' @param summarizedExperiment Object of type SummarizedExperiment (see 
#' \code{?SummarizedExperiment}).
#' @param assayName Name of assay to store in new ExpressionSet object.
#' @returns New object of type ExpressionSet.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase pData
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#'
#' @examples
#' summarizedExperiment <- sce_to_se(randomSingleCellExperiment())
#' se_to_eset(summarizedExperiment)
#' 
#' @export
se_to_eset <- function(summarizedExperiment, assayName="counts"){
	expressionSet <- ExpressionSet(
	  assayData=assays(summarizedExperiment)[[assayName]],
	  phenoData=
	    AnnotatedDataFrame(as.data.frame(colData(summarizedExperiment))))
	return(expressionSet)
}

#' get_eset_from_matrix
#' 
#' Makes an ExpressionSet from a matrix.
#'
#' @param inputMatrix User-specified expression matrix.
#' @param batchVariable Name of the batch variable.
#' 
#' @returns ExpressionSet.
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' 
#' @examples
#' exampleList <- getDeconvolutionExampleData()
#'
#' @export
get_eset_from_matrix <- function(inputMatrix, batchVariable="SampleName"){
  phenoData <- data.frame(new.variable=colnames(inputMatrix))
  colnames(phenoData) <- batchVariable
  rownames(phenoData) <- colnames(inputMatrix)
  expressionSet <- Biobase::ExpressionSet(
    assayData=inputMatrix, phenoData=Biobase::AnnotatedDataFrame(phenoData))
  return(expressionSet)
}