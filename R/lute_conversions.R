#!/usr/bin/env R

#
# Functions to convert between data object classes.
#

#' sce_to_eset
#' Convert SingleCellExperiment to ExpressionSet.
#' @param sce SingleCellExperiment.
#' @param assay.name Name of assay to store in new eset.
#' @returns ExpressionSet.
#' @examples
#' sce <- random_sce()
#' sce_to_eset(sce, "counts")
#' @export
sce_to_eset <- function(sce, assay.name = "counts"){
	eset <- ExpressionSet(assayData = assays(sce)[[assay.name]])
	pData(eset) <- as.data.frame(colData(sce))
	return(eset)
}

#' eset_to_sce
#' Convert ExpressionSet to SingleCellExperiment.
#' @param eset ExpressionSet.
#' @param assay.name Name of new assay in SingleCellExperiment.
#' @returns ExpressionSet.
#' @examples
#' eset <- lute:::.get_decon_example_data_bisque()$sc.eset
#' eset_to_sce(eset)
#' @export
eset_to_sce <- function(eset, assay.name = "counts"){
	sce.new <- SingleCellExperiment(assays = list(assay.name = exprs(eset)))
	colData(sce.new) <- DataFrame(as.matrix(pData(eset)))
	return(sce.new)
}

#' sce_to_se
#' Convert SingleCellExperiment to SummarizedExperiment.
#' @param sce SingleCellExperiment.
#' @returns SummarizedExperiment.
#' @examples
#' sce <- random_sce()
#' sce_to_se(sce)
#' @export
sce_to_se <- function(sce){
	se <- SummarizedExperiment(assays = assays(sce))
	colData(se) <- colData(sce)
	metadata(se) <- metadata(sce)
	return(se)
}

#' se_to_sce
#' Convert SummarizedExperiment to SingleCellExperiment.
#' @param se SummarizedExperiment.
#' @returns SingleCellExperiment.
#' @examples
#' se_to_sce(SummarizedExperiment())
#' @export
se_to_sce <- function(se){
	sce <- SingleCellExperiment(assays = assays(se))
	colData(sce) <- colData(se)
	metadata(sce) <- metadata(se)
	return(sce)
}

#' eset_to_se
#' Convert ExpressionSet to SummarizedExperiment.
#' @param eset ExpressionSet
#' @param assay.name Name of assay to store in new eset.
#' @returns SummarizedExperiment
#' @examples
#' eset <- lute:::.get_decon_example_data_bisque()$sc.eset
#' eset_to_se(eset, "counts")
#' @export
eset_to_se <- function(eset, assay.name = "counts"){
	se.new <- SummarizedExperiment(assays = list(assay.name = exprs(eset)))
	colData(se.new) <- DataFrame(as.matrix(pData(eset)))
	return(se.new)
}

#' se_to_eset
#' Convert SummarizedExperiment to ExpressionSet.
#' @param se SummarizedExperiment.
#' @param assay.name Name of assay to store in new eset.
#' @returns ExpressionSet.
#' @examples
#' se <- sce_to_se(random_sce())
#' se_to_eset(se)
#' @export
se_to_eset <- function(se, assay.name = "counts"){
	eset <- ExpressionSet(assayData = assays(se)[[assay.name]])
	pData(eset) <- as.data.frame(colData(se))
	return(eset)
}