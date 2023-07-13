#!/usr/bin/env R

#
# Functions to convert between data object classes.
#



sce_to_eset <- function(sce, assay.name = "counts"){
	eset <- ExpressionSet(assayData = assays(sce)[[assay.name]])
	pData(eset) <- as.data.frame(colData(sce))
	metadata(eset) <- metadata(sce)
	return(eset)
}

eset_to_sce <- function(eset, assay.name = "counts"){
	sce.new <- SingleCellExperiment(assays = list(assay.name = exprs(eset)))
	colData(sce.new) <- DataFrame(as.matrix(pData(eset)))
	metadata(sce.new) <- metadata(eset)
	return(sce.new)
}

sce_to_se <- function(sce){
	se <- SummarizedExperiment(assays = assays(sce))
	colData(se) <- colData(sce)
	metadata(se) <- metadata(sce)
	return(se)
}

se_to_sce <- function(se){
	sce <- SingleCellExperiment(assays = assays(se))
	colData(sce) <- colData(se)
	metadata(sce) <- metadata(se)
	return(sce)
}

eset_to_se <- function(eset, assay.name = "counts"){
	se.new <- SummarizedExperiment(assays = list(assay.name = exprs(sc.eset)))
	colData(se.new) <- DataFrame(as.matrix(pData(sc.eset)))
	return(se.new)
}

se_to_eset <- function(se, assay.name = "counts"){
	eset <- ExpressionSet(assayData = assays(se)[[assay.name]])
	pData(eset) <- as.data.frame(colData(se))
	return(eset)
}