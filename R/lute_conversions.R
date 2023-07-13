#!/usr/bin/env R

#
# Functions to convert between data object classes.
#

sce_to_eset <- function(sce){}

eset_to_sce <- function(eset, assay.name = "counts"){
	sce.new <- SingleCellExperiment(assays = list(assay.name = exprs(eset)))
	colData(sce.new) <- DataFrame(as.matrix(pData(eset)))
	metadata(sce.new) <- metadata(eset)
	return(sce.new)
}

sce_to_se <- function(sce){}

se_to_sce <- function(se){}

eset_to_se <- function(eset){}

se_to_eset <- function(se){}