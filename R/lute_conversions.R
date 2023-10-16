#!/usr/bin/env R

###
### Functions to convert between data object classes.
###

#' sce_to_eset
#' Convert SingleCellExperiment to ExpressionSet.
#' @param sce Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @param assay.name Name of assay to store in new eset.
#' @returns ExpressionSet.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase pData
#' @importFrom Biobase AnnotatedDataFrame
#'
#' @examples
#' sce <- random_sce()
#' sce_to_eset(sce, "counts")
#' 
#' @export
sce_to_eset <- function(sce, assay.name="counts"){
	eset <- ExpressionSet(assayData=assays(sce)[[assay.name]],
	                      phenoData=AnnotatedDataFrame(
	                        as.data.frame(colData(sce))))
	return(eset)
}

#' eset_to_sce
#' Convert ExpressionSet to SingleCellExperiment.
#' @param eset Object of type ExpressionSet (see \code{?ExpressionSet}).
#' @param assay.name Name of new assay in new SingleCellExperiment object.
#' @returns ExpressionSet.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase pData
#' @importFrom Biobase exprs
#' @importFrom SingleCellExperiment SingleCellExperiment
#' 
#' @examples
#' eset <- get_decon_example_data_bisque()$sc.eset
#' eset_to_sce(eset)
#' 
#' @export
eset_to_sce <- function(eset, assay.name="counts"){
  assays.list <- list(assay.name=exprs(eset))
  names(assays.list) <- assay.name
	sce.new <- SingleCellExperiment(assays=assays.list,
	                                colData=DataFrame(as.matrix(pData(eset))))
	return(sce.new)
}

#' sce_to_se
#' Convert SingleCellExperiment to SummarizedExperiment.
#' @param sce Object of type SingleCellExperiment (see 
#' \code{?SingleCellExperiment}).
#' @returns SummarizedExperiment.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors metadata
#' 
#' @examples
#' sce <- random_sce()
#' sce_to_se(sce)
#' @export
sce_to_se <- function(sce){
	se <- SummarizedExperiment(assays=assays(sce),
	                           colData=colData(sce),
	                           metadata=metadata(sce))
	return(se)
}

#' se_to_sce
#' Convert SummarizedExperiment to SingleCellExperiment.
#' @param se Object of type SummarizedExperiment (see 
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
se_to_sce <- function(se){
	sce <- SingleCellExperiment(assays=assays(se),
	                            colData=colData(se),
	                            metadata=metadata(se))
	return(sce)
}

#' eset_to_se
#' Convert ExpressionSet to SummarizedExperiment.
#' @param eset Object of type ExpressionSet (see \code{?ExpressionSet}).
#' @param assay.name Name of assay to store in new SummarizedExperiment object.
#' @returns New object of type SummarizedExperiment.
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase pData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' 
#' @examples
#' eset <- get_decon_example_data_bisque()$sc.eset
#' eset_to_se(eset, "counts")
#' 
#' @export
eset_to_se <- function(eset, assay.name="counts"){
  assays.list <- list(assay.name=exprs(eset))
  names(assays.list) <- assay.name
	se.new <- SummarizedExperiment(assays=assays.list,
	                               colData=DataFrame(as.matrix(pData(eset))))
	return(se.new)
}

#' se_to_eset
#' Convert SummarizedExperiment to ExpressionSet.
#' @param se Object of type SummarizedExperiment (see 
#' \code{?SummarizedExperiment}).
#' @param assay.name Name of assay to store in new ExpressionSet object.
#' @returns New object of type ExpressionSet.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom Biobase pData
#' @importFrom Biobase ExpressionSet
#' @importFrom Biobase AnnotatedDataFrame
#'
#' @examples
#' se <- sce_to_se(random_sce())
#' se_to_eset(se)
#' 
#' @export
se_to_eset <- function(se, assay.name="counts"){
	eset <- ExpressionSet(assayData=assays(se)[[assay.name]],
	                      phenoData=AnnotatedDataFrame(
	                        as.data.frame(colData(se))))
	return(eset)
}

#' get_eset_from_matrix
#' 
#' Makes an ExpressionSet from a matrix.
#'
#' @param mat User-specified expression matrix.
#' @param batch.variable Name of the batch variable.
#' @returns ExpressionSet.
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' 
#' @examples
#' example.data <- get_decon_example_data()
#'
#' @export
get_eset_from_matrix <- function(mat, batch.variable="SampleName"){
  pdata <- data.frame(new.variable=colnames(mat))
  colnames(pdata) <- batch.variable
  rownames(pdata) <- colnames(mat)
  eset <- Biobase::ExpressionSet(assayData=mat, phenoData=Biobase::AnnotatedDataFrame(pdata))
  return(eset)
}