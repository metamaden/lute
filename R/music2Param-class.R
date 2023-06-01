#!/usr/bin/env R

# Author: Sean Maden

#' music2Param-class
#' 
#' Main constructor for class to manage mappings to the deconvolution method 
#' function \code{MuSiC::music2_prop()}.
#' 
#' @include lute_generics.R
#' @include independentbulkParam-class.R
#' 
#' @details Main constructor for class \linkS4class{music2Param}.
#' 
#' @rdname music2Param-class
#' 
#' @seealso 
#' \linkS4class{deconParam}, \linkS4class{referencebasedParam}, 
#' \linkS4class{independentbulkParam}
#'
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
#' @aliases 
#' music2Param-class, MuSiC2Param-class, Music2Param-class
#' 
#' @references 
#' 
#' Fan, Jiaxin. MuSiC2: MuSiC2: cell type deconvolution for multi-condition bulk 
#' RNA-seq data. (2023) GitHub, R package version 0.1.0. URL: 
#' https://github.com/Jiaxin-Fan/MuSiC2
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#' 
#' Jiaxin Fan, Yafei Lyu, Qihuang Zhang, Xuran Wang, Mingyao Li, Rui Xiao, 
#' MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data, 
#' Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, bbac430, 
#' https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbac430
#'
setClass("music2Param", contains="independentbulkParam", slots=c(
	sc.eset = "ExpressionSet", assay.name = "character", batch.variable = "character", 
	celltype.variable = "character", iter.max = "numeric", nu = "numeric", epsilon = "numeric", 
	truep = "numeric", ct.cell.size = "numeric", celltype.subset = "character"))

#' Make an object of class music2Param
#'
#' Main function to make a new object of class \linkS4class{music2Param}, with
#' defaults for required arguments.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param yi Bulk mixed signals matrix of independent samples, which should not overlap samples in y.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param y.eset ExpressionSet of bulk mixed expression signals.
#' @param sc.eset ExpressionSet of single-cell transcriptomics data.
#' @param assay.name Expression data type (e.g. counts, logcounts, tpm, etc.).
#' @param batch.variable Name of variable identifying the batches in sc.eset pData/coldata.
#' @param celltype.variable Name of cell type labels variable in sc.eset pData/coldata.
#' @param condition.variable Name of variable in y.eset and sc.eset containing condition labels.
#' @param control.label Label of control condition samples in condition variable.
#' @param case.label Label of case condition samples in condition variable.
#' @param method.type Name of method source library to call for music2_prop; either "MuSiC" or "MuSiC2".
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @seealso \linkS4class{musicParam}
#'
#' @returns Object of class \linkS4class{music2Param}.
#' 
#' @export
music2Param <- function(y = NULL, yi = NULL, z = NULL, s = NULL, y.eset = NULL, sc.eset = NULL, 
	assay.name = "counts", batch.variable = "SubjectName", celltype.variable = "cellType", 
	condition.variable = "experiment.condition", control.label = "control", case.label = "case",
	method.type = "MuSiC2", return.info = FALSE) {

  new("music2Param", y = y, yi = yi, z = z, s = s, y.eset = y.eset, sc.eset = sc.eset, 
	assay.name = assay.name, batch.variable = batch.variable, celltype.variable = celltype.variable,
	condition.variable = condition.variable, control.label = control.label, case.label = case.label,
	method.type = method.type, return.info = return.info)
}

#' Deconvolution method for \linkS4class{music2Param}
#'
#' Main method to access the MuSiC2 deconvolution method from the main lute 
#' deconvolution genetic.
#'
#' @details Takes an object of class \linkS4class{music2Param} as input, 
#' returning a list or vector of predicted cell type proportions.
#'
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#' 
#' @references 
#' 
#' Fan, Jiaxin. MuSiC2: MuSiC2: cell type deconvolution for multi-condition bulk 
#' RNA-seq data. (2023) GitHub, R package version 0.1.0. URL: 
#' https://github.com/Jiaxin-Fan/MuSiC2
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#' 
#' Jiaxin Fan, Yafei Lyu, Qihuang Zhang, Xuran Wang, Mingyao Li, Rui Xiao, 
#' MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data, 
#' Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, bbac430, 
#' https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbac430
#'
#' @export
setMethod("deconvolution", signature(object = "music2Param"), function(object){
  require(Biobase)
  # load data
  lparam <- callNextMethod()
  object <- lparam[["object"]]

  # instantiate function objects
  y.eset <- lparam[["y.eset"]]
  sc.eset <- object@sc.eset
  celltype.subset <- object@celltype.subset
  batch.variable <- object@batch.variable
  iter.max <- object@iter.max
  nu <- object@nu
  epsilon <- object@epsilon
  truep <- object@truep
  s <- object@s
  weight.basis <- object@weight.basis
  transform.bisque <- object@transform.bisque

  # get result according to method type
  source.libarary <- object@method.type
  if(source.library %in% c("music", "MuSiC", "Music", "MUSIC")){
  	message("Using the MuSiC implementation of music2_prop()..."); require(MuSiC)
  	result <- MuSiC::music2_prop(bulk.control.mtx = y, bulk.case.mtx = yi, sc.sce = sc.sce, 
  		clusters = celltype.variable, samples = batch.variable, cell_size = cell_size, select.ct = NULL)
  } else{
  	message("Using the MuSiC2 implementation of music2_prop()..."); require(MuSiC2)
	result <- MuSiC2::music2_prop(bulk.eset = y.eset, sc.eset = sc.sce, condition = condition.variable,
	                      control = control.label, case = case.label, clusters = celltype.variable,
	                      samples = batch.variable, cell_size = cell_size, select.ct = unique.types)
  }
  # return results
  lr <- predictions <- result$bulk.props
  if(object[["return.info"]]){
    lr <- list(predictions = predictions, result.info = result, 
               metadata = list(lmd = lparam[["metadata"]], 
                y.eset = y.eset, sc.eset = sc.eset))}
  return(lr)
})