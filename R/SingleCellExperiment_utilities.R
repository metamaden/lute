#!/usr/bin/env R

# Author: Sean Maden
#
# Utilities to run deconvolution simulations, with support for 
# SingleCellExperiment objects.
#
#
#

#' kexpr_sce
#'
#' Get expression summaries by types from SingleCellExperiment or 
#' SummarizedExperiment objects.
#'
#' @param sce Either a SummarizedExperiment or SingleCellExperiment object. 
#' @param return.lgv Whether to return kexpr summaries as list. If False, returns
#' as a 2d matrix.
#' @param expr.scale Name of expression data obtained from assays(sce), e.g.
#' logcounts, counts, etc.
#' @param marker.varname Marker variable name from rowData(sce).
#' @param type.varname Type variable name from colData(sce).
#' @param marker.summary Summary function to use, e.g. mean, median, etc.
#' @param vebose Whether to return verbose status messages.
#' @returns Either a Z signature/reference matrix, or lgv list of marker expression
#' by type
#' @examples
#' # example
#' @seealso decon_analysis,
#' @export
kexpr_sce <- function(sce, return.lgv = TRUE, expr.scale = "logcounts", 
                      marker.varname = "marker", type.varname = "cellType", 
                      marker.summary = "mean", verbose = FALSE){
  # check sce type
  cond.sce <- is(sce, "SingleCellExperiment")|is(sce, "SummarizedExperiment")
  if(!cond.sce){
    stop("Error, sce not a SingleCellExperiment or SummarizedExperiment.")}
  # check marker variable
  cond.marker <- marker.varname %in% colnames(rowData(sce)); scef <- sce
  if(cond.marker){
    which.marker <- which(rowData(scef)[,marker.varname])
    scef <- scef[which.marker,]
  } else{
    if(verbose){message("marker.varname not in rowData; using all genes.")}
  }
  # check scale
  cond.scale <- expr.scale %in% names(assays(scef))
  if(!cond.scale){stop("Error, cond.scale not in sce assays.")}
  me <- eval(parse(text = paste0("assays(scef)$",expr.scale)))
  # check type variable
  cd <- colData(scef); cond.type <- type.varname %in% colnames(cd)
  if(!cond.type){stop("Error, type.varname not in colData.")}
  typevar <- cd[,type.varname]; typev <- unique(typevar)
  # parse marker summary options
  if(marker.summary == "mean"){
    if(verbose){message("Getting marker mean expression by type...")}
    z <- do.call(cbind, lapply(typev, function(typei){
      DelayedArray::rowMeans(me[,typevar==typei])
    }))
  } else if(marker.summary == "median"){
    if(verbose){message("Getting marker median expression by type...")}
    z <- do.call(cbind, lapply(typev, function(typei){
      matrixStats::rowMedians(me[,typevar==typei])
    }))
  } else{stop("Error, invalid marker.summary option.")}
  colnames(z) <- typev
  if(return.lgv){
    if(verbose){message("Returning lgv...")}
    lgv <- lapply(seq(ncol(z)), function(ci){z[,ci]})
    names(lgv) <- colnames(z); return(lgv)
  }
  if(verbose){message("Returning Z reference...")}; return(z)
}