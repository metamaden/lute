#!/usr/bin/env R

# Author: Sean Maden
#
# Defines functions used by NextFlow processes for lute workflow.
#
#
#
#

#-------------------------
# preprocessing -- filters
#-------------------------

#' filter_na_bytype
#'
#' In some SingleCellExperiment, filter the NA frequency by type.
#'
#' @param sce A SingleCellExperiment object.
#' @param type.variable Name of the type variable in sce colData.
#' @param na.freq The upper threshold for tolerated NA frequency. Cells with NA 
#' frequency exceeding this are excluded.
#' @param assayname Name of assay in assays(sce) to check for NAs.
#' @param verbose Whether to include verbose status messages.
#' @returns scef, filtered SingleCellExperiment with optional appended 
#' preprocessing metadata.
#' @examples 
#' sce <- random_sce(na.include = T, na.fract = 0.4)
#' scef <- filter_na_bytype(sce)
#'
#'
filter_na_bytype <- function(sce, type.variable = "celltype", na.freq = 0.25, 
                             assayname = "counts", verbose = FALSE){
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  cd <- colData(sce)
  if(!type.variable %in% colnames(cd)){
    stop("Error, type.variable must be a column name in sce colData.")
  }
  typev <- unique(cd[,typevar])
  
  
  
  # append metadata
}

#--------------
# normalization
#--------------

#-----------------
# bias adjustments
#-----------------