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

#' filter_na_cells
#'
#' Filter cells on NA frequency.
#'
#' @param sce A SingleCellExperiment object.
#' @param type.variable Name of the type variable in sce colData.
#' @param na.freq The upper threshold for tolerated NA frequency. Cells with NA 
#' frequency exceeding this are excluded.
#' @param assayname Name of assay in assays(sce) to check for NAs.
#' @param append.metadata Whether to append metadata summary of preprocessing 
#' parameters and results.
#' @param verbose Whether to include verbose status messages.
#' @returns scef, filtered SingleCellExperiment with optional appended 
#' preprocessing metadata.
#' @examples 
#' sce <- random_sce(na.include = T, na.fract = 0.4)
#' scef <- filter_na_bytype(sce)
#'
#'
filter_na_cells <- function(sce, type.variable = "celltype", na.freq = 0.25,
                          assayname = "counts", append.metadata = TRUE,
                          verbose = FALSE){
  if(!is(sce, "SingleCellExperiment")){
    stop("Error, sce must be a SingleCellExperiment.")}
  cd <- colData(sce)
  if(!type.variable %in% colnames(cd)){
    stop("Error, type.variable must be a column name in sce colData.")
  }
  typev <- unique(cd[,typevar])
  mexpr <- assays(sce)[[assayname]]
  # count nas
  na.freqv <- apply(mexpr, 2, function(ci){length(which(is.na(ci)))})
  
  
  
  # append metadata
}

#' filter_na_type
#'
#'
filter_na_type <- function(){
  
}

#' filter_na_sce
#'
#'
filter_na_sce <- function(){
  
}

#--------------
# normalization
#--------------

#-----------------
# bias adjustments
#-----------------