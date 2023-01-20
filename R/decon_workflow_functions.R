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
#' @param verbose Whether to include verbose status messages.
#' @returns scef, filtered SingleCellExperiment
#' @examples 
#' sce <- random_sce(na.include = T, na.fract = 0.4)
#' scef <- filter_na_bytype(sce)
#'
#'
filter_na_bytype <- function(sce, type.variable = "celltype", verbose = FALSE){
  
}

#--------------
# normalization
#--------------

#-----------------
# bias adjustments
#-----------------