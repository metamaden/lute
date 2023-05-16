#!/usr/bin/env R

# Author: Sean Maden
#
# Parses options to get cell scale factors from reference (package cellScaleFactors).
#

#'
#'
#'
get_csf_reference <- function(){
  require(cellScaleFactors)
  path <- system.file(
    file.path("rda", "cellScaleFactors.rda"), 
    package = "cellScaleFactors")
  csf <- get(load(path))
  return(csf)
}

#'
#'
#'
csf_options <- function(user.s, cell.types.vector = NULL){
  if(is(user.s, "character")){
    # s is an option to retrieve reference data
    csf.ref <- get_csf_reference()
    csf.return <- csf_filter(csf.ref, cell.types.vector)
  } else if(is(user.s, "numeric")|is(user.s, "matrix")){
    # s contains user-specified values. continue
  }
}

#'
#'
#'
csf_filter <- function(csf.ref, cell.types.vector, prefer.orthogonal = TRUE){
  # filter labels
  csf.ref.out <- csf.ref
  if(is(cell.types.vector, NULL)){
    return(csf.ref)
  } else{
    which.cell.types <- csf.ref$cell_type %in% cell.types.vector
    csf.ref.out <- csf.ref[which.cell.types,]
  }
  # filter available assay types
  if(prefer.orthogonal){
    assay.types <- unique(csf.ref.out)
  }
  
}

