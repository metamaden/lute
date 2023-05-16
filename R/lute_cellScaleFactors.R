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
#'
#' @examples
#' csf_filter(c("neuron", "glial"))
#'
csf_filter <- function(unique.cell.types, csf.ref = NULL, 
                       prefer.orthogonal = TRUE, summarize = "median"){
  if(is(csf.ref, "NULL")){csf.ref <- get_csf_reference()}
  
  # filter available cell type labels
  which.cell.types <- csf.ref$cell_type %in% unique.cell.types
  csf.ref.out <- csf.ref[which.cell.types,]
  # filter available assay types
  if(prefer.orthogonal){
    data.source.vector <- unique(csf.ref.out$scale.factor.data.source)
    orthogonal.sources <- data.source.vector[!data.source.vector == "expression"]
    if(length(orthogonal.sources) == 0){
      message("No orthogonal data sources. Returning expression scale factors.")
    } else{
      is.orthogonal <- csf.ref.out$scale.factor.data.source %in% orthogonal.sources
      csf.ref.out <- csf.ref.out[is.orthogonal,]
    }
  }
  return(csf.ref.out)
}



