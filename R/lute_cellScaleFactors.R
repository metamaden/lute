#!/usr/bin/env R

# Author: Sean Maden
#
# Parses options to get cell scale factors from reference (package cellScaleFactors).
#

#' get_csf_reference
#'
#' Retrieves the cell scale factors (csf) reference from the cellScaleFactors package.
#' 
#' @param user.celltypes.vector Vector of user-specified cell types.
#' @param prefer.orthogonal Whether to prefer expression-orthogonal values (if 
#' TRUE, removes expression-based values).
#' 
#' @details Returns a table of cell scale factors from various data sources.
#' 
#' Table has the following columns:
#' 
#' 1. cell_type : Label of the cell type for the scale factor ()
#' 2. tissue : Label of the tissue of origin (e.g. brain, blood, etc.)
#' 3. scale.factor.value : Point scale factor value prior to additional normalization
#' 4. scale.factor.type : Description of the type of scale factor (e.g. cell or nuclear area, etc.)
#' 
#' 
#' @returns Table of type "data.frame" or "tibble".
#'
#'
get_csf_reference <- function(user.celltypes.vector = NULL, 
                              prefer.orthogonal = TRUE){
  ref <- load_csf_rda()
  if(prefer.orthogonal){
    data.source.vector <- unique(ref$scale.factor.data.source)
    orthogonal.sources.vector <- 
      data.source.vector[!data.source.vector == "expression"]
    if(length(orthogonal.sources.vector) == 0){
      message("No orthogonal data sources. Returning expression scale factors.")
    } else{
      is.orthogonal <- 
        ref$scale.factor.data.source %in% orthogonal.sources.vector
      ref <- ref[is.orthogonal,]
    }
  }
  if(!is(user.celltypes.vector, "NULL")){
    ref <- csf_filter_labels(user.celltypes.vector, ref)
  }
  return(ref)
}

#'
#'
#'
load_csf_rda <- function(){
  require(cellScaleFactors)
  path <- system.file(
    file.path("rda", "cellScaleFactors.rda"), 
    package = "cellScaleFactors")
  return(get(load(path)))
}

#'
#'
#'
#'
csf_filter_labels <- function(labels, reference = NULL){
  require(dplyr)
  if(is(reference, "NULL")){reference <- get_csf_reference()}
  reference.labels <- reference$cell_type
  do.call(rbind, lapply(labels, function(label){
    filter1 <- grepl(label, reference.labels)
    filter2 <- grepl(toupper(label), toupper(reference.labels))
    filter3 <- grepl(label, gsub(" ", "",reference.labels))
    filter4 <- grepl(label, gsub(" ", "",reference.labels))
    filter5 <- grepl(gsub(" ", "", toupper(label)), 
                     gsub(" ", "", toupper(reference.labels)))
    label.filter <- filter1|filter2|filter3|filter4|filter5
    reference[label.filter,,drop=F]
  })) %>% as.data.frame()
}
