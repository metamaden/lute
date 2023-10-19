#!/usr/bin/env R

### Author: Sean Maden
###
### Parses options to get cell scale factors from reference (package cellScaleFactors).
###

#' get_csf_reference
#'
#' Retrieves the cell scale factors (csf) reference from the cellScaleFactors package.
#' 
#' @param userCellTypesVector Vector of user-specified cell types.
#' @param preferOrthogonal Whether to prefer expression-orthogonal values (if 
#' TRUE, removes expression-based values, but only if alternative value types 
#' are available).
#' 
#' @details Returns a table of cell scale factors from various data sources. The 
#' cell scale factors reference table has the following columns:
#' 
#' 1. cell_type : Label of the cell type for the scale factor (e.g. neuron, T cell, etc.)
#' 2. tissue : Label of the tissue of origin (e.g. brain, blood, etc.)
#' 3. scale.factor.value : Point scale factor value prior to additional normalization
#' 4. scale.factor.type : Label for scale factor type (e.g. cell or nuclear area, etc.)
#' 5. scale.factor.data.source : Label for scale factor source (e.g. osmFISH, 
#' housekeeping gene expression, etc.)
#' 6. citation.s : Citation(s) of source studies from which original measures or 
#' measure summaries were made.
#' 
#' Further details about the reference table can be found in the cellScaleFactors package.
#' 
#' @importFrom methods is
#' 
#' @examples
#' example.data <- get_decon_example_data()
#' 
#' @returns Table of type "data.frame" or "tibble".
#' @export
get_csf_reference <- function(userCellTypesVector=NULL, preferOrthogonal=TRUE){
  referenceReturn <- load_csf_rda()
  if(preferOrthogonal){
    dataSourceVector <- unique(referenceReturn$scale.factor.data.source)
    orthogonalSourcesVector <- 
      dataSourceVector[!dataSourceVector == "expression"]
    if(length(orthogonalSourcesVector) == 0){
      message("No orthogonal data sources. Returning expression scale factors.")
    } else{
      isOrthogonal <- 
        referenceReturn$scale.factor.data.source %in% orthogonalSourcesVector
      referenceReturn <- referenceReturn[isOrthogonal,]
    }
  }
  if(!is(userCellTypesVector, "NULL")){
    referenceReturn <- csf_filter_labels(userCellTypesVector, referenceReturn)
  }
  return(referenceReturn)
}


#'
load_csf_rda <- function(){
  path <- system.file(
    file.path("rda", "cellScaleFactors.rda"), package="cellScaleFactors")
  return(get(load(path)))
}

#'
csf_filter_labels <- function(labels, reference=NULL){
  if(is(reference, "NULL")){reference <- get_csf_reference()}
  referenceLabels <- reference$cell_type
  dfReturn <- do.call(rbind, lapply(labels, function(label){
    filter1 <- grepl(label, referenceLabels)
    filter2 <- grepl(toupper(label), toupper(referenceLabels))
    filter3 <- grepl(label, gsub(" ", "",referenceLabels))
    filter4 <- grepl(label, gsub(" ", "",referenceLabels))
    filter5 <- grepl(gsub(" ", "", toupper(label)), 
                     gsub(" ", "", toupper(referenceLabels)))
    labelFilter <- filter1|filter2|filter3|filter4|filter5
    reference[labelFilter,,drop=FALSE]
  }))
  dfReturn <- as.data.frame(dfReturn)
  dfReturn
}
