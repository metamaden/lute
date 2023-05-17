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
#' csf_filter_assays(c("neuron", "glial"))
#'
csf_filter_assays <- function(unique.cell.types, csf.ref = NULL,
                              prefer.orthogonal = TRUE, 
                              summarize = "median"){
  if(is(csf.ref, "NULL")){csf.ref <- get_csf_reference()}
  csf.ref.out <- csf_filter_labels(labels = unique.cell.types)
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

#'
#'
#'
#'
csf_filter_labels <- function(labels, reference = NULL){
  require(dplyr)
  if(is(reference, "NULL")){reference <- get_csf_reference()}
  reference.labels <- reference$cell_type
  do.call(rbind, lapply(labels, function(label){
    filter1 <- grepl(paste0("^", label, "$"), reference.labels)
    filter2 <- grepl(paste0("^", toupper(label), "$"), toupper(reference.labels))
    filter3 <- grepl(paste0("^", label, "$"), gsub(" ", "",reference.labels))
    filter4 <- grepl(paste0("^", gsub(" ", "", label), "$"), gsub(" ", "",reference.labels))
    filter5 <- grepl(paste0("^", gsub(" ", "", toupper(label)), "$"), 
                     gsub(" ", "", toupper(reference.labels)))
    label.filter <- filter1|filter2|filter3|filter4|filter5
    reference[label.filter,,drop=F]
  })) %>% as.data.frame()
}

#'
#'
#'
csf_median_norm_scales <- function(csf.reference, 
                                   seed.cell.type = "neuron", 
                                   unique.cell.types = c("neuron", "glial"),
                                   prefer.orthogonal = TRUE,
                                   summary.statistic = "median",
                                   return.type = c("data.frame", "vector")){
  require(dplyr)
  ref.filter <- csf.reference[csf.reference$cell_type %in% unique.cell.types,]
  if(prefer.orthogonal){ref.filter <- filter_reference_orthogonal(ref.filter)}
  
  # get summary set
  category.vector <- paste0(ref.filter$scale.factor.type, 
                            ";", ref.filter$scale.factor.data.source)
  list.norm.factors <- lapply(category.vector, function(category){
    ref.filter.iter <- ref.filter[ref.filter,]
  })
  
}

#'
#'
filter_reference_orthogonal <- function(ref, expression.string = "expression"){
  category.vector <- paste0(ref$scale.factor.type, ";", ref$scale.factor.data.source)
  unique.categories <- unique(category.vector)
  which.orthogonal <- !grepl(
    paste0(".*;",expression.string,"$"), unique.categories)
  unique.categories.ortho <- unique.categories[which.orthogonal]
  if(length(unique.categories.ortho) > 0){
    unique.categories.out <- category.vector[which.orthogonal]
  } else{
    unique.categories.out <- unique.categories
  }
  filter.ref <- category.vector %in% unique.categories.out
  return(ref[filter.ref,])
}

#'
#'
get_norm_factors <- function(ref.filter){
  
  
}