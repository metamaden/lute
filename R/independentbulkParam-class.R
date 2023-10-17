#!/usr/bin/env R

### Author: Sean Maden

#' independentbulkParam-class
#'
#' Class and methods for managing methods requiring independent bulk samples.
#'
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#'
#' @param yi Mixed signals matrix from bulk samples, independent from primary 
#' mixed signals matrix y.
#'
#' @details The main purpose of this class is to compare bulk sample data 
#' between the passed objects y and yi. Since we assume yi contains the 
#' independent bulk samples, it should not have overlapping sample IDs 
#' (colnames), and it should have overlapping marker IDs (rownames) compared to 
#' the reference bulk samples y.
#'
#' @seealso \linkS4class{deconParam}, \linkS4class{referencebasedParam}
#' 
#' @examples 
#' new("independentbulkParam")
#' 
#' @returns New object.
setClass("independentbulkParam", contains="referencebasedParam", 
         slots=c(yi="matrix"))

#' Make a new \linkS4class{independentbulkParam} object
#' 
#' Function to make a new object of class \linkS4class{independentbulkParam}
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples.
#' @param yi Bulk mixed signals matrix of independent samples, which should not 
#' overlap samples in y.
#' @param z Signature matrix of cell type-specific signals. If not provided, can 
#' be computed from a provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types 
#' to deconvolve.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @examples 
#' new("independentbulkParam")
#'
#' @returns New object.
#'
#' @export
independentbulkParam <- function(y=NULL, yi=NULL, z=NULL, s=NULL, 
                                 return.info=FALSE) {
  input_y <- y; input_yi <- yi; input_z <- z; input_s <- s
    if(is(input_y, "NULL")){input_y <- matrix(0)}
    if(is(input_z, "NULL")){input_z <- matrix(0)}
    if(is(input_yi, "NULL")){input_yi <- matrix(0)}
    if(is(input_s, "NULL")){input_s <- rep(1, ncol(input_z))}
    param <- new("independentbulkParam", 
                 y=input_y, yi=input_yi, z=input_z, s=input_s, 
                 return.info=return.info)
    return(param)
}

#' Deconvolution method for class \linkS4class{independentbulkParam}
#'
#' Function to perform standard operations prior to deconvolution (a.k.a. 
#' "deconvolution prep") for an object of class 
#' \linkS4class{independentbulkParam}.
#'
#' @param object An object of class \linkS4class{independentbulkParam}.
#'
#' @details Takes an object of \linkS4class{independentbulkParam} class as 
#' input, and returns a list with the filtered/checked/parsed experiment objects.
#' 
#' @examples 
#' new("independentbulkParam")
#'
#' @returns Method results.
#'
#' @export
setMethod("deconvolution", "independentbulkParam", function(object) {
    lparam <- callNextMethod()
    unique.marker.labels <- unique.sample.labels <- NULL
    overlapping.marker.labels <- overlapping.sample.labels <- NULL
    input_y <- lparam[["y"]]; input_yi <- lparam[["yi"]] # get bulk data
    markers.y <- rownames(input_y)
    markers.yi <- rownames(input_yi) # parse bulk marker IDs
    
    ## compare marker labels and subset yi on overlapping markers
    if(is(markers.y, "NULL")|is(markers.yi, "NULL")){
        message("Warning, no marker labels found in either y or yi.")
    } else{
        unique.marker.labels <- unique(markers.y, markers.yi)
        overlapping.marker.labels <- intersect(markers.y, markers.yi)
        if(length(overlapping.marker.labels) > 0){
          input_yi <- input_yi[overlapping.marker.labels,]
        }
    }
  
    ## compare sample labels and remove overlapping samples
    samples.y <- colnames(input_y)
    samples.yi <- colnames(input_yi) # parse bulk sample IDs
    ## compare sample IDs
    if(is(samples.y, "NULL")|is(samples.yi, "NULL")){
        message("Warning, no sample labels found in either y or yi.")
    } else{
        unique.sample.labels <- unique(samples.y, samples.yi)
        overlapping.sample.labels <- intersect(samples.y, samples.yi)
        if(length(overlapping.samples) > 0){
            filter <- !colnames(input_yi) %in% overlapping.sample.labels
            input_yi <- input_yi[, filter, drop=FALSE]
        }
    }
  
    ## parse return list
    ## get metadata to return
    list_metadata <- list(
      unique.marker.labels=unique.marker.labels,
                unique.sample.labels=unique.sample.labels,
                overlapping.marker.labels=overlapping.marker.labels,
                overlapping.sample.labels=overlapping.sample.labels)
    return_list <- list(y=input_y, yi=input_yi, object=object, 
                        metadata=list_metadata)
    return(return_list)
})

#' Method for \linkS4class{independentbulkParam}
#'
#' @param object An object of class \linkS4class{independentbulkParam} (see 
#' \code{?independentbulkParam}).
#' @details Display data summaries for an object of class 
#' \linkS4class{independentbulkParam}.
#' 
#' @examples 
#' new("independentbulkParam")
#'
#' @returns Shows object summaries.
#'
#' @export
setMethod("show", "independentbulkParam", function(object) {
  input_y <- object[["y"]]; input_yi <- object[["yi"]] # get bulk data
  samples.y <- colnames(input_y); samples.yi <- colnames(input_yi) # get samples
  markers.y <- rownames(input_y); markers.yi <- rownames(input_yi) # get markers
  ## print info summaries
  message("Summary of independentbulkParam data:")
  message("\tNumber of unique sample IDs : ", 
          length(unique(markers.y, markers.yi)), "\n")
  message("\tNumber of unique marker IDs : ", 
          length(unique(samples.y, samples.yi)), "\n")
  message("\tNumber of independent samples : ", 
          length(samples.yi[!samples.yi %in% samples.y]), "\n")
})
