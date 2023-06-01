#!/usr/bin/env R

# Author: Sean Maden

#' referencebasedParam-class
#'
#' Class and methods for managing reference-based deconvolution methods.
#' 
#' @include deconvolutionParam-class.R
#' 
#' @details This is a parent class to manage reference-based deconvolution 
#' algorithms. Child/sub-classes of this are distinguished by their use of
#' either an explicit or implied \code{z} signature matrix (i.e. Z[G,K] of
#' dimensions G markers by K cell types). These also have an implied cell size
#' transformation term, which is used for biases arising from systematic 
#' cell size differences. If no cell size transformation is intended, this is
#' the equivalent of passing equal size scales, (e.g. a K-length vector of equal 
#' values).
#' 
#' @examples 
#' lexample <- lute:::.get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#' 
setClass("referencebasedParam", contains="deconvolutionParam", 
         slots=c(z = "matrix", s = "numeric"))

#' Make new object of class referencebasedParam
#'
#' Main constructor for class \linkS4class{referencebasedParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types to deconvolve.
#' @param return.info Whether to return metadata and original method outputs with predicted proportions.
#'
#' @examples 
#' lexample <- lute:::.get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#'
#' @returns New object of class \linkS4class{referencebasedParam}.
#'
#' @details Takes standard inputs for reference-based deconvolution algorithms.
#' 
#' @export
referencebasedParam <- function(y, z, s, return.info = FALSE) {
  new("referencebasedParam", y = y, z = z, s = s, return.info = return.info)
}

#' Deconvolution generic behavior for object of class \linkS4class{referencebasedParam}
#' @param object An object of class \linkS4class{referencebasedParam}.
#' @details Method for behavior of deconvolution generic when called for object of class 
#' \linkS4class{referencebasedParam}
#' @examples 
#' lexample <- lute:::.get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#' @export
setMethod("deconvolution", "referencebasedParam", function(object) {
  # get metadata
  s <- object[["s"]]; y <- object[["y"]]; z <- object[["z"]]
  
  # cell types in z, s
  if(is(s, "NULL")){s <- rep(1, ncol(z))}
  unique.types <- try(colnames(object[["z"]]))
  condition.z.types <- is(unique.types, "NULL")|is(unique.types, "try-error")
  if(!condition.z.types){
    unique.types <- unique.types[order(unique.types)]
    z <- z[,order(colnames(z), unique.types)]
    condition.s.types <- is(names(s), "NULL")
    if(!condition.s.types){
      filter.s.types <- names(s) %in% unique.types
      s <- s[filter.s.types]
      s <- s[order(names(s), unique.types)]
    }
  }
  z <- .zstransform(z, s)
  
  # matching markers in y and z
  markers.y <- rownames(y)
  markers.z <- rownames(z)
  if(!is(markers.y, "NULL") & !is(markers.z, "NULL")){
    # markers.y <- rownames(y); markers.z <- rownames(z)
    unique.markers <- unique(c(markers.y, markers.z))
    overlapping.markers <- intersect(markers.y, markers.z)
    y.filter <- rownames(y) %in% overlapping.markers
    z.filter <- rownames(z) %in% overlapping.markers
    y <- y[y.filter,,drop=F]; z <- z[z.filter,,drop=F]
    y <- y[order(match(rownames(y), overlapping.markers)),]
    z <- z[order(match(rownames(z), overlapping.markers)),]
  } else{
    message("Warning, rownames not provided in both y and z. ",
            "Can't match marker labels.")
  }
  
  # parse additional warnings
  if(is(markers.y, "NULL")){cat("Warning, object 'y' has no marker labels (rownames)\n")}
  if(is(markers.z, "NULL")){cat("Warning, object 'z' has no marker labels (rownames)\n")}
  
  # get final metadata
  g <- nrow(z); j <- ncol(y); k <- ncol(z)
  lmd <- list(g = g, j = j, k = k, s = s, unique.types = unique.types, 
              markers.y = markers.y, marker.z = markers.z)
  # return list
  return(list(y = y, z = z, s = s, metadata = lmd))
})

#' Show generic behavior for object of class referencebasedParam
#' @param object Object of class \linkS4class{referencebasedParam}.
#' @examples 
#' lexample <- lute:::.get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#' @returns Prints data summary messages to console.
#' @export
setMethod("show", "referencebasedParam", function(object) {
  # get metadata
  s <- object[["s"]]; y <- object[["y"]]; z <- object[["z"]]
  unique.types <- try(colnames(object[["z"]]))
  markers.y <- rownames(y); markers.z <- rownames(z)
  unique.markers <- unique(c(markers.y, markers.z))
  overlapping.markers <- intersect(markers.y, markers.z)
  g <- nrow(z); j <- ncol(y); k <- ncol(z)
  lmd <- list(g = g, j = j, k = k, s = s, unique.types = unique.types, 
              markers.y = markers.y, marker.z = markers.z)
  # post console messages
  cat(paste0("class: ", class(object)[1], "\n\n"))
  cat("key deconvolution run info:\n")
  cat("\tmarker info:\n")
  cat("\tsignature markers (Gz): ", g, "\n")
  cat("\tunique marker labels (Gy | Gz): ", length(unique.markers), "\n")
  cat("\toverlapping marker labels (Gy & Gz): ", length(overlapping.markers), "\n\n")
  # bulk samples
  cat("\tsamples info:\n")
  cat("\tnumber of bulk samples (J): ", ncol(object[["y"]]), "\n")
  cat("\tsample labels: ", paste0(colnames(y), collapse = "; "), "\n")
  cat("\n")
  # cell size factors
  cat("\tcell size factor properties:\n")
  s <- object[["s"]]
  if(!is(s, "NULL")){
    for(type in names(object[["s"]])){
      cat("\tscale factor for type ", type, ": ", s[type], "\n")}
    if(length(s) == ncol(z)){z <- .zstransform(z, s)}
  }; cat("\n")
  # cell types
  cat("\ttypes info:\n")
  cat("\tnumber of types (K): ", ncol(object[["z"]]), "\n")
  if(!(is(unique.types, "NULL")|is(unique.types, "try-error"))){
    unique.types <- unique.types[order(unique.types)]
    cat("\tunique type labels: ", paste0(unique.types, collapse = ";"), "\n")
  } else{
    cat("\nWarning, object 'z' has no type labels (colnames)\n")
  }; cat("\n")
  # parse additional warnings
  if(is(markers.y, "NULL")){cat("Warning, object 'y' has no marker labels (rownames)\n\n")}
  if(is(markers.z, "NULL")){cat("Warning, object 'z' has no marker labels (rownames)\n\n")}
})
