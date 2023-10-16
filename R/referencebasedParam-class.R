#!/usr/bin/env R

### Author: Sean Maden

#' referencebasedParam-class
#'
#' Class and methods for managing reference-based deconvolution methods.
#' 
#' @include deconvolutionParam-class.R
#' 
#' @details This is a parent class to manage reference-based deconvolution 
#' algorithms. 
#' 
#' Child/sub-classes of this are distinguished by their use of
#' either an explicit or implied \code{z} signature matrix (i.e. Z[G,K] of
#' dimensions G markers by K cell types). These also have an implied cell size 
#' term for biases from systematic cell size differences. If no cell size 
#' transformation is intended, this is the equivalent of passing equal size 
#' scales, (e.g. a K-length vector of equal values). See 
#' `vignette(package="lute")` for details about experiment terms.
#' 
#' @examples 
#' lexample <- get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#'
#' @returns New object.
#' 
setClass("referencebasedParam", contains="deconvolutionParam", 
         slots=c(z = "matrix", s = "numeric"))

#' Make new object of class referencebasedParam
#'
#' Main constructor for class \linkS4class{referencebasedParam}.
#'
#' @param y Bulk mixed signals matrix of samples, which can be matched to 
#' single-cell samples.
#' @param z Signature matrix of cell type-specific signals. If not provided, can 
#' be computed from a provided ExpressionSet containing single-cell data.
#' @param s Cell size factor transformations of length equal to the K cell types 
#' to deconvolve.
#' @param return.info Whether to return metadata and original method outputs 
#' with predicted proportions.
#'
#' @examples 
#' lexample <- get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#'
#' @returns New object of class \linkS4class{referencebasedParam}.
#'
#' @details Takes standard inputs for reference-based deconvolution algorithms.
#'
#' @returns New object.
#' 
#' @export
referencebasedParam <- function(y, z, s, return.info = FALSE) {
  new("referencebasedParam", y = y, z = z, s = s, return.info = return.info)
}

#' Deconvolution generic behavior for object of class \linkS4class{referencebasedParam}
#' @param object An object of class \linkS4class{referencebasedParam} (see 
#' \code{?referencebasedParam}).
#' @details Method for behavior of deconvolution generic when called for object 
#' of class \linkS4class{referencebasedParam}.
#' @examples 
#' lexample <- get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#' @returns Method results.
#' @export
setMethod("deconvolution", "referencebasedParam", function(object) {
  ## get metadata
  input_s <- object[["s"]]; input_y <- object[["y"]]; input_z <- object[["z"]]
  
  ## cell types in z, s
  if(is(input_s, "NULL")){input_s <- rep(1, ncol(input_z))}
  unique.types <- try(colnames(object[["z"]]))
  condition.z.types <- is(unique.types, "NULL")|is(unique.types, "try-error")
  if(!condition.z.types){
    unique.types <- unique.types[order(unique.types)]
    input_z <- input_z[,order(colnames(input_z), unique.types)]
    condition.s.types <- is(names(input_s), "NULL")
    if(!condition.s.types){
      filter.s.types <- names(input_s) %in% unique.types
      input_s <- input_s[filter.s.types]
      input_s <- input_s[order(names(input_s), unique.types)]
    }
  }
  input_z <- .zstransform(input_z, input_s)
  ## matching markers in y and z
  markers.y <- rownames(input_y); markers.z <- rownames(input_z)
  if(!is(markers.y, "NULL") & !is(markers.z, "NULL")){
    ## markers.y <- rownames(y); markers.z <- rownames(z)
    unique.markers <- unique(c(markers.y, markers.z))
    overlapping.markers <- intersect(markers.y, markers.z)
    y.filter <- rownames(input_y) %in% overlapping.markers
    z.filter <- rownames(input_z) %in% overlapping.markers
    input_y <- input_y[y.filter,,drop=FALSE]; input_z <- input_z[z.filter,,drop=FALSE]
    input_y <- input_y[order(match(rownames(input_y), overlapping.markers)),]
    input_z <- input_z[order(match(rownames(input_z), overlapping.markers)),]
  } else{
    message("Warning, rownames not provided in both y and z. ",
            "Can't match marker labels.")
  }
  ## parse additional warnings
  if(is(markers.y, "NULL")){message("Warning, object 'y' has no marker labels (rownames)\n")}
  if(is(markers.z, "NULL")){message("Warning, object 'z' has no marker labels (rownames)\n")}
  ## get final metadata
  input_g <- nrow(input_z); input_j <- ncol(input_y); input_k <- ncol(input_z)
  metadata.list <- list(g = input_g, j = input_j, k = input_k, 
                        s = input_s, unique.types = unique.types,
                        markers.y = markers.y, marker.z = markers.z)
  ## return list
  return(
    list(y = as.matrix(input_y), z = as.matrix(input_z), 
         s = as.numeric(input_s), metadata = metadata.list)
    )
})

#' Show generic behavior for object of class referencebasedParam
#' @param object Object of class \linkS4class{referencebasedParam} (see 
#' \code{?referencebasedParam}).
#' @examples 
#' lexample <- get_decon_example_data()
#' referencebasedParam(y = lexample$y, z = lexample$z, s = lexample$s)
#' @returns Prints data summary messages to console.
#' @export
setMethod("show", "referencebasedParam", function(object) {
  ## get metadata
  input_s <- object[["s"]]; input_y <- object[["y"]]; input_z <- object[["z"]]
  unique.types <- try(colnames(object[["z"]]))
  markers.y <- rownames(input_y); markers.z <- rownames(input_z)
  unique.markers <- unique(c(markers.y, markers.z))
  overlapping.markers <- intersect(markers.y, markers.z)
  input_g <- nrow(input_z); input_j <- ncol(input_y); input_k <- ncol(input_z)
  lmd <- list(g = input_g, j = input_j, k = input_k, 
              s = input_s, unique.types = unique.types, 
              markers.y = markers.y, marker.z = markers.z)
  ## post console messages
  cat(paste0("class: ", class(object)[1], "\n\n"))
  cat("key deconvolution run info:\n")
  cat("\tmarker info:\n")
  cat("\tsignature markers (Gz): ", input_g, "\n")
  cat("\tunique marker labels (Gy | Gz): ", length(unique.markers), "\n")
  cat("\toverlapping marker labels (Gy & Gz): ", length(overlapping.markers), "\n\n")
  ## bulk samples
  cat("\tsamples info:\n")
  cat("\tnumber of bulk samples (J): ", ncol(object[["y"]]), "\n")
  cat("\tsample labels: ", paste0(colnames(input_y), collapse = "; "), "\n")
  cat("\n")
  ## cell size factors
  cat("\tcell size factor properties:\n")
  if(!is(input_s, "NULL")){
    for(type in names(input_s)){
      cat("\tscale factor for type ", type, ": ", input_s[type], "\n")}
    if(length(input_s) == ncol(input_z)){
      input_z <- .zstransform(input_z, input_s)}
  }; cat("\n")
  ## cell types
  cat("\ttypes info:\n")
  cat("\tnumber of types (K): ", ncol(object[["z"]]), "\n")
  if(!(is(unique.types, "NULL")|is(unique.types, "try-error"))){
    unique.types <- unique.types[order(unique.types)]
    cat("\tunique type labels: ", paste0(unique.types, collapse = ";"), "\n")
  } else{
    cat("\nWarning, object 'z' has no type labels (colnames)\n")
  }; cat("\n")
  ## parse additional warnings
  if(is(markers.y, "NULL")){cat("Warning, object 'y' has no marker labels (rownames)\n\n")}
  if(is(markers.z, "NULL")){cat("Warning, object 'z' has no marker labels (rownames)\n\n")}
})
