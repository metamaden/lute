#' independentbulkParam-class
#'
#' Class and methods for managing methods requiring independent bulk samples.
#'
#' @details The main purpose of this class is to compare bulk sample data between the passed objects y and yi.
#' Since we assume yi contains the independent bulk samples, it should not have overlapping sample IDs (colnames),
#' and it should have overlapping marker IDs (rownames) compared to the reference bulk samples y.
#' 
#' @include referencebasedParam-class.R
#'
#' @param yi Mixed signals matrix from bulk samples, independent from primary mixed signals matrix y.
#' 
#' @examples 
#' lexample <- .get_decon_example_data()
#' 
setClass("independentbulkParam", contains="referencebasedParam", slots = c(yi = "matrix"))

#' Function to get nnlsParam
#' @export
independentbulkParam <- function(y = NULL, yi = NULL, z = NULL, s = NULL, return.info = FALSE) {
    if(is(y, "NULL"))}{y <- matrix(0)}
    if(is(z, "NULL"){z <- matrix(0)})
    if(is(yi, "NULL"){yi <- matrix(0)})
    if(is(s, "NULL"){s <- 0})
  new("independentbulkParam", y = y, z = z, s = s, return.info = return.info)
}

#' @export
setMethod("deconvolution", "independentbulkParam", function(object) {
  # get bulk data
  y <- object[["y"]]; yi <- object[["yi"]]

  # compare bulk markers
  markers.y <- rownames(y)
  markers.yi <- rownames(yi)
  # compare bulk samples
  samples.y <- colnames(y)
  samples.yi <- colnames(yi)

  #  

  # return list
  return(list(y = y, z = z, s = s, metadata = lmd))
})

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
