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
    if(is(y, "NULL")){y <- matrix(0)}
    if(is(z, "NULL")){z <- matrix(0)}
    if(is(yi, "NULL")){yi <- matrix(0)}
    if(is(s, "NULL")){s <- 0}
  new("independentbulkParam", y = y, z = z, s = s, return.info = return.info)
}

#' @export
setMethod("deconvolution", "independentbulkParam", function(object) {
  # get bulk data
  y <- object[["y"]]; yi <- object[["yi"]]

  # parse bulk marker IDs
  markers.y <- rownames(y)
  markers.yi <- rownames(yi)
  # compare markers
  if(is(markers.y, "NULL")){
    cat("Warning, no marker labels found in y.")
    } else if(is(markers.yi, "NULL")){
        cat("Warning, no marker labels found in yi.")
    } else{
        cat("Found marker labels in y, yi. Comparing...")
        unique.markers <- unique(markers.y, markers.yi)
        overlapping.markers <- intersect(markers.y, markers.yi)
        cat("Found ", length(overlapping.markers), " overlapping markers.")
        if(length(overlapping.markers) > 0){
            cat("Subsetting yi on markers overlapping in y.")
            yi <- yi[overlapping.markers,]
        }
    }

  # parse bulk sample IDs
  samples.y <- colnames(y)
  samples.yi <- colnames(yi)
  # compare sample IDs
  if(is(samples.y, "NULL")){
    cat("Warning, no sample labels found in y.")
    } else if(is(samples.yi, "NULL")){
        cat("Warning, no sample labels found in yi.")
    } else{
        cat("Found sample labels in y, yi. Comparing...")
        unique.samples <- unique(samples.y, samples.yi)
        overlapping.samples <- intersect(samples.y, samples.yi)
        cat("Found ", length(overlapping.markers), " overlapping samples.")
        if(length(overlapping.samples) > 0){
            cat("Removing overlapping samples from yi.")
            yi <- yi[,!overlapping.samples,drop=F]
        }
    }

    # get return data
    lmd <- list(unique.markers = unique.markers,
                unique.samples = unique.samples,
                overlapping.markers = overlapping.markers,
                overlapping.samples = overlapping.samples)
  # return list
  return(list(y = y, yi = yi, metadata = lmd))
})

#' @export
setMethod("show", "independentbulkParam", function(object) {
  lmd <- object[["metadata"]]
  cat("Data summaries for class `independentbulkParam`:")
  cat("\tUnique sample IDs : ", lmd[["unique.samples"]])
  cat("\tUnique marker IDs : ", lmd[["unique.markers"]])
  cat("\tTotal independent samples : ", ncol(object[["yi"]]))
})
