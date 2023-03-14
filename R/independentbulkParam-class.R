#' independentbulkParam-class
#'
#' Class and methods for managing methods requiring independent bulk samples.
#'
#' @include lute_generics.R
#' @include referencebasedParam-class.R
#'
#' @param yi Mixed signals matrix from bulk samples, independent from primary mixed signals matrix y.
#'
#' @details The main purpose of this class is to compare bulk sample data between the passed objects y and yi.
#' Since we assume yi contains the independent bulk samples, it should not have overlapping sample IDs (colnames),
#' and it should have overlapping marker IDs (rownames) compared to the reference bulk samples y.
#'
#' @seealso \linkS4class{deconParam}, \linkS4class{referencebasedParam}
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
    if(is(s, "NULL")){s <- rep(1, ncol(z))}
  param <- new("independentbulkParam", y = y, yi = yi, z = z, s = s, return.info = return.info)
  # check inputs with parse
  param <- deconvolution(param)
  return(param)
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
    cat("Warning, no marker labels found in y.\n")
    } else if(is(markers.yi, "NULL")){
        cat("Warning, no marker labels found in yi.\n")
    } else{
        cat("Found marker labels in y, yi. Comparing...\n")
        unique.markers <- unique(markers.y, markers.yi)
        overlapping.markers <- intersect(markers.y, markers.yi)
        cat("Found ", length(overlapping.markers), " overlapping markers.\n")
        if(length(overlapping.markers) > 0){
            cat("Subsetting yi on markers overlapping in y.\n")
            yi <- yi[overlapping.markers,]
        }
    }

  # parse bulk sample IDs
  samples.y <- colnames(y)
  samples.yi <- colnames(yi)
  # compare sample IDs
  if(is(samples.y, "NULL")){
    cat("Warning, no sample labels found in y.\n")
    } else if(is(samples.yi, "NULL")){
        cat("Warning, no sample labels found in yi.\n")
    } else{
        cat("Found sample labels in y, yi. Comparing...\n")
        unique.samples <- unique(samples.y, samples.yi)
        overlapping.samples <- intersect(samples.y, samples.yi)
        cat("Found ", length(overlapping.markers), " overlapping samples.\n")
        if(length(overlapping.samples) > 0){
            cat("Removing overlapping samples from yi.\n")
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
    if("metadata" %in% names(object)){
      lmd <- object[["metadata"]]
      cat("Data summaries for class `independentbulkParam`:\n")
      if("unique.samples" %in% names(lmd)){
        cat("\tUnique sample IDs : ", lmd[["unique.samples"]], "\n")  
      } else{
        cat("Warning, didn't find object 'unique.samples'.\n")
      }
      if("unique.markers" %in% names(lmd)){
        cat("\tUnique marker IDs : ", lmd[["unique.markers"]], "\n")  
      } else{
        cat("Warning, didn't find object 'unique.markers'.\n")
      }
      if("yi" %in% names(lmd)){
        cat("\tTotal independent samples : ", ncol(object[["yi"]]), "\n")
        } else{
            cat("Warning, didn't find object `yi`.\n")
        } 
    }
})
