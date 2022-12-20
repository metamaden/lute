#!/usr/bin/env R

# Author: Sean Maden
#
# Define set class, validity tests, maker and conversion functions, etc.
#


setClass(
  "SummarizedExperimentSets",
  representation(typeSummary = "character", annotation = "character"),
  contains = "SummarizedExperiment"
)

SummarizedExperimentTypes <- function(assays, 
                                      type_summary,
                                      gr = GRanges(), 
                                      annotation = "", ...) {
  assays <- Assays(assays, as.null.if.no.assay=TRUE)
  new("SummarizedExperimentTypes",
      SummarizedExperiment(
        assays = assays,
        rowRanges = as(gr, "GRanges"),
        ...),
      annotation = annotation
  )
}

#---------------------
# define methods
#---------------------

setMethod("show", signature(object = "GenomicMethylSet"), function(object) {
  callNextMethod()
  .show.annotation(annotation(object))
  .show.preprocessMethod(preprocessMethod(object))
})