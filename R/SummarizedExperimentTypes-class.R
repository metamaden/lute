#!/usr/bin/env R

# Define SummarizedExperimentSets class, validity tests, maker and conversion 
# functions, etc. Borrows heavily from definitions of established Bioconductor 
# classes for inspiration, including from packages minfi and 
# SummarizedExperiment.
#


setClass(
  "SummarizedExperimentTypes",
  representation(type_summary = "character", annotation = "character"),
  contains = "SummarizedExperiment"
)

SummarizedExperimentTypes <- function(assays, 
                                      type_summary,
                                      gr = GRanges(), 
                                      annotation = "", 
                                      ...) {
  #assays <- Assays(SimpleList(assays), as.null.if.no.assay=TRUE)
  new("SummarizedExperimentTypes",
      SummarizedExperiment(assays = assays, ...), 
      annotation = annotation
  )
}

SummarizedExperimentTypesRanges <- function(assays, 
                                      type_summary,
                                      gr = GRanges(), 
                                      annotation = "", 
                                      ...) {
  #assays <- Assays(SimpleList(assays), as.null.if.no.assay=TRUE)
  new("SummarizedExperimentTypes",
      SummarizedExperiment(assays = assays, 
                           rowRanges = as(gr, "GRanges"),
                           ...), 
      annotation = annotation
  )
}

#---------------
# define methods
#---------------
# setMethod("show", 
#          signature(object = "SummarizedExperimentTypes"), 
#          function(object) {
#            callNextMethod()
#            .show.annotation(annotation(object))
#            .show.type_summary(type_summary(object))
#            }
#          )









