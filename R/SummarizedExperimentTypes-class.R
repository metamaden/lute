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

RangedSummarizedExperimentTypes <- function(assays, 
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
setMethod("show",
          "SummarizedExperimentTypes",
          function(sce) {
            sce <- object
            nmarker <- nrow(sce)
            ntype <- ncol(sce)
            mdnames <- names(metadata(sce))
            assaynames <- names(assays(sce))
            cat("nmarkers: ", nmarker, "\n")
            cat("ntypes: ", ntype, "\n")
            cat("assaynames: ", 
                paste0(assaynames, collapse = ";"), "\n")
            cat("mdnames: ", 
                paste0(mdnames, collapse = ";"), "\n")
          }
)







