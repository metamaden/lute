#!/usr/bin/env R

# Author: Sean Maden
#
# Defines methods for object class to contain deconvolution predictions, 
# metadata, and benchmarking data.
#

setClass(
  "deconvolution.results",
  contains = "list"
)

setMethod("show",
          "deconvolution.results",
          function(object) {
            ldecon <- object
            lmd <- ldecon[["metadata"]]
            ntypes <- nmarkers <- method <- NA
            cat("object of type 'deconvolution.results'\nobject summary:\n")
            cat("method: ", lmd[["method"]], "\n")
            cat("number of markers:", lmd[["number.of.markers"]], "\n")
            cat("number of types:", lmd[["number.of.types"]], "\n")
            cat("type labels:", paste0(lmd[["types"]], collapse = "; "), "\n")
          }
)