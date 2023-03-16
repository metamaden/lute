#!/usr/bin/env R

# Author: Sean Maden
#
# Defines the donor.data.frame class, for tables containing donor-level signals 
# by marker and type. 
#
#

setClass(
  "donor.data.frame",
  #slots = list(donordf = "data.frame", ndonor = "integer", 
  #             ntype = "integer", nmarker = "integer"),
  contains = "data.frame"
)

setMethod("show",
          "donor.data.frame",
          function(object) {
            df <- object
            # parse metadata
            nmarker <- ndonor <- ntype <- 0
            typev <- "NA"
            cnv <- colnames(df)
            filt.donor <- grepl("^donor\\d", cnv)
            ndonor <- as.integer(length(which(filt.donor)))
            if(length(which(filt.donor))==0){df$donor1 <- NA}
            if("marker" %in% cnv){
              nmarker <- as.integer(length(unique(df[,"marker"])))
            }
            if("type" %in% cnv){
              ntype <- as.integer(length(unique(df[,"type"])))
              typev <- unique(df[,"type"])
            }
            cat("donor.data.frame head:\n")
            print(head(df))
            cat("ndonor:", ndonor, "\n")
            cat("ntype:", ntype, "\n")
            cat("nmarker:", nmarker, "\n")
            cat("unique_types: ", 
                paste0(unique(typev), collapse = ";"), 
                "\n")
          }
)
