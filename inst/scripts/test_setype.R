#!/usr/bin/env R

# Author: Sean Maden
#
# Test methods for SummarizedExperimentTypes class
#
#

libv <- c("SingleCellExperiment", "lute")
sapply(libv, library, character.only = TRUE)

#----------------
# make sce object
#----------------
sce <- random_sce()
# assign new group var "donor"
# two donors type1, one donor type2:
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
# make set from sce
set <- set_from_sce(sce, typevar = "celltype", method = "mean")
