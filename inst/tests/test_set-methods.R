#!/usr/bin/env R

# Author: Sean Maden
#
# Testing methods for SummarizedExperimentTypes objects
#

require(lute)

#-------------
# set_from_sce
#-------------

sce = random_sce()
colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))

group.variable = "donor"
method = "mean"
type.variable = "celltype"
assayname = "counts"
make.set.plots = TRUE
verbose = FALSE


typev <- unique(sce[[type.variable]])

expr.set <- do.call(cbind, lapply(typev, function(typei){
  if(verbose){message("Summarizing type: ", typei, "...")}
  type.filt <- sce[[type.variable]]==typei
  scef <- sce[,type.filt]
  exprf <- assays(scef)[[assayname]]
  ma <- matrix(rowMeans(exprf), ncol = 1)
  colnames(ma) <- typei
  return(ma)
}))



