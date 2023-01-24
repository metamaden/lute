#!/usr/bin/env R

# Author: Sean Maden
#
# Test methods for SingleCellExperiments
#

#---------------------
# test sce_groupstat()
#---------------------
scef = random_sce()
group.variable
ugroupv
assayname = "counts"
summarytype = "rowData"
groupstat = c("numzero", "var")
verbose = FALSE

colData(sce)$donor <- c(rep("donor1", 7), rep("donor2", 3))
