#!/usr/bin/env R

# Author: Sean Maden
#
#
#

require(lute)

# test random data
sce <- random_sce()
sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
sce[["typevar"]] <- paste0(sce[["celltype"]], ";", sce[["donor"]])
set <- set_from_sce(sce, group.variable = "donor", type.variable = "typevar")
set2 <- set_from_set(set, typevar = "celltype", groupvar = "donor")

# test heatmap
# get randomized singlecellexperiment
sce <- random_sce()
sce[["donor"]] <- c(rep("donor1", 2), rep("donor2", 8))
sce[["typevar"]] <- paste0(sce[["celltype"]], ";", sce[["donor"]])
# make new set from sce
set1 <- set_from_sce(sce, type.variable = "typevar")
set1[["donor"]] <- gsub(".*;", "", set1[["type"]])
set1[["typevar"]] <- gsub(";.*", "", set1[["type"]])
rowData(set1)$marker_type <- c(rep("type1", 10), rep("type2", 10))
# heatmap of type;donor assays
hm1 <- get_set_heatmap(set1, type.variable = "typevar", 
                       group.variable = "donor", mtype.variable = "marker_type")
# heatmap of type assays
set2 <- set_from_set(set1, typevar = "typevar")
hm2 <- get_set_heatmap(set2, type.variable = "type", 
                       mtype.variable = "marker_type")

