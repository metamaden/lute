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
set <- set_from_sce(sce, type.variable = "typevar")
set[["donor"]] <- gsub(".*;", "", set[["type"]])
set[["typevar"]] <- gsub(";.*", "", set[["type"]])
rowData(set)$marker_type <- c(rep("type1", 10), rep("type2", 10))

# heatmap of type;donor assays
get_set_heatmap(set, type.variable = "typevar", group.variable = "donor",
                mtype.variable = "marker_type")

# heatmap of type assays
set2 <- set_from_set(set, typevar = "typevar")
get_set_heatmap(set2, type.variable = "type", mtype.variable = "marker_type")

