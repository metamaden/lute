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
set <- set_from_sce(sce, group.variable = "donor", type.variable = "typevar")
metadata(set)[["set_plots"]]$heatmap

get_set_pca(set, assayname = "summarized_counts", type.variable = "type")

