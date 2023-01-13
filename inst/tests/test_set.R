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
set <- set_from_sce(sce, groupvar = "donor", typevar = "typevar")
set2 <- set_from_set(set, typevar = "celltype", groupvar = "donor")

# test heatmap
set[["donor"]] <- gsub(".*;", "", set[["type"]])
set[["typevar"]] <- gsub(";.*", "", set[["type"]])

get_set_heatmap(set)