#!/usr/bin/env R

# Author: Sean Maden
#
# Get diagram of deconvolution classes from lute.
#

libv <- c("DiagrammeR", "lute")
sapply(libv, library, character.only = T)

# load data table
csv.name <- "lute-deconvolution_transfer-learning-table.csv"
path <- system.file("csv", package = "lute")
path <- file.path(path, csv.name)
csv <- read.csv(path)

#-------------------
# get chart elements
#-------------------
# parse csv data
# parse parent node edges
input <- c("deconParam", "deconParam", 
           "referencebasedParam", 
           "referencebasedParam")
output <- c("referencebasedParam",
            "referencefreeParam",
            "independentbulkParam",
            "matchedbulkParam")
# parse final child class edges
paramv <- c("independentbulkParam", "matchedbulkParam")
for(param in paramv){
  filtv <- grepl(param, csv$parent_classes)
  methodv <- csv[filtv,]$method_shortname
  input <- c(input, rep(param, length(methodv)))
  output <- c(output, paste0(tolower(methodv), "Param"))
  nodev <- unique(c(input, output))
}

# get chart strings
# get node strings
node.string <- paste0("node [shape = box, ",
                      "fontname = Helvectica]\n",
                      paste0(nodev, collapse = ";"),
                      "\n")
# get edge strings
edge.string <- paste0(input, "->", output, collapse = " ")
# get chart strings
chart.string <- paste0("digraph new_flowchart {\n",
                       "graph [overl = true, fonsize = 10]",
                       node.string, "\n", 
                       edge.string, "\n",
                       "}")

#-----------
# make chart
#-----------
grViz(chart.string)

