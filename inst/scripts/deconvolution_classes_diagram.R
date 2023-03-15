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

#-----------
# parameters
#-----------
append.string <- "Param"
method.colname <- "method_shortname"
parent.class.colname <- "parent_classes"

#-----------------
# helper functions
#-----------------
add_method_edges <- function(filter.string, csv, input, output,
                             filter.type = "in", end.append.name = "Param",
                             input.node.string = "referencebasedParam",
                             filter.column = "parent_classes",
                             method.column = "method_class"){
  # adds edges for methods, where originating node is called node.string
  method.vector <- csv[,method.column]
  filter <- grepl(filter.string, csv[,filter.column])
  if(!filter.type == "in"){filter <- !filter}
  new.edge.out <- unique(method.vector[filter])
  # get formatted end node names/edges
  # new.edge.out <- tolower(unique(method.vector[filter]))
  # new.edge.out <- paste0(new.edge.out, end.append.name) # append string
  new.edge.in <- rep(input.node.string, length(new.edge.out))
  return(list(input = c(input, new.edge.in), 
              output = c(output, new.edge.out)))
}

#----------------
# get chart edges
#----------------
# parse csv data
# initialize parent class nodes and edges
input <- c("deconParam", "deconParam", "referencebasedParam")
output <- c("referencebasedParam", "referencefreeParam", "independentbulkParam")

# parse final child class edges
# append independent bulk methods
lchart <- add_method_edges(filter.string = "independentbulkParam",
                           filter.type = "in", csv = csv,
                           input = input, output = output, 
                           input.node.string = "independentbulkParam")
# append remaining methods
lchart <- add_method_edges(filter.string = "independentbulkParam", 
                           filter.type = "out", csv = csv,
                           lchart[["input"]], lchart[["output"]], 
                           input.node.string = "referencebasedParam")

# get edge strings
edge.string <- paste0(lchart[["input"]], "->", 
                      lchart[["output"]], collapse = " ")

#----------------
# get chart nodes
#----------------
# get chart strings
methods <- paste0(unique(csv[,method.colname]), append.string)
node.vector <- unique(c(lchart[["input"]], lchart[["output"]]))
# get node strings
filter <- node.vector %in% methods
# set nodes for parent classes
nodev.filter <- node.vector[!filter]
node.string.parent <- paste0("node [shape = box, color = 'grey',",
                      "fontname = Helvectica]\n",
                      paste0(nodev.filter, collapse = ";"),
                      "\n")
# set nodes for methods
nodev.filter <- node.vector[filter]
node.string.method <- paste0("node [shape = cricle, fillcolor = grey, ",
                      "fontname = Helvectica]\n",
                      paste0(nodev.filter, collapse = ";"),
                      "\n")
# final node.string
node.string <- paste0("\n\n",node.string.parent,"\n\n", node.string.method)

#----------------------------------
# make chart from component strings
#----------------------------------
# get chart strings
chart.string <- paste0("digraph new_flowchart {\n",
                       "graph [overl = true, fonsize = 10]",
                       node.string, "\n", 
                       edge.string, "\n",
                       "}")

grViz(chart.string)
