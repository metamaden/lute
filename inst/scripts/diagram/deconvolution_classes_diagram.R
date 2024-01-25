#!/usr/bin/env R

### Author: Sean Maden
###
### Get diagram of deconvolution classes from lute.
###
### 1. A. deconvolutionParam
### 2. B. referencebasedParam
### 3. C. independentbulkParam
### 4. 1. nnlsParam
### 5. 2. musicParam
### 6. 3. epicParam
### 7. 4. deconrnaseqParam
### 8. 5. music2Param
### 9. 6. music2Param
### 10. 7. bisqueParam
### 11. 8. scdcParam
###
###
###


libv <- c("DiagrammeR", "lute")
sapply(libv, library, character.only=T)

## load data table
csv.name <- "lute-deconvolution_transfer-learning-table.csv"
path <- system.file("csv", package="lute")
path <- file.path(path, csv.name)
csv <- read.csv(path)

###-----------
### parameters
###-----------
append.string <- "Param"
method.colname <- "method_class"
parent.class.colname <- "parent_classes"

###-----------------
### helper functions
###-----------------
add_method_edges <- function(filter.string, csv, input, output,
                             filter.type="in", end.append.name="Param",
                             input.node.string="referencebasedParam",
                             filter.column="parent_classes",
                             method.column="method_class"){
  ## adds edges for methods, where originating node is called node.string
  method.vector <- csv[,method.column]
  filter <- grepl(filter.string, csv[,filter.column])
  if(!filter.type == "in"){filter <- !filter}
  new.edge.out <- unique(method.vector[filter])
  ## get formatted end node names/edges
  ## new.edge.out <- tolower(unique(method.vector[filter]))
  ## new.edge.out <- paste0(new.edge.out, end.append.name) # append string
  new.edge.in <- rep(input.node.string, length(new.edge.out))
  return(list(input=c(input, new.edge.in), 
              output=c(output, new.edge.out)))
}

###----------------
### get chart edges
###----------------
## parse csv data
## initialize parent class nodes and edges

#input <- c("deconParam", "deconParam", "referencebasedParam")
#output <- c("referencebasedParam", "referencefreeParam", "independentbulkParam")

input <- c("deconParam", "deconParam", "referencebasedParam")
output <- c("referencebasedParam", "referencefreeParam", "independentbulkParam")

## parse final child class edges
## append independent bulk methods
lchart <- add_method_edges(filter.string="independentbulkParam",
                           filter.type="in", csv=csv,
                           input=input, output=output, 
                           input.node.string="independentbulkParam")
## append remaining methods
lchart <- add_method_edges(filter.string="independentbulkParam", 
                           filter.type="out", csv=csv,
                           lchart[["input"]], lchart[["output"]], 
                           input.node.string="referencebasedParam")

## get edge strings
edge.string <- paste0(lchart[["input"]], "->", 
                      lchart[["output"]], collapse=" ")

###----------------
### get chart nodes
###----------------
## get chart strings
node.vector <- unique(c(lchart[["input"]], lchart[["output"]]))
filter <- node.vector %in% csv[,method.colname]
## get node strings
## set nodes for parent classes
nodev.filter <- node.vector[!filter]
node.string.parent <- paste0("node [shape=box, color='black',",
                      "fontname=Courier]\n",
                      paste0(nodev.filter, collapse=";"),
                      "\n")
## set nodes for methods
nodev.filter <- node.vector[filter]
node.string.method <- paste0("node [shape=oval, color='black', fillcolor='lightgray', ",
                      "fontname=Courier, style=filled]\n",
                      paste0(nodev.filter, collapse=";"),
                      "\n")
## final node.string
node.string <- paste0("\n\n",node.string.parent,"\n\n", node.string.method)

###----------------------------------
### make chart from component strings
#----------------------------------
## get chart strings
chart.string <- paste0("digraph new_flowchart {\n",
                       "graph [overlap=false, fonsize=10]",
                       node.string, "\n", 
                       edge.string, "\n",
                       "}")

grViz(chart.string)

message(chart.string)


###-----------------------
### define chart functions
###-----------------------
add_method_edges <- function(filter.string, csv, input, output,
                             filter.type="in", end.append.name="Param",
                             input.node.string="referencebasedParam",
                             filter.column="parent_classes",
                             method.column="method_class"){
  ## adds edges for methods, where originating node is called node.string
  method.vector <- csv[,method.column]
  filter <- grepl(filter.string, csv[,filter.column])
  if(!filter.type == "in"){filter <- !filter}
  new.edge.out <- unique(method.vector[filter])
  new.edge.in <- rep(input.node.string, length(new.edge.out))
  return(list(input=c(input, new.edge.in), 
              output=c(output, new.edge.out)))
}

deconvolution_classes_chart <- function(csv.filename=paste0("lute-deconvolution",
                                                              "_transfer-learning-table.csv"),
                                        method.colname="method_class",
                                        parent.class.colname="parent_classes",
                                        #input.node.start=c("deconvolutionParam", 
                                        #                     "deconvolutionParam", 
                                        #                     "referencebasedParam"),
                                        #output.node.start=c("referencebasedParam", 
                                        #                      "referencefreeParam", 
                                        #                      "independentbulkParam"),
                                        input.node.start=c("A", "A", "B"),
                                        output.node.start=c("B", "B", "C"),
                                        method.node.fillcolor="lightgray",
                                        method.node.outlinecolor="black",
                                        method.node.shape="oval",
                                        parent.node.fillcolor="white",
                                        parent.node.outlinecolor="black",
                                        parent.node.shape="box",
                                        method.node.font="Courier",
                                        parent.node.font="Courier"){
  require(DiagrammeR)
  ## load csv
  csv.name <- "lute-deconvolution_transfer-learning-table.csv"
  path <- file.path(system.file("csv", package="lute"), csv.name)
  csv <- read.csv(path)
  
  ## set chart edges
  input <- input.node.start; output <- output.node.start
  ## append independent bulk methods
  #lchart <- add_method_edges(filter.string="independentbulkParam",
  #                           filter.type="in", csv=csv,
  #                           input=input, output=output, 
  #                           input.node.string="independentbulkParam")
  lchart <- add_method_edges(filter.string="C",
                             filter.type="in", csv=csv,
                             input=input, output=output, 
                             input.node.string="C")
  ## append remaining methods
  #lchart <- add_method_edges(filter.string="independentbulkParam", 
  #                           filter.type="out", csv=csv,
  #                           lchart[["input"]], lchart[["output"]], 
  #                           input.node.string="referencebasedParam")
  lchart <- add_method_edges(filter.string="C", 
                             filter.type="out", csv=csv,
                             lchart[["input"]], lchart[["output"]], 
                             input.node.string="B")
  ## set edge strings
  edge.string <- paste0(lchart[["input"]], "->", lchart[["output"]], 
                        collapse=" ")
  
  ## set node labels and formats
  node.vector <- unique(c(lchart[["input"]], lchart[["output"]]))
  filter <- node.vector %in% csv[,method.colname]
  ## set nodes for parent classes
  nodev.filter <- node.vector[!filter]
  node.string.parent <- paste0("node [shape=", parent.node.shape, ", ",
                               "color='", parent.node.outlinecolor, "',",
                               "fillcolor='", parent.node.fillcolor, "',",
                               "fontname=", parent.node.font, ",",
                               "style=filled]\n",
                               paste0(nodev.filter, collapse=";"), "\n")
  ## set nodes for methods
  nodev.filter <- node.vector[filter]
  node.string.method <- paste0("node [shape=", method.node.shape, ", ",
                               "color='", method.node.outlinecolor,"', ",
                               "fillcolor='", method.node.fillcolor,"', ",
                               "fontname=", method.node.font, ", ",
                               "style=filled]\n",
                               paste0(nodev.filter, collapse=";"), "\n")
  ## get final node.string
  node.string <- paste0("\n\n",node.string.parent,"\n\n", node.string.method)
  
  ## make chart from component strings
  ## get chart strings
  chart.string <- paste0("digraph new_flowchart {\n",
                         "graph [overlap=false, fonsize=10]",
                         node.string, "\n", edge.string, "\n", "}")
  chart <- grViz(chart.string)
  return(list(chart=chart, chart.string=chart.string))
}

## example
new.chart <- deconvolution_classes_chart()
new.chart$chart

##-----
## save
##-----
jpeg("deconvolutionParam_hierarchy_diagram.jpg",
     res = 400, units = "in", width = 4, height = 4)
new.chart$chart
dev.off()
