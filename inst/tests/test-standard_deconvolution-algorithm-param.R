#!/usr/bin/env R

# Author: Sean Maden
#
# Standard algorithm test. This is for testing a new algorithm param class object.
#

require(lute)
# require(argparse)

lparam_from_lexample <- function(lexample){
  
  lparam <- referencebasedParam(
    y = lexample[["y"]], z = lexample[["z"]], s = lexample[["s"]]
  )
  
  return(lparam)
}

lparam <- lparam_from_lexample(lute:::.get_decon_example_data())

output.from.one.bulk.sample <- function(lexample,
                                        algo.index.from.mappings = "nnls", 
                                        file.name = "lute-deconvolution_transfer-learning-table.csv"){
  mappings.data <- read.csv(
    file.path(
      system.file(package = "lute"), "csv", file.name))
  mappings.algo <- mappings.data[mappings.data[,"method_shortname"]==algo.index.from.mappings,]
  function.string <- mappings.algo["function."] %>% as.character()
  lute_y = mappings.algo["lute_y"] %>% as.character()
  lute_z = mappings.algo["lute_z"] %>% as.character()
  eval(
    parse(
      text = 
        paste0(
          function.string, "(", lute_z, " = lexample[['z']], ", lute_y, " = lexample[['y']][,1])")))
}

#------------
# default use
#------------
# run new param class with any applicable generics

#------------------------
# test embedded algorithm
#------------------------
# test algorithm contained by the new param class