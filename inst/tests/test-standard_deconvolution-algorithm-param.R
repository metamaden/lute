#!/usr/bin/env R

# Author: Sean Maden
#
# Standard algorithm test. This is for testing a new algorithm param class object.
#

require(lute)

lparam_from_lexample <- function(lexample){
  
  lparam <- referencebasedParam(
    
    y = lexample[["y"]], z = lexample[["z"]], s = lexample[["s"]]
    
    )
  
  return(lparam)
}

lparam <- lparam_from_lexample(lute:::.get_decon_example_data())

#------------
# default use
#------------
# run new param class with any applicable generics

#------------------------
# test embedded algorithm
#------------------------
# test algorithm contained by the new param class