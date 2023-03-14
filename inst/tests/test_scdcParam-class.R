#!/usr/bin/env R

library(SCDC)
library(lute)


# source("./R/scdcParam-class.R")
# source("./R/lute_utilities.R")
#

# try example
lexample <- .get_decon_example_data_bisque()
param <- bisqueParam(s = lexample[["s"]], 
                     y = lexample[["y"]], 
                     z = lexample[["z"]])