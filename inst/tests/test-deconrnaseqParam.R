#!/usr/bin/env R

# Author: Sean Maden
#
# Tests for the DeconRNASeq deconvolution method parameters class, 
# deconrnaseqParam.
#


source("./R/lute_generics.R")
source("./R/lute_utilities.R")
source("./R/deconParam-class.R")
source("./R/referencebasedParam-class.R")
source("./R/deconrnaseqParam-class.R")

# example
lexample <- .get_decon_example_data()
param <- deconrnaseqParam(s = lexample[["s"]], 
                          y = lexample[["y"]], 
                          z = lexample[["z"]])

# return only predicted proportions
deconvolution(param)
# type1     type2 
# 0.9819837 0.0180163 

# return full results
param@return.info <- T
names(deconvolution(param))
# [1] "predictions" "result.info" "metadata"