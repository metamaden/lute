source("lute_generics.R")
source("lute_utilities.R")
source("deconParam-class.R")
source("deconrnaseqParam-class.R")
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