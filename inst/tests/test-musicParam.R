source("lute_generics.R")
source("lute_utilities.R")
source("deconParam-class.R")
source("musicParam-class.R")
# example
lexample <- .get_decon_example_data()
param <- musicParam(s = lexample[["s"]], 
                   y = lexample[["y"]], 
                   z = lexample[["z"]])

# return only predicted proportions
deconvolution(param)
# type1     type2 
# 0.6770833 0.3229167

# return full results
param@return.info <- T
names(deconvolution(param))
# [1] "predictions" "result.info" "metadata"