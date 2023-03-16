require(lute)

# example
lexample <- lute:::.get_decon_example_data()
param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])

# return only predicted proportions
deconvolution(param)
# type1      type2 
# 0.48908543 0.05896868
# return full results

param@return.info <- T
names(deconvolution(param))
# [1] "predictions" "result.info" "metadata"