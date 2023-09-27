require(lute)
# example
lexample <- lute:::.get_decon_example_data()
param <- epicParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
# return only predicted proportions
results <- deconvolution(param)
# return full results
param@return.info <- T
results <- deconvolution(param)
