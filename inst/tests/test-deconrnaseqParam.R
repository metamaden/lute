require(lute)
# example
lexample <- lute:::.get_decon_example_data()
param <- deconrnaseqParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
# return only predicted proportions
results <- deconvolution(param)
# return full results
param@return.info <- T
results <- deconvolution(param)
