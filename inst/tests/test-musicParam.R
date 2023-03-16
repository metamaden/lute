require(lute)
# example
lexample <- lute:::.get_decon_example_data()
param <- musicParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
# return only predicted proportions
deconvolution(param)
# return full results
param@return.info <- T
deconvolution(param)
