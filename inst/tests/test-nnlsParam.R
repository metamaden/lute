require(lute)
# example
lexample <- lute:::.get_decon_example_data()
param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
# return only predicted proportions
result <- deconvolution(param)
# return detailed results
param[["return.info"]] <- T
result <- deconvolution(param)