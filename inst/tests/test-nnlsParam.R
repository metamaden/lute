require(lute)
# example
lexample <- lute:::.get_decon_example_data()
param <- nnlsParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
# return only predicted proportions
result <- deconvolution(param)
# return detailed results
param[["return.info"]] <- T
result <- deconvolution(param)

#------------------------
# test embedded algorithm
#------------------------
require(nnls)
# lparam <- callNextMethod()
lparam <- referencebasedParam(y = lexample[["y"]],
                              z = lexample[["z"]],
                              s = lexample[["s"]])
y <- lparam[["y"]]
z <- lparam[["z"]]
s <- lparam[["s"]]
bulk.samples.index.vector <- seq(ncol(y))
result <- lapply(
  bulk.samples.index.vector, function(index){
    
    nnls::nnls(A = z, b = y[,index])
    
  }
)

names(result) <- colnames(y)
predictions <- lapply(result, function(iter){iter$x})
lr <- .parse_deconvolution_predictions_results(predictions, 
                                               colnames(z), 
                                               colnames(y))



