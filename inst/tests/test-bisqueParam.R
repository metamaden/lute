require(lute)
# get data
lexample <- .get_decon_example_data_bisque()
sc.eset <- lexample[["sc.eset"]]
y.eset <- lexample[["y.eset"]]
# example params
batch.variable <- "SubjectName"
celltype.variable <- "cellType"
# get param object
param <- bisqueParam(y.eset = y.eset, sc.eset = sc.eset, 
                     batch.variable = "SubjectName",
                     celltype.variable = "cellType")
# get just predictions
results <- deconvolution(param)
# get full results
param@return.info <- TRUE
results <- deconvolution(param)
