require(lute)
# get data
lexample <- lute:::.get_decon_example_data_bisque()
sc.data <- lexample[["sc.eset"]]
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

# these evaluate to false -- lute does not currently run with independent bulk samples 
z <- lute:::.get_z_from_sce(sc.eset, "counts", "celltype")
sce <- random_sce()
lute(y.se = SummarizedExperiment(y.eset), 
     sce = SummarizedExperiment(sc.eset), 
     typemarker.algorithm = NULL,
     celltype.variable = "cellType")
lute(y.se = SummarizedExperiment(y.eset), 
     sce = sce, 
     typemarker.algorithm = NULL,
     celltype.variable = "cellType")



