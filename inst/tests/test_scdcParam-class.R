library(lute)
# source("./R/lute_utilities.R")

#--------------------
# final class example
#--------------------
# get data
lexample <- lute:::.get_decon_example_data_bisque()
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
res <- deconvolution(param)

# get full results
param@return.info <- TRUE
res <- deconvolution(param)

#------------------
# load example data
#------------------
# example params
batch.variable <- "SubjectName"
celltype.variable <- "cellType"

# get data
lexample <- lute:::.get_decon_example_data_bisque()
sc.eset <- lexample[["sc.eset"]]
y.eset <- lexample[["y.eset"]]

#-----------------------------------
# test bisqueParam with example data
#-----------------------------------
param <- scdcParam(y.eset = y.eset, sc.eset = sc.eset,
                     batch.variable = "SubjectName",
                     celltype.variable = "cellType")
res <- deconvolution(param)