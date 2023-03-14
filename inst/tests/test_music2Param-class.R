library(lute)
# source("./R/lute_utilities.R")

#--------------------
# final class example
#--------------------
# example params
batch.variable <- "SubjectName"
celltype.variable <- "cellType"

# get data
lexample <- lute:::.get_decon_example_data_scdc()
sc.eset <- lexample[["sc.eset"]]
y.eset <- lexample[["y.eset"]]

# get yi
id.sc <- unique(sc.eset[[batch.variable]])
id.bulk <- colnames(y.eset)
filter <- id.bulk[!id.bulk %in% id.sc]
yi <- exprs(y.eset)[,filter]

# get param object
param <- scdcParam(yi = yi, y.eset = y.eset, sc.eset = sc.eset,
                   batch.variable = "SubjectName",
                   celltype.variable = "cellType")

# get just predictions
res <- deconvolution(param)

# get full results
param@return.info <- TRUE
res <- deconvolution(param)

#---------------------
# test music2 function
#---------------------

result <- MuSiC::music2_prop()
