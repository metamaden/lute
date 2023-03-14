library(lute)
# source("./R/lute_utilities.R")

#--------------------
# final class example
#--------------------
# example params
batch.variable <- "SubjectName"
celltype.variable <- "cellType"

# get data
lexample <- lute:::.get_decon_example_data_bisque()
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

#------------------
# test scdc function
#------------------
result <- SCDC::SCDC_prop(bulk.eset = y.eset, 
                          sc.eset = sc.eset,
                          ct.varname = celltype.variable, 
                          sample = batch.variable,
                          iter.max = iter.max, 
                          nu = nu, epsilon = epsilon, 
                          truep = truep,
                          ct.cell.size = s, 
                          ct.sub = celltype.subset)