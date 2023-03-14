library(lute)
# source("./R/lute_utilities.R")

#-----------------
# get example data
#-----------------
lexample <- lute:::.get_decon_example_data_bisque()

#-----------------------------------
# test bisqueParam with example data
#-----------------------------------
param <- bisqueParam(y.eset = lexample[["y.eset"]],
                     sc.eset = lexample[["sc.eset"]])
res <- deconvolution(param)

#-------------------------------------
# test example data on original method
#-------------------------------------
results <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset = lexample[["y.eset"]], 
                                                  sc.eset = lexample[["sc.eset"]])
# proportions <- results$bulk.props
# proportions.sample1 <- proportions$bulk.props[,1]
