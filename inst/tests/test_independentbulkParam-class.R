#!/usr/bin/env R

library(lute)

# source scripts
#source("./R/independentbulkParam-class.R")
#source("./R/lute_utilities.R")

# try example
lexample <- lute:::.get_decon_example_data_bisque()

# get example objects
y.eset <- lexample[["y.eset"]]
z.eset <- lexample[["sc.eset"]]
# check sample ids
colnames(y.eset)
unique(z.eset$SubjectName)
# set objects
yi <- exprs(y.eset[,grepl("bulk", colnames(y.eset))])
y <- exprs(y.eset[,grepl("sample", colnames(y.eset))])
z <- exprs(z.eset[, grepl(paste0('sample', seq(2)), z.eset[["SubjectName"]])])

param <- independentbulkParam(y = y, z = z, yi = yi)
res <- deconvolution(param)
res
