#!/usr/bin/env R

# Author: Sean Maden
#
# Standard algorithm test. This is for testing a new algorithm param class object.
#

#------------
# default use
#------------
# run new param class with any applicable generics



require(lute)
# example
lexample <- lute:::.get_decon_example_data()
param <- musicParam(s = lexample[["s"]], y = lexample[["y"]], z = lexample[["z"]])
# return only predicted proportions
deconvolution(param)

# return full results
param@return.info <- T
deconvolution(param)






#------------------------
# test embedded algorithm
#------------------------
# test algorithm contained by the new param class



require(MuSiC); require(dplyr)
lparam <- callNextMethod()
# instantiate objects
nu <- object[["nu"]]
iter.max <- object[["iter.max"]]
eps <- object[["eps"]]
sigma <- object[["sigma"]]
y <- lparam[["y"]]
z <- lparam[["z"]]
s <- lparam[["s"]]
sigma <- as.matrix(sigma)
bulk.samples.index.vector <- seq(ncol(y))
result <- lapply(bulk.samples.index.vector, function(index){
  MuSiC::music.basic(X = z, 
                     Y = y[,index,drop=F], 
                     S = s, 
                     Sigma = sigma, 
                     nu = nu, 
                     iter.max = iter.max, 
                     eps = eps)
})
names(result) <- colnames(y)
predictions <- lapply(result, function(iter){iter$p.weight})
predictions <- do.call(rbind, predictions)
predictions <- apply(predictions, 1, function(ri){ri/sum(ri)}) %>% t()
colnames(predictions) <- colnames(z)
rownames(predictions) <- colnames(y)
lr <- predictions


