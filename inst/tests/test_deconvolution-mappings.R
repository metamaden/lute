#!/usr/bin/env R

# Test mappings to deconvolution functions in lute.

require(lute)

sce <- random_sce()
typev <- unique(sce[["celltype"]])
Z <- do.call(cbind, lapply(typev, function(typei){
  rowMeans(counts(sce[,sce[["celltype"]]==typei]))
}))
Y <- matrix(rowMeans(counts(sce)), ncol = 1)

# run nnls
ldecon <- run_deconvolution(method = "nnls", Y = Y, Z = Z)

## run music
# ldecon <- run_deconvolution(method = "music", Y = Y, Z = Z)

# inspect results
names(ldecon)

# try music
method <- "music"
arguments <- list("sigma" = NULL)
sigma <- matrix(0, ncol = 1, nrow = nrow(Z))

arguments <- list("sigma = NULL")

command.list <- map_music(arguments)

